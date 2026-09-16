#!/usr/bin/env python3
"""Targeted Ag3.0 site-level external validation from public Sanger all-sites VCFs.

Scientific role
---------------
EXTERNAL VALIDATION ONLY. The target set must already have been locked by
R/05b_freeze_external_validation_targets.R from Phase-2 discovery outputs.
Ag3 never selects targets, changes discovery weights, or reranks genes.

Bandwidth design
----------------
MalariaGEN publishes one ~3-GB *all-sites* bgzip/tabix VCF per Ag3 specimen at
https://vo_agam_output.cog.sanger.ac.uk/.  This script does NOT download those
whole files.  HTSlib/pysam performs indexed HTTPS range requests for only the
small genomic windows covering the locked 23-bp targets.  A small taxon-balanced
pilot is mandatory before --final can run.

Fail-closed rules
-----------------
* HTTP byte-range support is checked before genomic access.
* Every target must have exactly 23 genomic records in an all-sites VCF before
  it can be called for that sample.
* Coordinate drift, malformed ploidy, excessive sample failures, or loss of
  the immutable target-lock checksum abort the run rather than being ignored.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, NoReturn, Sequence

import numpy as np
import pandas as pd
import requests
import yaml
try:
    import pysam
except ImportError as exc:
    raise SystemExit(
        "Missing Python dependency 'pysam'. Run scripts/setup_hybrid_python.ps1 first."
    ) from exc

SCHEMA_VERSION = "ag3-external-public-all-sites-vcf-v2"
CANONICAL = ("2R", "2L", "3R", "3L", "X")
FIXED_TAXA = ("gambiae", "coluzzii")


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def log(msg: str) -> None:
    print(f"{now()}\t{msg}", flush=True)


def fail(msg: str, code: int = 2) -> NoReturn:
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(code)


def locate_root(explicit: str | None) -> Path:
    root = Path(explicit).expanduser().resolve() if explicit else Path(__file__).resolve().parents[1]
    if not (root / "analysis_config.yml").exists():
        fail(f"Project root does not contain analysis_config.yml: {root}")
    return root


def load_cfg(root: Path) -> dict[str, Any]:
    with (root / "analysis_config.yml").open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent)); os.close(fd)
    p = Path(tmp)
    try:
        p.write_text(text, encoding="utf-8", newline="")
        os.replace(p, path)
    finally:
        p.unlink(missing_ok=True)


def atomic_df(df: pd.DataFrame, path: Path, sep: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent)); os.close(fd)
    p = Path(tmp)
    try:
        df.to_csv(p, sep=sep, index=False)
        os.replace(p, path)
    finally:
        p.unlink(missing_ok=True)


def atomic_json(path: Path, obj: Any) -> None:
    atomic_text(path, json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")


def get_lines(url: str, timeout: int = 90) -> list[str]:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return [x.strip() for x in r.text.splitlines() if x.strip()]


def normalize_taxon(x: str) -> str:
    z = str(x).strip().lower().replace("an. ", "").replace("anopheles ", "")
    return {
        "gambiae": "gambiae", "gambiae s.s.": "gambiae", "gambiae_ss": "gambiae",
        "coluzzii": "coluzzii", "m": "coluzzii", "s": "gambiae",
    }.get(z, z)


def sample_id_from_feature_url(url: str) -> str:
    b = url.rstrip("/").split("/")[-1]
    for suf in (".gatk.zarr.zip", ".zarr.zip", ".vcf.gz"):
        if b.endswith(suf):
            return b[:-len(suf)]
    return b.split(".")[0]


def load_public_manifest(features_url: str, labels_url: str, vcf_base_url: str) -> pd.DataFrame:
    # The public classifier manifests are used only as a stable sample-ID/taxon
    # catalogue. Genomic calls are read from MalariaGEN's official all-sites VCFs.
    urls = get_lines(features_url)
    labels = get_lines(labels_url)
    if len(urls) != len(labels):
        fail(f"Ag3 public feature/label manifests have different lengths: {len(urls)} vs {len(labels)}")
    if len(urls) < 1000:
        fail(f"Ag3 public feature manifest unexpectedly contains only {len(urls)} rows.")
    d = pd.DataFrame({"feature_url": urls, "taxon_raw": labels})
    d["sample_id"] = d.feature_url.map(sample_id_from_feature_url)
    d["taxon"] = d.taxon_raw.map(normalize_taxon)
    if d.sample_id.duplicated().any():
        fail(f"Duplicate sample IDs in Ag3 public manifest: {d.loc[d.sample_id.duplicated(),'sample_id'].head().tolist()}")
    base = vcf_base_url.rstrip("/")
    d["vcf_url"] = d.sample_id.map(lambda s: f"{base}/{s}.vcf.gz")
    return d


def deterministic_select(md: pd.DataFrame, taxa: Sequence[str], per_taxon: int) -> pd.DataFrame:
    # Deterministic rather than adaptive selection prevents cherry-picking.
    chunks = []
    for taxon in taxa:
        d = md[md.taxon == taxon].sort_values("sample_id").copy()
        if d.empty:
            fail(f"No samples labelled {taxon!r} in public Ag3 manifest.")
        if per_taxon > 0:
            d = d.head(per_taxon)
        chunks.append(d)
    return pd.concat(chunks, ignore_index=True)


def head_range_support(url: str, timeout: int = 30) -> tuple[bool, int | None, str]:
    try:
        r = requests.head(url, allow_redirects=True, timeout=timeout)
        if r.status_code >= 400:
            return False, None, f"HEAD HTTP {r.status_code}"
        size_s = r.headers.get("Content-Length")
        size = int(size_s) if size_s and size_s.isdigit() else None
        ar = (r.headers.get("Accept-Ranges") or "").lower()
        if "bytes" in ar:
            return True, size, "Accept-Ranges: bytes"
        q = requests.get(url, headers={"Range": "bytes=0-0"}, stream=True, timeout=timeout)
        ok = q.status_code == 206 and any(k.lower() == "content-range" for k in q.headers)
        status = q.status_code
        q.close()
        return ok, size, f"range-probe HTTP {status}"
    except Exception as exc:
        return False, None, repr(exc)


def merge_windows(targets: pd.DataFrame, gap: int = 500, max_span: int = 25000) -> list[tuple[str,int,int]]:
    out: list[tuple[str,int,int]] = []
    for contig, d in targets.groupby("genomic_seqid", sort=False):
        iv = sorted(zip(d.genomic_start.astype(int), d.genomic_end.astype(int)))
        if not iv: continue
        s0,e0 = iv[0]
        for s,e in iv[1:]:
            ne = max(e0,e)
            if s <= e0 + gap + 1 and ne - s0 + 1 <= max_span:
                e0 = ne
            else:
                out.append((str(contig), int(s0), int(e0)))
                s0,e0 = s,e
        out.append((str(contig), int(s0), int(e0)))
    return out


def pam_positions(start: int, end: int, strand: str) -> set[int]:
    if strand == "+": return {end-2, end-1, end}
    if strand == "-": return {start, start+1, start+2}
    fail(f"Invalid guide strand: {strand}")


def genotype_state(rec) -> tuple[bool, bool, int, int]:
    """Return called?, reference-only?, called-allele count, nonref-allele count."""
    samples = list(rec.samples.values())
    if len(samples) != 1:
        fail(f"Expected one sample in per-sample VCF record, observed {len(samples)} at {rec.contig}:{rec.pos}")
    gt = samples[0].get("GT")
    if gt is None or len(gt) == 0:
        return False, False, 0, 0
    alleles = [a for a in gt if a is not None and int(a) >= 0]
    if not alleles:
        return False, False, 0, 0
    return True, all(int(a) == 0 for a in alleles), len(alleles), sum(int(a) != 0 for a in alleles)


def sample_target_states(vcf_url: str, sample_id: str, taxon: str, targets: pd.DataFrame, windows: list[tuple[str,int,int]]) -> list[dict[str,Any]]:
    # Collect only records returned by the small locked windows.
    by_pos: dict[tuple[str,int], Any] = {}
    try:
        vf = pysam.VariantFile(vcf_url)
    except Exception as exc:
        raise RuntimeError(f"Could not open remote tabix VCF: {exc}") from exc
    try:
        header_samples = list(vf.header.samples)
        if len(header_samples) != 1:
            raise RuntimeError(f"Expected one VCF sample; header has {len(header_samples)}: {header_samples[:5]}")
        for contig, start, end in windows:
            try:
                it = vf.fetch(contig, start-1, end)
            except Exception as exc:
                raise RuntimeError(f"tabix fetch failed for {contig}:{start}-{end}: {exc}") from exc
            for rec in it:
                p = int(rec.pos)
                if start <= p <= end:
                    by_pos[(contig,p)] = rec
        rows: list[dict[str,Any]] = []
        for _, t in targets.iterrows():
            contig = str(t.genomic_seqid); st = int(t.genomic_start); en = int(t.genomic_end)
            ppam = pam_positions(st,en,str(t.guide_genomic_strand))
            records = [by_pos.get((contig,p)) for p in range(st,en+1)]
            present = sum(r is not None for r in records)
            called_all = present == 23
            exact = True
            pam_called = True
            pam_intact = True
            called_alleles = 0
            nonref_alleles = 0
            missing_positions = []
            for p, rec in zip(range(st,en+1), records):
                if rec is None:
                    called_all = False; exact = False; missing_positions.append(p)
                    if p in ppam: pam_called = False; pam_intact = False
                    continue
                called, refonly, ncall, nnonref = genotype_state(rec)
                called_alleles += ncall; nonref_alleles += nnonref
                if not called:
                    called_all = False; exact = False
                    if p in ppam: pam_called = False; pam_intact = False
                elif not refonly:
                    exact = False
                    if p in ppam: pam_intact = False
            if not called_all: exact = False
            if not pam_called: pam_intact = False
            rows.append({
                "sample_id": sample_id, "taxon": taxon,
                "population_site_id": str(t.population_site_id), "contig": contig,
                "records_present_23bp": present,
                "callable_23bp": bool(called_all), "exact_23bp": bool(exact) if called_all else np.nan,
                "pam_callable": bool(pam_called), "pam_intact": bool(pam_intact) if pam_called else np.nan,
                "called_alleles_23bp": int(called_alleles), "nonref_alleles_23bp": int(nonref_alleles),
                "missing_positions": ";".join(map(str,missing_positions)),
            })
        return rows
    finally:
        vf.close()


def aggregate(states: pd.DataFrame, targets: pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame]:
    def one(d: pd.DataFrame) -> dict[str,Any]:
        called = d[d.callable_23bp == True]
        pamc = d[d.pam_callable == True]
        ca = int(d.called_alleles_23bp.sum()); na = int(d.nonref_alleles_23bp.sum())
        return {
            "n_samples_total": int(len(d)),
            "n_samples_called_23bp": int(len(called)),
            "n_samples_exact_23bp": int(called.exact_23bp.fillna(False).astype(bool).sum()),
            "n_samples_pam_called": int(len(pamc)),
            "n_samples_pam_intact": int(pamc.pam_intact.fillna(False).astype(bool).sum()),
            "callable_sample_fraction": float(len(called)/len(d)) if len(d) else np.nan,
            "exact_23bp_fraction": float(called.exact_23bp.astype(float).mean()) if len(called) else np.nan,
            "pam_intact_fraction": float(pamc.pam_intact.astype(float).mean()) if len(pamc) else np.nan,
            "called_alleles_23bp": ca, "nonreference_alleles_23bp": na,
            "nonreference_allele_fraction_23bp": float(na/ca) if ca else np.nan,
        }
    arows=[]; brows=[]
    for sid,d in states.groupby("population_site_id",sort=False):
        arows.append({"population_site_id":sid,**one(d)})
    for (sid,tax),d in states.groupby(["population_site_id","taxon"],sort=False):
        brows.append({"population_site_id":sid,"taxon":tax,**one(d)})
    ann=targets[["population_site_id","gene_id","selection_reason"]].drop_duplicates()
    return pd.DataFrame(arows).merge(ann,on="population_site_id",how="left"), pd.DataFrame(brows).merge(ann,on="population_site_id",how="left")


def main() -> int:
    ap=argparse.ArgumentParser(description="Locked Ag3 external validation via public all-sites tabix VCFs.")
    ap.add_argument("--project-root",default=None)
    m=ap.add_mutually_exclusive_group(); m.add_argument("--pilot",action="store_true"); m.add_argument("--final",action="store_true")
    ap.add_argument("--dry-run",action="store_true")
    ap.add_argument("--force",action="store_true")
    args=ap.parse_args()

    root=locate_root(args.project_root); os.chdir(root); cfg=load_cfg(root)
    ev=cfg.get("external_validation",{}) or {}
    if not bool(ev.get("enabled",True)): fail("external_validation.enabled is false")
    if bool(ev.get("require_no_discovery_retuning",True)) is not True: fail("Validation isolation must remain enabled.")
    lock_file=root/"data_processed"/"05b_external_validation_targets.csv"
    if not lock_file.exists(): fail("Missing 05b external target lock. Run Phase-2 discovery through R Step 05b first.")
    targets=pd.read_csv(lock_file)
    required={"population_site_id","gene_id","genomic_seqid","genomic_start","genomic_end","guide_genomic_strand","selection_locked_before_ag3"}
    miss=sorted(required-set(targets.columns))
    if miss: fail(f"External target lock missing columns: {miss}")
    if targets.population_site_id.duplicated().any(): fail("Duplicate site IDs in target lock.")
    if not targets.selection_locked_before_ag3.astype(bool).all(): fail("Target lock is not immutable for every row.")
    if not targets.genomic_seqid.astype(str).isin(CANONICAL).all(): fail("Noncanonical target in external lock.")
    if not ((targets.genomic_end.astype(int)-targets.genomic_start.astype(int)+1)==23).all(): fail("External targets must all be 23 bp.")
    lock_hash=sha256_file(lock_file)

    features_url=str(ev.get("features_url")); labels_url=str(ev.get("labels_url"))
    vcf_base=str(ev.get("public_vcf_base_url","https://vo_agam_output.cog.sanger.ac.uk"))
    md=load_public_manifest(features_url,labels_url,vcf_base)
    taxa=[str(x).lower() for x in ev.get("primary_taxa",list(FIXED_TAXA))]
    md=md[md.taxon.isin(taxa)].copy()
    if md.empty: fail("No configured primary taxa in Ag3 public manifest.")

    run_mode="final" if args.final else "pilot"
    if args.final:
        pilot_ok=root/str(ev.get("output_dir","data_raw/ag3_external"))/"pilot"/"ACQUISITION_COMPLETE.ok"
        if not pilot_ok.exists(): fail("--final is blocked until the small Ag3 pilot completes successfully.")
        per=int(ev.get("final_max_samples_per_taxon",250) or 0)
    else:
        per=int(ev.get("pilot_samples_per_taxon",25) or 25)
    selected=deterministic_select(md,taxa,per)
    out_dir=root/str(ev.get("output_dir","data_raw/ag3_external"))/run_mode
    out_dir.mkdir(parents=True,exist_ok=True)
    atomic_df(selected,out_dir/"selected_samples.tsv")
    atomic_df(targets,out_dir/"locked_targets.tsv")
    gap=int(ev.get("region_merge_gap_bp",500) or 500); maxspan=int(ev.get("max_region_span_bp",25000) or 25000)
    windows=merge_windows(targets,gap,maxspan)
    window_bp=sum(e-s+1 for _,s,e in windows)
    log(f"Project: {root}")
    log(f"Locked external targets: {len(targets):,} sites / {targets.gene_id.nunique():,} genes; SHA256={lock_hash[:16]}...")
    log(f"Ag3 {run_mode}: selected {len(selected):,} samples ({', '.join(f'{t}={int((selected.taxon==t).sum())}' for t in taxa)})")
    log(f"Per-sample tabix plan: {len(windows):,} merged small windows / {window_bp:,} logical bp")
    log("Bandwidth guard: public all-sites VCFs remain remote; only indexed HTTPS ranges for locked windows are read")
    log("Validation-isolation guard: Ag3 cannot change Phase-2 discovery scores/ranks/thresholds")
    if args.dry_run:
        log("DRY RUN complete: tiny public manifests read; no genomic VCF opened")
        return 0

    # Confirm byte-range serving for first 3 VCFs and tiny tabix indexes.
    for _,r in selected.head(min(3,len(selected))).iterrows():
        ok,size,detail=head_range_support(str(r.vcf_url))
        if not ok: fail(f"VCF server lacks reliable byte ranges; refusing possible whole download: {r.vcf_url} ({detail})")
        idx_url=str(r.vcf_url)+".tbi"; ir=requests.head(idx_url,allow_redirects=True,timeout=30)
        if ir.status_code>=400: fail(f"Missing tabix index for {r.sample_id}: HTTP {ir.status_code} {idx_url}")
        log(f"Range access confirmed: {r.sample_id}; remote VCF size={size if size is not None else 'unknown'} bytes")

    states_rows=[]; errors=[]; nsel=len(selected); max_err=float(ev.get("max_sample_failure_fraction",0.10) or 0.10)
    for j,(_,r) in enumerate(selected.iterrows(),start=1):
        sid,tax,url=str(r.sample_id),str(r.taxon),str(r.vcf_url)
        log(f"[{j}/{nsel}] Ag3 targeted all-sites VCF: {sid} ({tax})")
        try:
            if sha256_file(lock_file)!=lock_hash: fail("External target-lock checksum changed during acquisition; aborting.")
            states_rows.extend(sample_target_states(url,sid,tax,targets,windows))
        except SystemExit: raise
        except Exception as exc:
            errors.append({"sample_id":sid,"taxon":tax,"vcf_url":url,"error":repr(exc)})
            log(f"WARNING: sample failed: {sid}: {exc}")
            if len(errors)/j>max_err and j>=10:
                atomic_df(pd.DataFrame(errors),out_dir/"sample_errors.tsv")
                fail(f"Sample failure fraction exceeded {max_err:.0%}; stopping rather than silently reducing external cohort.")
    states=pd.DataFrame(states_rows); err=pd.DataFrame(errors,columns=["sample_id","taxon","vcf_url","error"])
    atomic_df(err,out_dir/"sample_errors.tsv")
    if states.empty: fail("No external sample-site states generated.")
    good=states.sample_id.nunique()
    if good<max(10,int(0.8*nsel)): fail(f"Only {good}/{nsel} selected samples succeeded; external cohort too incomplete.")
    expected=good*len(targets)
    if len(states)!=expected: fail(f"External state matrix incomplete: rows={len(states)}, expected={expected}.")
    site,tax=aggregate(states,targets)
    atomic_df(states,out_dir/"sample_site_states.tsv")
    atomic_df(site,out_dir/"site_summary.tsv")
    atomic_df(tax,out_dir/"site_taxon_summary.tsv")
    manifest={
        "schema_version":SCHEMA_VERSION,"run_mode":run_mode,
        "external_dataset":"Ag3.0 public Sanger per-sample all-sites VCF",
        "vcf_base_url":vcf_base,"feature_manifest_url":features_url,"label_manifest_url":labels_url,
        "target_lock_file":str(lock_file.relative_to(root)),"target_lock_sha256":lock_hash,
        "n_locked_targets":int(len(targets)),"n_target_windows":int(len(windows)),"target_window_logical_bp":int(window_bp),
        "n_selected_samples":int(nsel),"n_successful_samples":int(good),"n_failed_samples":int(len(err)),
        "taxa":taxa,"whole_vcf_files_downloaded":0,"access_method":"pysam/HTSlib tabix HTTPS byte ranges",
        "discovery_retuning_allowed":False,
    }
    atomic_json(out_dir/"run_manifest.json",manifest)
    atomic_text(out_dir/"ACQUISITION_COMPLETE.ok",SCHEMA_VERSION+"\n")
    log(f"Ag3 {run_mode} external validation acquisition completed: {good:,} samples x {len(targets):,} locked targets")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
