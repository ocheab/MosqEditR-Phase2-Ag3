#!/usr/bin/env python3
"""Targeted Ag3.0 site-level external validation from public Sanger per-sample Zarr.

Scientific role
---------------
This branch is EXTERNAL VALIDATION ONLY. It consumes the immutable target lock
created by R/05b_freeze_external_validation_targets.R. It never selects targets,
changes thresholds, or modifies Phase-2 discovery rankings.

Bandwidth design
----------------
* The 2,784-sample public URL/label manifests are tiny text files.
* Each sample remains remote as ``*.gatk.zarr.zip``.
* The script requires HTTP byte-range support and opens the ZIP/Zarr remotely.
* Only Zarr chunks covering the locked 23-bp target intervals are read.
* No per-sample whole-genome VCF/Zarr archive is intentionally downloaded.
* A small taxon-balanced pilot is the default. ``--final`` must be explicit.

Validation rule
---------------
A target is considered callable in a sample only when all 23 reference genomic
positions are present in the all-site coordinate axis and have at least one
called genotype allele. If the public Zarr is variant-only rather than all-site,
the script aborts instead of assuming missing positions are reference.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import fsspec
import numpy as np
import pandas as pd
import requests
import yaml
import zarr

SCHEMA_VERSION = "ag3-external-public-zarr-v1"
CANONICAL = ("2R", "2L", "3R", "3L", "X")
FIXED_TAXA = ("gambiae", "coluzzii")


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def log(msg: str) -> None:
    print(f"{now()}\t{msg}", flush=True)


def fail(msg: str, code: int = 2) -> "NoReturn":
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(code)


def locate_root(explicit: str | None) -> Path:
    if explicit:
        p = Path(explicit).expanduser().resolve()
        if not (p / "analysis_config.yml").exists():
            fail(f"Project root does not contain analysis_config.yml: {p}")
        return p
    p = Path(__file__).resolve().parent.parent
    if not (p / "analysis_config.yml").exists():
        fail(f"Could not identify project root from script location: {p}")
    return p


def load_cfg(root: Path) -> dict[str, Any]:
    with (root / "analysis_config.yml").open("r", encoding="utf-8") as f:
        x = yaml.safe_load(f)
    return x or {}


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


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def get_lines(url: str, timeout: int = 90) -> list[str]:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return [x.strip() for x in r.text.splitlines() if x.strip()]


def normalize_taxon(x: str) -> str:
    z = str(x).strip().lower().replace("an. ", "").replace("anopheles ", "")
    mapping = {
        "gambiae": "gambiae", "gambiae s.s.": "gambiae", "gambiae_ss": "gambiae",
        "coluzzii": "coluzzii", "m": "coluzzii", "s": "gambiae",
    }
    return mapping.get(z, z)


def sample_id_from_url(url: str) -> str:
    b = url.rstrip("/").split("/")[-1]
    for suf in (".gatk.zarr.zip", ".zarr.zip", ".vcf.gz"):
        if b.endswith(suf):
            return b[: -len(suf)]
    return b.split(".")[0]


def load_public_manifest(features_url: str, labels_url: str) -> pd.DataFrame:
    urls = get_lines(features_url)
    labels = get_lines(labels_url)
    if len(urls) != len(labels):
        fail(f"Ag3 public feature/label manifests have different lengths: {len(urls)} vs {len(labels)}")
    if len(urls) < 1000:
        fail(f"Ag3 public feature manifest unexpectedly contains only {len(urls)} rows.")
    out = pd.DataFrame({"url": urls, "taxon_raw": labels})
    out["sample_id"] = out.url.map(sample_id_from_url)
    out["taxon"] = out.taxon_raw.map(normalize_taxon)
    if out.sample_id.duplicated().any():
        d = out.loc[out.sample_id.duplicated(), "sample_id"].head().tolist()
        fail(f"Duplicate sample IDs in Ag3 public manifest: {d}")
    return out


def deterministic_select(md: pd.DataFrame, taxa: Sequence[str], per_taxon: int) -> pd.DataFrame:
    chunks = []
    for taxon in taxa:
        d = md[md.taxon == taxon].copy().sort_values("sample_id")
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
        ar = (r.headers.get("Accept-Ranges") or "").lower()
        size = r.headers.get("Content-Length")
        n = int(size) if size and size.isdigit() else None
        if "bytes" in ar:
            return True, n, "Accept-Ranges: bytes"
        # Some servers omit the header but honor Range; verify with a 1-byte GET.
        q = requests.get(url, headers={"Range": "bytes=0-0"}, stream=True, timeout=timeout)
        ok = q.status_code == 206 and "content-range" in {k.lower() for k in q.headers.keys()}
        q.close()
        return ok, n, f"range-probe HTTP {q.status_code}"
    except Exception as exc:
        return False, None, repr(exc)


class RemoteZipZarr:
    """Keep the HTTP handle and ZipFileSystem alive while Zarr arrays are used."""
    def __init__(self, url: str, block_size: int = 4 * 1024 * 1024):
        self.url = url
        self.block_size = block_size
        self.http = None
        self.zipfs = None
        self.group = None

    def __enter__(self):
        # HTTPFile uses ranged reads and a local readahead cache; no explicit full download.
        self.http = fsspec.open(
            self.url,
            mode="rb",
            block_size=self.block_size,
            cache_type="readahead",
        ).open()
        self.zipfs = fsspec.filesystem("zip", fo=self.http)
        mapper = self.zipfs.get_mapper("")
        try:
            self.group = zarr.open_consolidated(store=mapper, mode="r")
        except Exception:
            self.group = zarr.open_group(store=mapper, mode="r")
        return self.group

    def __exit__(self, exc_type, exc, tb):
        try:
            if self.zipfs is not None:
                self.zipfs.close()
        except Exception:
            pass
        try:
            if self.http is not None:
                self.http.close()
        except Exception:
            pass
        return False


def array_paths(group) -> list[str]:
    out: list[str] = []
    try:
        def visitor(name, obj):
            if isinstance(obj, zarr.Array):
                out.append(name)
        group.visititems(visitor)
    except Exception:
        pass
    return out


def find_contig_array(group, contig: str, kind: str):
    if kind == "pos":
        candidates = [f"{contig}/variants/POS", f"{contig}/variants/position", f"{contig}/POS"]
        suffixes = ("/variants/POS", "/variants/position", "/POS")
    elif kind == "gt":
        candidates = [f"{contig}/calldata/GT", f"{contig}/calldata/genotype", f"{contig}/GT"]
        suffixes = ("/calldata/GT", "/calldata/genotype", "/GT")
    else:
        raise ValueError(kind)
    for p in candidates:
        try:
            return group[p], p
        except Exception:
            pass
    paths = array_paths(group)
    hits = [p for p in paths if contig in p.split("/") and p.endswith(suffixes)]
    if len(hits) == 1:
        return group[hits[0]], hits[0]
    fail(f"Could not uniquely locate {kind} array for {contig}. Candidates={candidates}; example arrays={paths[:80]}")


def scalar_int(arr, i: int) -> int:
    return int(np.asarray(arr[i]).reshape(-1)[0])


def lower_bound(arr, x: int) -> int:
    lo, hi = 0, int(arr.shape[0])
    while lo < hi:
        mid = (lo + hi) // 2
        v = scalar_int(arr, mid)
        if v < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


def map_targets_to_rows(pos_arr, targets: pd.DataFrame, contig: str) -> tuple[np.ndarray, dict[str, list[int]], dict[int, int]]:
    """Find exact all-site rows for every locked 23-bp target on one contig."""
    all_rows: list[int] = []
    site_rows: dict[str, list[int]] = {}
    expected_by_row: dict[int, int] = {}
    for r in targets.itertuples(index=False):
        start, end = int(r.genomic_start), int(r.genomic_end)
        if end - start + 1 != 23:
            fail(f"External target {r.population_site_id} is not 23 bp: {contig}:{start}-{end}")
        lo = lower_bound(pos_arr, start)
        hi = lower_bound(pos_arr, end + 1)
        if hi <= lo:
            fail(
                f"Public Ag3 Zarr contains no coordinate records for {r.population_site_id} ({contig}:{start}-{end}). "
                "Site-level validation requires all-site coordinates; refusing to assume absent records are reference."
            )
        block = np.asarray(pos_arr[lo:hi], dtype=np.int64).reshape(-1)
        need = np.arange(start, end + 1, dtype=np.int64)
        if len(block) != 23 or not np.array_equal(block, need):
            missing = sorted(set(need.tolist()) - set(block.tolist()))
            fail(
                f"Public Ag3 Zarr is not all-site across locked target {r.population_site_id}: "
                f"found {len(block)}/23 positions; missing example={missing[:10]}. "
                "No false reference assumption was made."
            )
        rows = list(range(lo, hi))
        site_rows[str(r.population_site_id)] = rows
        for rr, pp in zip(rows, need.tolist()):
            expected_by_row[int(rr)] = int(pp)
        all_rows.extend(rows)
    unique_rows = np.asarray(sorted(set(all_rows)), dtype=np.int64)
    row_to_offset = {int(rr): i for i, rr in enumerate(unique_rows.tolist())}
    site_offsets = {sid: [row_to_offset[int(rr)] for rr in rows] for sid, rows in site_rows.items()}
    return unique_rows, site_offsets, expected_by_row


def slice_gt_rows(gt_arr, rows: np.ndarray) -> np.ndarray:
    if len(rows) == 0:
        return np.empty((0, 2), dtype=np.int16)
    ndim = len(gt_arr.shape)
    try:
        if ndim == 2:
            x = gt_arr.get_orthogonal_selection((rows, slice(None)))
        elif ndim == 3:
            # Per-sample stores normally have a singleton sample dimension.
            x = gt_arr.get_orthogonal_selection((rows, slice(None), slice(None)))
        else:
            fail(f"Unsupported GT array shape {gt_arr.shape}; expected 2-D or 3-D per-sample genotype array.")
    except Exception as exc:
        fail(f"Targeted GT row selection failed for shape {gt_arr.shape}: {exc}")
    x = np.asarray(x)
    if ndim == 3:
        if x.shape[1] != 1:
            fail(f"Per-sample Zarr GT array unexpectedly has sample dimension {x.shape}: expected 1.")
        x = x[:, 0, :]
    if x.ndim != 2 or x.shape[0] != len(rows):
        fail(f"Unexpected GT slice shape {x.shape} for {len(rows)} requested rows.")
    return x.astype(np.int16, copy=False)


def target_state(gt23: np.ndarray, strand: str) -> dict[str, Any]:
    if gt23.ndim != 2 or gt23.shape[0] != 23:
        raise ValueError(f"Expected 23 x ploidy genotype matrix, got {gt23.shape}")
    called = gt23 >= 0
    base_called = called.any(axis=1)
    callable23 = bool(base_called.all())
    nonref_called = (gt23 > 0) & called
    base_nonref = nonref_called.any(axis=1)
    exact23 = bool(callable23 and not base_nonref.any())
    if strand == "+":
        pam_idx = np.arange(20, 23)
    elif strand == "-":
        pam_idx = np.arange(0, 3)
    else:
        raise ValueError(f"Unexpected guide strand {strand!r}")
    pam_callable = bool(base_called[pam_idx].all())
    pam_intact = bool(pam_callable and not base_nonref[pam_idx].any())
    called_alleles = int(called.sum())
    nonref_alleles = int(nonref_called.sum())
    return {
        "callable_23bp": callable23,
        "exact_23bp": exact23 if callable23 else np.nan,
        "pam_callable": pam_callable,
        "pam_intact": pam_intact if pam_callable else np.nan,
        "called_alleles_23bp": called_alleles,
        "nonref_alleles_23bp": nonref_alleles,
    }


def aggregate(states: pd.DataFrame, targets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    tax_rows = []
    for keys, d in states.groupby(["population_site_id"], sort=False):
        sid = keys if isinstance(keys, str) else keys[0]
        called = d[d.callable_23bp == True]
        pamc = d[d.pam_callable == True]
        rows.append({
            "population_site_id": sid,
            "n_samples_total": int(len(d)),
            "n_samples_called_23bp": int(len(called)),
            "n_samples_exact_23bp": int(called.exact_23bp.fillna(False).astype(bool).sum()),
            "n_samples_pam_called": int(len(pamc)),
            "n_samples_pam_intact": int(pamc.pam_intact.fillna(False).astype(bool).sum()),
            "callable_sample_fraction": float(len(called)/len(d)) if len(d) else np.nan,
            "exact_23bp_fraction": float(called.exact_23bp.astype(float).mean()) if len(called) else np.nan,
            "pam_intact_fraction": float(pamc.pam_intact.astype(float).mean()) if len(pamc) else np.nan,
            "called_alleles_23bp": int(d.called_alleles_23bp.sum()),
            "nonreference_alleles_23bp": int(d.nonref_alleles_23bp.sum()),
            "nonreference_allele_fraction_23bp": float(d.nonref_alleles_23bp.sum()/d.called_alleles_23bp.sum()) if d.called_alleles_23bp.sum() else np.nan,
        })
    for (sid, tax), d in states.groupby(["population_site_id", "taxon"], sort=False):
        called = d[d.callable_23bp == True]
        pamc = d[d.pam_callable == True]
        tax_rows.append({
            "population_site_id": sid,
            "taxon": tax,
            "n_samples_total": int(len(d)),
            "n_samples_called_23bp": int(len(called)),
            "n_samples_exact_23bp": int(called.exact_23bp.fillna(False).astype(bool).sum()),
            "n_samples_pam_called": int(len(pamc)),
            "n_samples_pam_intact": int(pamc.pam_intact.fillna(False).astype(bool).sum()),
            "callable_sample_fraction": float(len(called)/len(d)) if len(d) else np.nan,
            "exact_23bp_fraction": float(called.exact_23bp.astype(float).mean()) if len(called) else np.nan,
            "pam_intact_fraction": float(pamc.pam_intact.astype(float).mean()) if len(pamc) else np.nan,
            "called_alleles_23bp": int(d.called_alleles_23bp.sum()),
            "nonreference_alleles_23bp": int(d.nonref_alleles_23bp.sum()),
            "nonreference_allele_fraction_23bp": float(d.nonref_alleles_23bp.sum()/d.called_alleles_23bp.sum()) if d.called_alleles_23bp.sum() else np.nan,
        })
    a = pd.DataFrame(rows)
    b = pd.DataFrame(tax_rows)
    ann = targets[["population_site_id", "gene_id", "selection_reason"]].drop_duplicates()
    return a.merge(ann, on="population_site_id", how="left"), b.merge(ann, on="population_site_id", how="left")


def main() -> int:
    ap = argparse.ArgumentParser(description="Targeted public Ag3 external validation without whole-genome downloads.")
    ap.add_argument("--project-root", default=None)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--pilot", action="store_true", help="Use configured small taxon-balanced pilot (default).")
    mode.add_argument("--final", action="store_true", help="Use the full configured external-validation cohort.")
    ap.add_argument("--dry-run", action="store_true", help="Read no genomic Zarr data; print target/sample plan only.")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    root = locate_root(args.project_root)
    os.chdir(root)
    cfg = load_cfg(root)
    ev = cfg.get("external_validation", {}) or {}
    if not bool(ev.get("enabled", True)):
        fail("external_validation.enabled is false in analysis_config.yml")
    if bool(ev.get("require_no_discovery_retuning", True)) is not True:
        fail("Validation-isolation guard requires external_validation.require_no_discovery_retuning: true")

    lock_file = root / "data_processed" / "05b_external_validation_targets.csv"
    if not lock_file.exists():
        fail("Missing data_processed/05b_external_validation_targets.csv. Run R/05b_freeze_external_validation_targets.R before ANY Ag3 validation access.")
    targets = pd.read_csv(lock_file)
    required = {"population_site_id","gene_id","genomic_seqid","genomic_start","genomic_end","guide_genomic_strand","selection_locked_before_ag3"}
    miss = sorted(required - set(targets.columns))
    if miss: fail(f"External target lock missing columns: {miss}")
    if targets.population_site_id.duplicated().any(): fail("Duplicate site IDs in external-validation target lock.")
    if not targets.selection_locked_before_ag3.astype(bool).all(): fail("Target-lock validation flag is not TRUE for every row.")
    if not targets.genomic_seqid.astype(str).isin(CANONICAL).all(): fail("External target lock contains a noncanonical contig.")
    if not ((targets.genomic_end.astype(int)-targets.genomic_start.astype(int)+1)==23).all(): fail("Every external-validation target must be exactly 23 bp.")

    features_url = str(ev.get("features_url"))
    labels_url = str(ev.get("labels_url"))
    md = load_public_manifest(features_url, labels_url)
    taxa = [str(x).lower() for x in ev.get("primary_taxa", list(FIXED_TAXA))]
    md = md[md.taxon.isin(taxa)].copy()
    if md.empty: fail("No configured primary taxa remain in Ag3 public manifest.")

    final = bool(args.final)
    if final:
        nmax = int(ev.get("final_max_samples_per_taxon", 0) or 0)
        selected = deterministic_select(md, taxa, nmax)
        run_mode = "final"
    else:
        npilot = int(ev.get("pilot_samples_per_taxon", 25) or 25)
        selected = deterministic_select(md, taxa, npilot)
        run_mode = "pilot"

    out_dir = root / str(ev.get("output_dir", "data_raw/ag3_external")) / run_mode
    out_dir.mkdir(parents=True, exist_ok=True)
    atomic_df(selected, out_dir / "selected_samples.tsv")
    atomic_df(targets, out_dir / "locked_targets.tsv")

    log(f"Project: {root}")
    log(f"External-validation lock: {len(targets):,} Phase-2-selected target sites across {targets.gene_id.nunique():,} genes")
    log(f"Ag3 public manifest: {len(md):,} configured primary-taxon samples; selected {len(selected):,} for {run_mode}")
    for t in taxa:
        log(f"  {t}: available={int((md.taxon==t).sum()):,}; selected={int((selected.taxon==t).sum()):,}")
    log("Validation-isolation guard: Ag3 data cannot modify discovery target selection, weights, classes, or ranks")
    log("Bandwidth guard: remote ZIP/Zarr byte-range reads only; no whole per-sample archive is intentionally downloaded")
    if args.dry_run:
        log("DRY RUN complete: small public manifests were read, but no genomic Zarr was opened")
        return 0

    # Require HTTP byte ranges on the first few selected archives before any Zarr access.
    for url in selected.url.head(min(3, len(selected))):
        ok, size, detail = head_range_support(url)
        if not ok:
            fail(f"Sanger object does not provide reliable HTTP byte-range access; refusing a possible full download: {url} ({detail})")
        log(f"Range access confirmed: {sample_id_from_url(url)}; remote archive size={size if size is not None else 'unknown'} bytes")

    # Build a coordinate-index template from the first selected sample. All other
    # samples are required to match those coordinates at every queried row.
    first = selected.iloc[0]
    index_template: dict[str, dict[str, Any]] = {}
    with RemoteZipZarr(first.url) as group:
        paths = array_paths(group)
        atomic_text(out_dir / "first_sample_array_paths.txt", "\n".join(paths) + "\n")
        for contig in CANONICAL:
            tt = targets[targets.genomic_seqid.astype(str) == contig]
            if tt.empty: continue
            pos_arr, pos_path = find_contig_array(group, contig, "pos")
            gt_arr, gt_path = find_contig_array(group, contig, "gt")
            rows, site_offsets, expected = map_targets_to_rows(pos_arr, tt, contig)
            # Verify GT has a compatible row axis before committing template.
            if int(gt_arr.shape[0]) <= int(rows.max(initial=-1)):
                fail(f"GT/POS row-axis mismatch on {contig}: POS={pos_arr.shape}, GT={gt_arr.shape}")
            index_template[contig] = {
                "rows": rows.tolist(),
                "site_offsets": site_offsets,
                "expected_positions": {str(k): int(v) for k,v in expected.items()},
                "pos_path": pos_path,
                "gt_path": gt_path,
                "pos_shape": list(pos_arr.shape),
                "gt_shape": list(gt_arr.shape),
            }
            log(f"{contig}: locked-coordinate template covers {len(site_offsets):,} targets / {len(rows):,} all-site rows")
    atomic_json(out_dir / "target_row_index_template.json", index_template)

    sample_rows: list[dict[str, Any]] = []
    error_rows: list[dict[str, str]] = []
    n_selected = len(selected)
    max_error_fraction = 0.10

    for i, sm in selected.iterrows():
        sid, tax, url = str(sm.sample_id), str(sm.taxon), str(sm.url)
        log(f"[{i+1}/{n_selected}] Ag3 targeted validation: {sid} ({tax})")
        try:
            with RemoteZipZarr(url) as group:
                for contig, templ in index_template.items():
                    tt = targets[targets.genomic_seqid.astype(str) == contig]
                    rows = np.asarray(templ["rows"], dtype=np.int64)
                    pos_arr, _ = find_contig_array(group, contig, "pos")
                    gt_arr, _ = find_contig_array(group, contig, "gt")
                    if list(pos_arr.shape) != templ["pos_shape"]:
                        fail(f"{sid} {contig}: POS shape differs from index-template sample ({pos_arr.shape} vs {templ['pos_shape']}).")
                    # Validate coordinates at queried rows; this is strict enough to
                    # prevent row-index drift without materializing the whole POS array.
                    observed = np.asarray(pos_arr.get_orthogonal_selection((rows,)), dtype=np.int64).reshape(-1)
                    expected = np.asarray([int(templ["expected_positions"][str(int(r))]) for r in rows], dtype=np.int64)
                    if not np.array_equal(observed, expected):
                        mismatch = np.where(observed != expected)[0][:5]
                        fail(f"{sid} {contig}: coordinate axis differs at locked rows; mismatch offsets={mismatch.tolist()}")
                    gt = slice_gt_rows(gt_arr, rows)
                    offsets = templ["site_offsets"]
                    tmap = tt.set_index("population_site_id")
                    for target_id, oo in offsets.items():
                        g23 = gt[np.asarray(oo, dtype=int), :]
                        st = target_state(g23, str(tmap.loc[target_id, "guide_genomic_strand"]))
                        sample_rows.append({
                            "sample_id": sid,
                            "taxon": tax,
                            "population_site_id": target_id,
                            "contig": contig,
                            **st,
                        })
        except SystemExit:
            raise
        except Exception as exc:
            error_rows.append({"sample_id": sid, "taxon": tax, "url": url, "error": repr(exc)})
            log(f"WARNING: sample failed: {sid}: {exc}")
            if len(error_rows) / max(1, i + 1) > max_error_fraction and (i + 1) >= 10:
                atomic_df(pd.DataFrame(error_rows), out_dir / "sample_errors.tsv")
                fail(f"Sample failure fraction exceeded {max_error_fraction:.0%}; stopping rather than silently reducing external cohort.")

    states = pd.DataFrame(sample_rows)
    errors = pd.DataFrame(error_rows, columns=["sample_id","taxon","url","error"])
    atomic_df(errors, out_dir / "sample_errors.tsv")
    if states.empty:
        fail("No external-validation sample-site states were generated.")

    good_samples = states.sample_id.nunique()
    if good_samples < max(10, int(0.8 * n_selected)):
        fail(f"Only {good_samples}/{n_selected} selected samples produced targeted states; external validation is too incomplete.")
    expected_rows = good_samples * len(targets)
    if len(states) != expected_rows:
        counts = states.groupby("sample_id").population_site_id.nunique()
        fail(f"External state matrix is incomplete: rows={len(states)}, expected={expected_rows}; per-sample site count range={counts.min()}-{counts.max()}")

    site_summary, taxon_summary = aggregate(states, targets)
    atomic_df(states, out_dir / "sample_site_states.tsv")
    atomic_df(site_summary, out_dir / "site_summary.tsv")
    atomic_df(taxon_summary, out_dir / "site_taxon_summary.tsv")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "run_mode": run_mode,
        "external_dataset": "Ag3.0 public Sanger per-sample GATK Zarr",
        "feature_manifest_url": features_url,
        "label_manifest_url": labels_url,
        "target_lock_file": str(lock_file.relative_to(root)),
        "target_lock_sha256": sha256_file(lock_file),
        "n_locked_targets": int(len(targets)),
        "n_selected_samples": int(n_selected),
        "n_successful_samples": int(good_samples),
        "n_failed_samples": int(len(errors)),
        "taxa": taxa,
        "whole_genome_archives_downloaded": 0,
        "access_method": "HTTP byte-range -> remote ZIP -> Zarr chunk reads",
        "discovery_retuning_allowed": False,
    }
    atomic_json(out_dir / "run_manifest.json", manifest)
    atomic_text(out_dir / "ACQUISITION_COMPLETE.ok", SCHEMA_VERSION + "\n")
    log(f"Ag3 {run_mode} external validation acquisition completed: {good_samples:,} samples x {len(targets):,} locked targets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
