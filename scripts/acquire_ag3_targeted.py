#!/usr/bin/env python3
"""
Targeted local acquisition for MosqEditR Manuscript 2.

This script intentionally DOES NOT download chromosome-scale VCF files.
It uses the official malariagen_data Ag3 API to fetch only:
  1) Ag3 release sample metadata;
  2) phased haplotype SNPs overlapping the frozen 23-bp CRISPR targets; and
  3) accessibility/site-mask values for small windows covering those targets.

Outputs are plain local CSV/TSV files consumed by R Steps 01-02.

Recommended Python: 3.10-3.12. The pinned malariagen_data release used by this
project does not support Python 3.13.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

try:
    import numpy as np
    import pandas as pd
    import yaml
    import malariagen_data
except ImportError as exc:
    raise SystemExit(
        "Missing Python dependency: %s\n\n"
        "Create the project Python environment first. On Windows PowerShell run:\n"
        "  powershell -ExecutionPolicy Bypass -File scripts/setup_ag3_python.ps1\n"
        "On Git Bash/WSL/Linux/macOS run:\n"
        "  bash scripts/setup_ag3_python.sh\n" % exc
    ) from exc

ACQUISITION_SCHEMA_VERSION = "2"
FIXED_GT_COLUMNS = ["variant_id", "contig", "position", "ref", "alt", "filter"]


@dataclass(frozen=True)
class Interval:
    contig: str
    start: int
    end: int

    @property
    def span(self) -> int:
        return self.end - self.start + 1

    @property
    def region(self) -> str:
        return f"{self.contig}:{self.start}-{self.end}"


def log(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    print(f"{stamp}\t{msg}", flush=True)


def fail(msg: str, code: int = 2) -> "None":
    raise SystemExit(f"ERROR: {msg}")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        tmp.write_text(text, encoding="utf-8", newline="")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)


def atomic_write_dataframe(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        df.to_csv(tmp, index=False, sep=sep)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)


def json_safe(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def atomic_write_json(path: Path, obj: dict) -> None:
    text = json.dumps(obj, indent=2, sort_keys=True, default=json_safe) + "\n"
    atomic_write_text(path, text)


def locate_root(explicit: str | None) -> Path:
    if explicit:
        root = Path(explicit).expanduser().resolve()
    else:
        root = Path(__file__).resolve().parents[1]
    if not (root / "analysis_config.yml").exists():
        fail(f"Project root does not contain analysis_config.yml: {root}")
    return root


def load_config(root: Path) -> dict:
    with (root / "analysis_config.yml").open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        fail("analysis_config.yml did not parse to a mapping.")
    return cfg


def cfg_get(cfg: dict, *keys, default=None):
    cur = cfg
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def normalize_bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    x = s.astype(str).str.strip().str.lower()
    return x.isin({"true", "t", "1", "yes", "y"})


def merge_intervals(
    rows: pd.DataFrame,
    gap_bp: int,
    max_span_bp: int,
) -> list[Interval]:
    out: list[Interval] = []
    if rows.empty:
        return out
    for contig, d in rows.groupby("genomic_seqid", sort=False):
        intervals = sorted(
            [(int(a), int(b)) for a, b in zip(d.genomic_start, d.genomic_end)],
            key=lambda z: (z[0], z[1]),
        )
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            proposed_e = max(cur_e, e)
            can_merge = s <= cur_e + gap_bp + 1
            within_span = proposed_e - cur_s + 1 <= max_span_bp
            if can_merge and within_span:
                cur_e = proposed_e
            else:
                out.append(Interval(str(contig), cur_s, cur_e))
                cur_s, cur_e = s, e
        out.append(Interval(str(contig), cur_s, cur_e))
    return out


def make_plan_df(intervals: Sequence[Interval], kind: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kind": kind,
            "window_id": [f"{kind}_{i:04d}" for i in range(1, len(intervals) + 1)],
            "contig": [x.contig for x in intervals],
            "start": [x.start for x in intervals],
            "end": [x.end for x in intervals],
            "span_bp": [x.span for x in intervals],
            "region": [x.region for x in intervals],
        }
    )


def batched(seq: Sequence, n: int) -> Iterable[Sequence]:
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


def materialize(x):
    # xarray DataArray -> Dask-backed DataArray -> numpy
    if hasattr(x, "compute"):
        x = x.compute()
    if hasattr(x, "values"):
        x = x.values
    return np.asarray(x)


def get_ds_array(ds, candidates: Sequence[str], label: str):
    for name in candidates:
        try:
            if name in ds:
                return ds[name]
        except Exception:
            pass
        try:
            if hasattr(ds, "coords") and name in ds.coords:
                return ds.coords[name]
        except Exception:
            pass
    available = []
    try:
        available = list(ds.variables)
    except Exception:
        pass
    fail(f"Could not locate {label} in Ag3 dataset. Candidates={candidates}; available={available}")


def string_array(x) -> np.ndarray:
    arr = materialize(x)
    flat = arr.astype(object, copy=False)
    conv = np.vectorize(
        lambda z: z.decode("utf-8") if isinstance(z, (bytes, bytearray, np.bytes_)) else str(z),
        otypes=[object],
    )
    return conv(flat)


def normalize_allele(z) -> str:
    if isinstance(z, (bytes, bytearray, np.bytes_)):
        z = z.decode("utf-8")
    if z is None:
        return ""
    s = str(z)
    if s in {"None", "nan", "NaN", "."}:
        return ""
    return s


def extract_haplotype_arrays(ds):
    pos_da = get_ds_array(ds, ["variant_position", "POS"], "variant positions")
    allele_da = get_ds_array(ds, ["variant_allele", "alleles"], "variant alleles")
    sample_da = get_ds_array(ds, ["sample_id", "samples"], "sample IDs")
    gt_da = get_ds_array(ds, ["call_genotype", "genotype"], "phased genotype calls")

    # Prefer named dimension order where xarray supplies it.
    try:
        dims = list(gt_da.dims)
        target = []
        for canonical, aliases in [
            ("variants", {"variants", "variant"}),
            ("samples", {"samples", "sample"}),
            ("ploidy", {"ploidy", "haplotypes", "haplotype"}),
        ]:
            hit = next((d for d in dims if d.lower() in aliases), None)
            if hit is not None:
                target.append(hit)
        if len(target) == 3 and target != dims:
            gt_da = gt_da.transpose(*target)
    except Exception:
        pass

    pos = materialize(pos_da).astype(np.int64).reshape(-1)
    samples = string_array(sample_da).reshape(-1).astype(str)
    gt = materialize(gt_da)
    alleles = materialize(allele_da)

    if gt.ndim != 3:
        fail(f"Expected call_genotype to be 3-D (variants,samples,ploidy); observed shape {gt.shape}")
    if gt.shape[0] != len(pos):
        fail(f"Genotype variant axis {gt.shape[0]} != position count {len(pos)}")
    if gt.shape[1] != len(samples):
        fail(f"Genotype sample axis {gt.shape[1]} != sample count {len(samples)}")
    if gt.shape[2] < 1 or gt.shape[2] > 2:
        fail(f"Expected phased ploidy of 1 or 2; observed {gt.shape[2]}")
    if alleles.ndim != 2 or alleles.shape[0] != len(pos):
        fail(f"Unexpected variant_allele shape {alleles.shape}; positions={len(pos)}")

    refs: list[str] = []
    alts: list[str] = []
    for row in alleles:
        vals = [normalize_allele(z) for z in row]
        if not vals or not vals[0]:
            fail("A haplotype variant has a missing reference allele.")
        refs.append(vals[0])
        alt = [z for z in vals[1:] if z]
        alts.append(",".join(alt))

    return pos, np.asarray(refs, dtype=object), np.asarray(alts, dtype=object), samples, gt


def reorder_samples(gt: np.ndarray, observed: np.ndarray, desired: np.ndarray) -> np.ndarray:
    if np.array_equal(observed, desired):
        return gt
    if len(observed) != len(desired) or set(observed.tolist()) != set(desired.tolist()):
        missing = sorted(set(desired.tolist()) - set(observed.tolist()))[:20]
        extra = sorted(set(observed.tolist()) - set(desired.tolist()))[:20]
        fail(f"Phased-panel sample membership changed across target batches. missing={missing}; extra={extra}")
    index = {s: i for i, s in enumerate(observed.tolist())}
    idx = [index[s] for s in desired.tolist()]
    return gt[:, idx, :]


def gt_token(pair: np.ndarray) -> str:
    vals = []
    for z in pair.tolist():
        try:
            iz = int(z)
        except Exception:
            iz = -1
        vals.append("." if iz < 0 else str(iz))
    if len(vals) == 1:
        return vals[0]
    return "|".join(vals[:2])


def write_genotype_tsv(
    path: Path,
    contig: str,
    positions: np.ndarray,
    refs: np.ndarray,
    alts: np.ndarray,
    sample_ids: np.ndarray,
    gt: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with tmp.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(FIXED_GT_COLUMNS + sample_ids.tolist())
            for i in range(len(positions)):
                ref = str(refs[i])
                alt = str(alts[i])
                vid = f"{contig}:{int(positions[i])}:{ref}:{alt}"
                row = [vid, contig, int(positions[i]), ref, alt, "PHASED_PANEL"]
                row.extend(gt_token(gt[i, j, :]) for j in range(gt.shape[1]))
                w.writerow(row)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)


def normalized_taxon_for_filter(x: pd.Series) -> pd.Series:
    s = x.astype(str).str.strip().str.lower()
    s = s.str.replace(r"^an[.]?\s*", "", regex=True)
    s = s.str.replace("anopheles ", "", regex=False)
    s = s.str.replace("a. ", "", regex=False)
    s = s.str.replace("an. ", "", regex=False)
    s = s.str.replace("gambiae s.s.", "gambiae", regex=False)
    return s


def resolve_column(df: pd.DataFrame, candidates: Sequence[str], required: bool = True) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        fail(f"None of the required columns {candidates} were found. Available={list(df.columns)}")
    return None


def create_client():
    try:
        return malariagen_data.Ag3()
    except Exception as exc:
        fail(
            "Could not initialize malariagen_data.Ag3().\n"
            f"Original error: {exc}\n\n"
            "If this is a Google Cloud 401/403 error, authenticate once outside R:\n"
            "  gcloud auth application-default login\n"
            "Then rerun this acquisition script. No large VCF download is required."
        )


def acquire_metadata(ag3, release: str, out_dir: Path, exclude_sets: set[str]) -> pd.DataFrame:
    out = out_dir / "metadata" / f"ag3_{release}_sample_metadata.csv"
    sets_out = out_dir / "metadata" / f"ag3_{release}_sample_sets.csv"
    log(f"Fetching Ag3 {release} sample metadata (small table, not sequence data)")
    try:
        md = ag3.sample_metadata(sample_sets=release)
    except Exception as exc:
        fail(
            f"Ag3 sample metadata request failed: {exc}\n"
            "If the message contains 401/403, run 'gcloud auth application-default login' and retry."
        )
    if not isinstance(md, pd.DataFrame) or md.empty:
        fail("Ag3 sample_metadata() returned no rows.")
    sample_col = resolve_column(md, ["sample_id"])
    if md[sample_col].astype(str).duplicated().any():
        fail("Ag3 metadata contains duplicate sample_id values.")
    set_col = resolve_column(md, ["sample_set"], required=False)
    if set_col and exclude_sets:
        md = md[~md[set_col].astype(str).isin(exclude_sets)].copy()
    atomic_write_dataframe(md, out)

    try:
        ss = ag3.sample_sets(release=release)
        if isinstance(ss, pd.DataFrame):
            if "sample_set" in ss.columns and exclude_sets:
                ss = ss[~ss["sample_set"].astype(str).isin(exclude_sets)].copy()
            atomic_write_dataframe(ss, sets_out)
    except Exception as exc:
        log(f"WARNING: sample_sets(release={release}) could not be saved: {exc}")

    log(f"Saved metadata: {out} ({len(md):,} samples)")
    return md



def discover_variant_bearing_windows(
    ag3,
    contig: str,
    windows: list[Interval],
    analysis: str,
    max_regions_per_request: int,
) -> tuple[list[Interval], int]:
    """Use the sites-only hierarchy first, avoiding genotype reads for invariant targets."""
    bearing: list[Interval] = []
    unique_positions: set[int] = set()
    batches = list(batched(windows, max_regions_per_request))
    for bi, batch in enumerate(batches, start=1):
        regions = [x.region for x in batch]
        region_arg = regions[0] if len(regions) == 1 else regions
        log(f"{contig}: sites-only prefilter batch {bi}/{len(batches)}; {len(regions)} target window(s)")
        try:
            arr = ag3.haplotype_sites(
                region=region_arg,
                field="POS",
                analysis=analysis,
                inline_array=True,
                chunks="native",
            )
            pos = materialize(arr).astype(np.int64).reshape(-1)
        except Exception as exc:
            fail(
                f"Sites-only phased-SNP prefilter failed for {contig}, batch {bi}: {exc}\n"
                "The workflow will NOT silently fall back to fetching genotypes for every target window, "
                "because that would increase data transfer. Resolve access/authentication and rerun."
            )
        # Integrity guard: a multi-region API response must contain only sites
        # inside the exact regions requested in this batch. This also protects the
        # transfer budget from an unexpected whole-contig response.
        if len(pos):
            inside = np.zeros(len(pos), dtype=bool)
            for w in batch:
                inside |= (pos >= w.start) & (pos <= w.end)
            if not bool(np.all(inside)):
                bad = pos[~inside][:20].tolist()
                fail(
                    f"Sites-only API returned positions outside the requested target windows for {contig}: {bad}. "
                    "Aborting rather than risk an unexpectedly broad data request."
                )

        unique_positions.update(int(x) for x in pos.tolist())
        if len(pos):
            for w in batch:
                if np.any((pos >= w.start) & (pos <= w.end)):
                    bearing.append(w)
    # Windows are non-overlapping; preserve deterministic order and de-duplicate defensively.
    seen = set()
    out = []
    for w in bearing:
        key = (w.contig, w.start, w.end)
        if key not in seen:
            seen.add(key)
            out.append(w)
    return out, len(unique_positions)


def acquire_haplotypes_for_contig(
    ag3,
    contig: str,
    windows: list[Interval],
    release: str,
    analysis: str,
    sample_query: str,
    max_regions_per_request: int,
    out_dir: Path,
    force: bool,
    signature: str,
    sites_first_prefilter: bool,
) -> dict:
    hap_dir = out_dir / "haplotypes"
    gt_path = hap_dir / f"{contig}.phased_genotypes.tsv"
    done_path = hap_dir / f"{contig}.done.json"
    if not force and gt_path.exists() and done_path.exists():
        try:
            done = json.loads(done_path.read_text(encoding="utf-8"))
            if done.get("signature") == signature and done.get("sha256") == sha256_file(gt_path):
                log(f"Reusing validated cached haplotypes for {contig}: {gt_path.name}")
                return done
        except Exception:
            pass

    if not windows:
        fail(f"No haplotype target windows were planned for contig {contig}.")

    all_pos: list[np.ndarray] = []
    all_ref: list[np.ndarray] = []
    all_alt: list[np.ndarray] = []
    all_gt: list[np.ndarray] = []
    master_samples: np.ndarray | None = None

    if sites_first_prefilter:
        genotype_windows, n_prefilter_positions = discover_variant_bearing_windows(
            ag3, contig, windows, analysis, max_regions_per_request
        )
        log(
            f"{contig}: sites-only prefilter found {n_prefilter_positions:,} phased SNP position(s) "
            f"across {len(genotype_windows):,}/{len(windows):,} target window(s)"
        )
    else:
        genotype_windows = windows
        n_prefilter_positions = -1

    # If every frozen target on the contig is invariant in the phased-sites index,
    # query one 23-bp window only to establish the phased-panel sample axis.
    query_windows = genotype_windows if genotype_windows else [windows[0]]
    batches = list(batched(query_windows, max_regions_per_request))
    for bi, batch in enumerate(batches, start=1):
        regions = [x.region for x in batch]
        region_arg = regions[0] if len(regions) == 1 else regions
        log(f"{contig}: phased haplotypes batch {bi}/{len(batches)}; {len(regions)} target window(s); {sum(x.span for x in batch):,} logical bp")
        try:
            ds = ag3.haplotypes(
                region=region_arg,
                analysis=analysis,
                sample_sets=release,
                sample_query=sample_query,
                inline_array=True,
                chunks="native",
            )
            pos, ref, alt, samples, gt = extract_haplotype_arrays(ds)
        except Exception as exc:
            fail(
                f"Targeted haplotype request failed for {contig}, batch {bi}.\n"
                f"Regions (first 10): {regions[:10]}\nOriginal error: {exc}\n\n"
                "No chromosome-scale VCF was requested. If this is 401/403, run "
                "'gcloud auth application-default login' and retry; completed contigs are cached."
            )

        if master_samples is None:
            master_samples = samples
        else:
            gt = reorder_samples(gt, samples, master_samples)
        all_pos.append(pos)
        all_ref.append(ref)
        all_alt.append(alt)
        all_gt.append(gt)

    assert master_samples is not None
    if all_pos:
        pos = np.concatenate(all_pos) if all_pos else np.asarray([], dtype=np.int64)
        ref = np.concatenate(all_ref) if all_ref else np.asarray([], dtype=object)
        alt = np.concatenate(all_alt) if all_alt else np.asarray([], dtype=object)
        gt = np.concatenate(all_gt, axis=0) if all_gt else np.empty((0, len(master_samples), 2), dtype=np.int8)
    else:
        pos = np.asarray([], dtype=np.int64)
        ref = np.asarray([], dtype=object)
        alt = np.asarray([], dtype=object)
        gt = np.empty((0, len(master_samples), 2), dtype=np.int8)

    # Sort deterministically and remove exact duplicate variants that could only arise
    # from an upstream region-concatenation implementation detail.
    if len(pos):
        order = np.lexsort((alt.astype(str), ref.astype(str), pos))
        pos, ref, alt, gt = pos[order], ref[order], alt[order], gt[order, :, :]
        keep = np.ones(len(pos), dtype=bool)
        for i in range(1, len(pos)):
            same = pos[i] == pos[i - 1] and str(ref[i]) == str(ref[i - 1]) and str(alt[i]) == str(alt[i - 1])
            if same:
                if not np.array_equal(gt[i], gt[i - 1]):
                    fail(f"Duplicate variant {contig}:{pos[i]} has conflicting genotype arrays.")
                keep[i] = False
        pos, ref, alt, gt = pos[keep], ref[keep], alt[keep], gt[keep, :, :]

    # Integrity guard: every returned genotype variant must lie inside one of
    # the target windows actually requested after the sites-only prefilter.
    if len(pos):
        inside = np.zeros(len(pos), dtype=bool)
        for w in query_windows:
            inside |= (pos >= w.start) & (pos <= w.end)
        if not bool(np.all(inside)):
            bad = pos[~inside][:20].tolist()
            fail(
                f"Haplotype API returned positions outside requested target windows for {contig}: {bad}. "
                "Aborting rather than accepting a broader-than-planned response."
            )

    # Guardrail: phased panel for this workflow is SNP-only.
    if any(len(str(x)) != 1 for x in ref.tolist()):
        fail(f"Non-SNP REF allele detected in targeted phased data for {contig}.")
    for x in alt.tolist():
        if x and any(len(a) != 1 for a in str(x).split(",")):
            fail(f"Non-SNP ALT allele detected in targeted phased data for {contig}.")

    write_genotype_tsv(gt_path, contig, pos, ref, alt, master_samples, gt)
    done = {
        "schema_version": ACQUISITION_SCHEMA_VERSION,
        "signature": signature,
        "contig": contig,
        "n_windows": len(windows),
        "logical_bp": int(sum(x.span for x in windows)),
        "sites_first_prefilter": bool(sites_first_prefilter),
        "sample_query": sample_query,
        "n_variant_bearing_windows": int(len(genotype_windows)),
        "n_prefilter_variant_positions": int(n_prefilter_positions),
        "genotype_query_logical_bp": int(sum(x.span for x in query_windows)),
        "n_variants": int(len(pos)),
        "n_samples": int(len(master_samples)),
        "path": str(gt_path.relative_to(Path.cwd())),
        "bytes": int(gt_path.stat().st_size),
        "sha256": sha256_file(gt_path),
    }
    atomic_write_json(done_path, done)
    log(f"{contig}: saved {len(pos):,} targeted phased variant(s) for {len(master_samples):,} panel sample(s); {gt_path.stat().st_size/1024/1024:.2f} MiB")
    return done


def acquire_accessibility_for_contig(
    ag3,
    contig: str,
    windows: list[Interval],
    sites_contig: pd.DataFrame,
    site_mask: str,
    out_dir: Path,
    force: bool,
    signature: str,
) -> dict:
    acc_dir = out_dir / "accessibility"
    out = acc_dir / f"{contig}.target_accessibility.tsv"
    done_path = acc_dir / f"{contig}.done.json"
    if not force and out.exists() and done_path.exists():
        try:
            done = json.loads(done_path.read_text(encoding="utf-8"))
            if done.get("signature") == signature and done.get("sha256") == sha256_file(out):
                log(f"Reusing validated cached accessibility for {contig}: {out.name}")
                return done
        except Exception:
            pass

    rows: list[dict] = []
    seen_sites: set[str] = set()
    for wi, window in enumerate(windows, start=1):
        log(f"{contig}: accessibility window {wi}/{len(windows)} {window.region} ({window.span:,} bp)")
        try:
            arr = materialize(ag3.is_accessible(region=window.region, site_mask=site_mask, inline_array=True, chunks="native")).astype(bool).reshape(-1)
        except Exception as exc:
            fail(
                f"Targeted accessibility request failed for {window.region}: {exc}\n"
                "No chromosome-scale object was requested. If this is 401/403, authenticate with "
                "'gcloud auth application-default login' and rerun."
            )
        if len(arr) != window.span:
            fail(f"Accessibility length mismatch for {window.region}: expected {window.span}, observed {len(arr)}")

        hit = sites_contig[(sites_contig.genomic_start >= window.start) & (sites_contig.genomic_end <= window.end)]
        for _, s in hit.iterrows():
            sid = str(s.population_site_id)
            if sid in seen_sites:
                continue
            start = int(s.genomic_start)
            end = int(s.genomic_end)
            v = arr[(start - window.start) : (end - window.start + 1)]
            if len(v) != 23:
                fail(f"Expected 23 accessibility positions for {sid}; observed {len(v)}")
            strand = str(s.guide_genomic_strand)
            pam_idx = np.zeros(23, dtype=bool)
            if strand == "+":
                pam_idx[-3:] = True
            elif strand == "-":
                pam_idx[:3] = True
            else:
                fail(f"Invalid guide strand for {sid}: {strand}")
            rows.append(
                {
                    "population_site_id": sid,
                    "contig": contig,
                    "accessibility_records_expected": 23,
                    "accessibility_records_found": 23,
                    "accessibility_complete": True,
                    "accessibility_missing_positions": "",
                    "accessibility_fraction_23bp": float(v.mean()),
                    "pam_accessibility_fraction": float(v[pam_idx].mean()),
                    "protospacer_accessibility_fraction": float(v[~pam_idx].mean()),
                }
            )
            seen_sites.add(sid)

    expected = set(sites_contig.population_site_id.astype(str).tolist())
    missing = sorted(expected - seen_sites)
    if missing:
        fail(f"Accessibility plan failed to cover {len(missing)} target(s) on {contig}; examples={missing[:20]}")
    df = pd.DataFrame(rows).sort_values(["population_site_id"]).reset_index(drop=True)
    atomic_write_dataframe(df, out, sep="\t")
    done = {
        "schema_version": ACQUISITION_SCHEMA_VERSION,
        "signature": signature,
        "contig": contig,
        "n_windows": len(windows),
        "logical_bp": int(sum(x.span for x in windows)),
        "n_sites": int(len(df)),
        "path": str(out.relative_to(Path.cwd())),
        "bytes": int(out.stat().st_size),
        "sha256": sha256_file(out),
    }
    atomic_write_json(done_path, done)
    return done


def main() -> int:
    p = argparse.ArgumentParser(description="Fetch only Ag3 data required by the frozen MosqEditR target panel.")
    p.add_argument("--project-root", default=None)
    p.add_argument("--dry-run", action="store_true", help="Build/validate request plan but perform no network access.")
    p.add_argument("--force", action="store_true", help="Ignore validated local cache and reacquire.")
    p.add_argument("--metadata-only", action="store_true", help="Acquire metadata only; do not fetch target haplotypes/accessibility.")
    args = p.parse_args()

    if not ((3, 10) <= sys.version_info[:2] < (3, 13)):
        fail(
            f"Python {sys.version.split()[0]} is not supported by the pinned malariagen_data environment. "
            "Use Python 3.10, 3.11, or 3.12 (3.12 recommended)."
        )

    root = locate_root(args.project_root)
    os.chdir(root)
    cfg = load_config(root)
    local_cfg = cfg_get(cfg, "local_acquisition", default={}) or {}

    frozen = root / "data_processed" / "00_frozen_target_sites.csv"
    if not frozen.exists():
        fail("Missing data_processed/00_frozen_target_sites.csv. Run R Step 00 first: source('R/00_preflight_and_freeze.R')")

    release = str(local_cfg.get("release", cfg_get(cfg, "ag3_release", default="3.0")))
    analysis = str(local_cfg.get("phasing_analysis", cfg_get(cfg, "phasing_panel", default="gamb_colu")))
    site_mask = str(local_cfg.get("site_mask", cfg_get(cfg, "site_filter", default="gamb_colu")))
    canonical = [str(x) for x in cfg_get(cfg, "canonical_contigs", default=["2R", "2L", "3R", "3L", "X"])]
    exclude_sets = set(str(x) for x in cfg_get(cfg, "exclude_sample_set", default=[]))
    primary_taxa = [str(x).strip().lower() for x in cfg_get(cfg, "primary_taxa", default=["gambiae", "coluzzii"]) if str(x).strip()]
    if not primary_taxa:
        fail("analysis_config.yml primary_taxa is empty; refusing an unfiltered genotype request.")

    hap_gap = int(local_cfg.get("haplotype_merge_gap_bp", 0))
    acc_gap = int(local_cfg.get("accessibility_merge_gap_bp", 500))
    max_span = int(local_cfg.get("max_window_span_bp", 25000))
    batch_n = int(local_cfg.get("max_regions_per_haplotype_request", 50))
    max_hap_bp = int(local_cfg.get("max_haplotype_logical_bp", 250000))
    max_acc_bp = int(local_cfg.get("max_accessibility_logical_bp", 2500000))
    max_output_gb = float(local_cfg.get("max_local_output_gb", 2.0))
    sites_first_prefilter = bool(local_cfg.get("sites_first_prefilter", True))

    if hap_gap < 0 or acc_gap < 0 or max_span < 23 or batch_n < 1:
        fail("Invalid local_acquisition interval settings in analysis_config.yml.")

    sites = pd.read_csv(frozen)
    required_cols = {
        "population_site_id", "gene_id", "genomic_seqid", "genomic_start", "genomic_end",
        "guide_genomic_strand", "ag3_phased_panel_assessable",
    }
    missing_cols = sorted(required_cols - set(sites.columns))
    if missing_cols:
        fail(f"Frozen target table is missing columns: {missing_cols}")
    sites["ag3_phased_panel_assessable"] = normalize_bool_series(sites["ag3_phased_panel_assessable"])
    query = sites[(sites.genomic_seqid.astype(str).isin(canonical)) & sites.ag3_phased_panel_assessable].copy()
    query["genomic_seqid"] = query.genomic_seqid.astype(str)
    query["genomic_start"] = query.genomic_start.astype(int)
    query["genomic_end"] = query.genomic_end.astype(int)
    if query.empty:
        fail("No canonical frozen targets are assessable.")
    if ((query.genomic_end - query.genomic_start + 1) != 23).any():
        fail("Target acquisition requires all queried frozen targets to be exactly 23 bp.")

    hap_windows = merge_intervals(query, gap_bp=hap_gap, max_span_bp=max_span)
    acc_windows = merge_intervals(query, gap_bp=acc_gap, max_span_bp=max_span)
    hap_bp = sum(x.span for x in hap_windows)
    acc_bp = sum(x.span for x in acc_windows)
    if hap_bp > max_hap_bp:
        fail(f"Safety budget exceeded: haplotype plan requests {hap_bp:,} logical bp > configured maximum {max_hap_bp:,}.")
    if acc_bp > max_acc_bp:
        fail(f"Safety budget exceeded: accessibility plan requests {acc_bp:,} logical bp > configured maximum {max_acc_bp:,}.")

    out_dir = root / str(local_cfg.get("output_dir", "data_raw/ag3_local"))
    out_dir.mkdir(parents=True, exist_ok=True)
    plan = pd.concat([make_plan_df(hap_windows, "haplotype"), make_plan_df(acc_windows, "accessibility")], ignore_index=True)
    atomic_write_dataframe(plan, out_dir / "target_request_plan.tsv", sep="\t")

    target_hash = sha256_file(frozen)
    signature_payload = {
        "schema_version": ACQUISITION_SCHEMA_VERSION,
        "target_sha256": target_hash,
        "release": release,
        "phasing_analysis": analysis,
        "site_mask": site_mask,
        "canonical_contigs": canonical,
        "exclude_sample_sets": sorted(exclude_sets),
        "primary_taxa": primary_taxa,
        "haplotype_gap_bp": hap_gap,
        "accessibility_gap_bp": acc_gap,
        "max_window_span_bp": max_span,
        "sites_first_prefilter": sites_first_prefilter,
    }
    signature = sha256_text(json.dumps(signature_payload, sort_keys=True))

    log(f"Project: {root}")
    log(f"Frozen targets: {len(sites):,}; canonical targets to query: {len(query):,}")
    log(f"Haplotype request plan: {len(hap_windows):,} small window(s), {hap_bp:,} logical bp total")
    log(f"Accessibility request plan: {len(acc_windows):,} small window(s), {acc_bp:,} logical bp total")
    log("Safety design: no whole chromosome/arm VCF is requested or downloaded")
    if args.dry_run:
        log("DRY RUN complete: no network access performed")
        return 0

    ag3 = create_client()
    try:
        phasing_ids = tuple(str(x) for x in ag3.phasing_analysis_ids)
        if analysis not in phasing_ids:
            fail(f"Configured phasing analysis '{analysis}' is unavailable. Available={phasing_ids}")
    except SystemExit:
        raise
    except Exception as exc:
        fail(f"Could not inspect Ag3 phasing analyses: {exc}")
    try:
        mask_ids = tuple(str(x) for x in ag3.site_mask_ids)
        if site_mask not in mask_ids:
            fail(f"Configured site mask '{site_mask}' is unavailable. Available={mask_ids}")
    except SystemExit:
        raise
    except Exception as exc:
        fail(f"Could not inspect Ag3 site masks: {exc}")

    md = acquire_metadata(ag3, release, out_dir, exclude_sets)

    # Restrict genotype payload to the prespecified primary taxa. The Ag3 haplotype
    # API evaluates this query against sample metadata before materializing genotype
    # arrays, so non-primary samples are not transferred unnecessarily.
    md_taxon_col = resolve_column(md, ["taxon", "aim_species", "species_gambcolu_arabiensis", "species"], required=True)
    if md_taxon_col != "taxon":
        fail(
            f"Targeted genotype transfer requires the standard Ag3 metadata column 'taxon'; "
            f"observed taxon-like column '{md_taxon_col}'. Update the acquisition adapter before continuing rather than fetching a broader sample panel."
        )
    sample_query = "taxon in " + repr(primary_taxa)
    if exclude_sets and "sample_set" in md.columns:
        sample_query += " and sample_set not in " + repr(sorted(exclude_sets))
    log(f"Genotype sample query: {sample_query}")

    if args.metadata_only:
        log("Metadata-only acquisition completed successfully")
        return 0

    hap_done = []
    acc_done = []
    for contig in canonical:
        hws = [x for x in hap_windows if x.contig == contig]
        aws = [x for x in acc_windows if x.contig == contig]
        scontig = query[query.genomic_seqid == contig].copy()
        if scontig.empty:
            continue
        contig_signature = sha256_text(signature + "|" + contig)
        hap_done.append(
            acquire_haplotypes_for_contig(
                ag3, contig, hws, release, analysis, sample_query, batch_n, out_dir, args.force, contig_signature,
                sites_first_prefilter
            )
        )
        acc_done.append(
            acquire_accessibility_for_contig(
                ag3, contig, aws, scontig, site_mask, out_dir, args.force, contig_signature
            )
        )

    # Combine contig accessibility files into the single local handoff table used by R.
    acc_parts = []
    for d in acc_done:
        pth = root / d["path"]
        acc_parts.append(pd.read_csv(pth, sep="\t"))
    acc_all = pd.concat(acc_parts, ignore_index=True)
    if acc_all.population_site_id.astype(str).duplicated().any():
        fail("Combined accessibility output contains duplicate population_site_id values.")
    if set(acc_all.population_site_id.astype(str)) != set(query.population_site_id.astype(str)):
        fail("Combined accessibility output does not exactly match canonical frozen target IDs.")
    acc_all = acc_all.sort_values(["contig", "population_site_id"]).reset_index(drop=True)
    combined_acc = out_dir / "accessibility" / "target_accessibility.tsv"
    atomic_write_dataframe(acc_all, combined_acc, sep="\t")

    # Cross-contig sample-axis integrity from TSV headers.
    sample_axes: dict[str, list[str]] = {}
    for d in hap_done:
        pth = root / d["path"]
        with pth.open("r", encoding="utf-8") as f:
            header = f.readline().rstrip("\n\r").split("\t")
        sample_axes[d["contig"]] = header[len(FIXED_GT_COLUMNS):]
    first_contig = next(iter(sample_axes))
    master = sample_axes[first_contig]
    for contig, axis in sample_axes.items():
        if axis != master:
            if set(axis) != set(master):
                fail(f"Phased-panel sample membership differs between {first_contig} and {contig}.")
            fail(f"Phased-panel sample order differs between {first_contig} and {contig}; acquisition should have been deterministic.")

    total_bytes = sum((root / d["path"]).stat().st_size for d in hap_done + acc_done)
    total_bytes += combined_acc.stat().st_size
    metadata_path = out_dir / "metadata" / f"ag3_{release}_sample_metadata.csv"
    total_bytes += metadata_path.stat().st_size
    if total_bytes > max_output_gb * (1024**3):
        fail(
            f"Local output size {total_bytes/1024**3:.3f} GiB exceeded configured safety limit {max_output_gb:.3f} GiB. "
            "Files were targeted, but inspect before continuing."
        )

    qc_rows = [
        ("acquisition_schema_version", ACQUISITION_SCHEMA_VERSION),
        ("ag3_release", release),
        ("phasing_analysis", analysis),
        ("site_mask", site_mask),
        ("frozen_targets_total", len(sites)),
        ("canonical_targets_queried", len(query)),
        ("haplotype_windows", len(hap_windows)),
        ("haplotype_logical_bp", hap_bp),
        ("accessibility_windows", len(acc_windows)),
        ("accessibility_logical_bp", acc_bp),
        ("phased_panel_samples", len(master)),
        ("sites_first_prefilter", str(sites_first_prefilter).upper()),
        ("variant_bearing_haplotype_windows", sum(int(d.get("n_variant_bearing_windows", 0)) for d in hap_done)),
        ("genotype_query_logical_bp", sum(int(d.get("genotype_query_logical_bp", 0)) for d in hap_done)),
        ("local_output_bytes", total_bytes),
        ("local_output_mib", round(total_bytes / 1024**2, 3)),
        ("whole_chromosome_vcf_downloaded", "NO"),
    ]
    qc = pd.DataFrame(qc_rows, columns=["metric", "value"])
    atomic_write_dataframe(qc, out_dir / "acquisition_qc.tsv", sep="\t")

    manifest = {
        **signature_payload,
        "signature": signature,
        "python": sys.version,
        "platform": platform.platform(),
        "malariagen_data_version": getattr(malariagen_data, "__version__", "unknown"),
        "n_metadata_samples_after_exclusions": int(len(md)),
        "n_phased_panel_samples": int(len(master)),
        "haplotypes": hap_done,
        "accessibility": acc_done,
        "combined_accessibility": {
            "path": str(combined_acc.relative_to(root)),
            "sha256": sha256_file(combined_acc),
            "bytes": combined_acc.stat().st_size,
        },
        "local_output_bytes": total_bytes,
        "whole_chromosome_vcf_downloaded": False,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    atomic_write_json(out_dir / "acquisition_manifest.json", manifest)
    atomic_write_text(out_dir / "ACQUISITION_COMPLETE.ok", signature + "\n")

    log(f"Targeted acquisition complete. Local handoff size: {total_bytes/1024**2:.2f} MiB")
    log("No 50-GB chromosome VCF was downloaded.")
    log("Next: source('R/01_acquire_ag3_metadata.R') in R; Steps 01-02 are now local-only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
