#!/usr/bin/env python3
"""
Targeted open-access Ag1000G Phase 2 AR1 acquisition for MosqEditR Manuscript 2.

Design principles
-----------------
* Discovery/development data only: Ag1000G Phase 2 AR1 phased haplotypes.
* No Ag3 credentials required.
* No chromosome-scale VCF download.
* Prefer the public legacy GCS Zarr representation, which permits chunked reads.
* Use Sanger public HTTP only as a fallback/discovery endpoint.
* Download the Phase-2 accessibility HDF5 only if it is below the configured
  hard size ceiling; only the 23-bp frozen target intervals are exported.
* Outputs are local, plain CSV/TSV files consumed by R Steps 01-02.

This script is intentionally defensive. Unexpected schema, coordinates outside
frozen targets, duplicate sample IDs, or an unexpectedly large object abort the
run rather than silently broadening data transfer.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence
from urllib.parse import urljoin

try:
    import numpy as np
    import pandas as pd
    import yaml
    import requests
    import fsspec
    import gcsfs
    import zarr
    import h5py
except ImportError as exc:
    raise SystemExit(
        f"Missing Python dependency: {exc}\nRun scripts/setup_hybrid_python.ps1 (Windows) "
        "or scripts/setup_hybrid_python.sh (Linux/WSL/macOS)."
    ) from exc

SCHEMA_VERSION = "phase2-targeted-v1"
FIXED_GT_COLUMNS = ["variant_id", "contig", "position", "ref", "alt", "filter"]
CANONICAL = ("2R", "2L", "3R", "3L", "X")


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
    print(f"{time.strftime('%Y-%m-%dT%H:%M:%S%z')}\t{msg}", flush=True)


def fail(msg: str) -> "None":
    raise SystemExit(f"ERROR: {msg}")


def locate_root(explicit: str | None) -> Path:
    root = Path(explicit).expanduser().resolve() if explicit else Path(__file__).resolve().parents[1]
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


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    os.close(fd)
    p = Path(tmp)
    try:
        p.write_text(text, encoding="utf-8", newline="")
        os.replace(p, path)
    finally:
        p.unlink(missing_ok=True)


def atomic_df(df: pd.DataFrame, path: Path, sep: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    os.close(fd)
    p = Path(tmp)
    try:
        df.to_csv(p, index=False, sep=sep)
        os.replace(p, path)
    finally:
        p.unlink(missing_ok=True)


def atomic_json(path: Path, obj: dict) -> None:
    atomic_text(path, json.dumps(obj, indent=2, sort_keys=True, default=lambda x: int(x) if isinstance(x, np.integer) else x) + "\n")


def merge_intervals(rows: pd.DataFrame, gap_bp: int, max_span_bp: int) -> list[Interval]:
    out: list[Interval] = []
    for contig, d in rows.groupby("genomic_seqid", sort=False):
        iv = sorted(zip(d.genomic_start.astype(int), d.genomic_end.astype(int)))
        if not iv:
            continue
        s0, e0 = iv[0]
        for s, e in iv[1:]:
            pe = max(e0, e)
            if s <= e0 + gap_bp + 1 and pe - s0 + 1 <= max_span_bp:
                e0 = pe
            else:
                out.append(Interval(str(contig), int(s0), int(e0)))
                s0, e0 = s, e
        out.append(Interval(str(contig), int(s0), int(e0)))
    return out


def make_plan(intervals: Sequence[Interval], kind: str) -> pd.DataFrame:
    return pd.DataFrame({
        "kind": kind,
        "window_id": [f"{kind}_{i:04d}" for i in range(1, len(intervals)+1)],
        "contig": [x.contig for x in intervals],
        "start": [x.start for x in intervals],
        "end": [x.end for x in intervals],
        "span_bp": [x.span for x in intervals],
        "region": [x.region for x in intervals],
    })


def get_text(url: str, timeout: int = 60) -> str:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text


def read_phase2_metadata(url: str) -> pd.DataFrame:
    log(f"Downloading small Phase 2 sample metadata table: {url}")
    r = requests.get(url, timeout=90)
    r.raise_for_status()
    from io import StringIO
    md = pd.read_csv(StringIO(r.text), sep="\t")
    if md.empty:
        fail("Phase 2 metadata table is empty.")
    if "ox_code" not in md.columns:
        fail(f"Phase 2 metadata lacks ox_code. Columns={list(md.columns)}")
    if md["ox_code"].astype(str).duplicated().any():
        fail("Phase 2 metadata contains duplicate ox_code sample IDs.")
    return md


def map_phase2_taxon(md: pd.DataFrame) -> pd.Series:
    # Historical Phase 2 M/S nomenclature: M form -> A. coluzzii; S form -> A. gambiae.
    m = md.get("m_s", pd.Series([None] * len(md), index=md.index)).astype(str).str.strip().str.upper()
    out = pd.Series(np.where(m.eq("M"), "coluzzii", np.where(m.eq("S"), "gambiae", np.nan)), index=md.index, dtype="object")
    # Population codes commonly end in M/S and provide a conservative fallback.
    if "population" in md.columns:
        pop = md["population"].astype(str).str.strip().str.upper()
        out = out.where(out.notna(), np.where(pop.str.endswith("M"), "coluzzii", np.where(pop.str.endswith("S"), "gambiae", np.nan)))
    return out


def prepare_metadata(md: pd.DataFrame, primary_taxa: set[str]) -> pd.DataFrame:
    out = md.copy()
    out["sample_id"] = out["ox_code"].astype(str).str.strip()
    out["sample_set"] = "AG1000G-PHASE2-AR1"
    out["taxon"] = map_phase2_taxon(out)
    out["country"] = out["country"].astype(str).str.strip() if "country" in out else np.nan
    out["location"] = out["region"].astype(str).str.strip() if "region" in out else np.nan
    out["year"] = out["year"] if "year" in out else np.nan
    out["phase2_population"] = out["population"].astype(str).str.strip() if "population" in out else np.nan
    out["primary_taxon"] = out["taxon"].astype(str).str.lower().isin(primary_taxa)
    return out


def _candidate_zarr_roots(bucket: str, prefix: str) -> list[str]:
    base = f"{bucket}/{prefix}/haplotypes/main"
    stem = "ag1000g.phase2.ar1.haplotypes"
    return [
        f"{base}/zarr2/zstd/{stem}",
        f"{base}/zarr2/blosc_zstd/{stem}",
        f"{base}/zarr/{stem}",
        f"{base}/zarr/{stem}.zarr",
        f"{base}/zarr2/{stem}",
    ]


def discover_public_gcs_haplotype_root(bucket: str, prefix: str) -> tuple[gcsfs.GCSFileSystem, str]:
    fs = gcsfs.GCSFileSystem(token="anon", cache_timeout=300)
    for root in _candidate_zarr_roots(bucket, prefix):
        try:
            if fs.exists(root + "/.zgroup") or fs.exists(root + "/.zmetadata"):
                log(f"Found public Phase 2 phased Zarr: gs://{root}")
                return fs, root
        except Exception:
            continue
    # Metadata-only listing, bounded depth; do not recursively list chunk objects.
    probe = f"{bucket}/{prefix}/haplotypes/main"
    try:
        found = fs.find(probe, maxdepth=6)
    except Exception as exc:
        fail(f"Could not list the open Phase-2 GCS haplotype prefix: {exc}")
    roots = []
    for p in found:
        if p.endswith("/.zgroup") or p.endswith("/.zmetadata"):
            r = p.rsplit("/", 1)[0]
            if "haplot" in r.lower():
                roots.append(r)
    roots = sorted(set(roots), key=lambda x: ("zstd" not in x.lower(), len(x)))
    if not roots:
        fail(
            "Could not discover a public Phase-2 phased Zarr hierarchy. The script did not fall back to large VCF downloads. "
            "Check public Sanger/GCS availability or use --manual-haplotype-zarr with a local/downloaded Zarr root."
        )
    log(f"Discovered public Phase 2 phased Zarr: gs://{roots[0]}")
    return fs, roots[0]


def open_zarr_group(fs, root: str):
    store = fs.get_mapper(root)
    try:
        # Consolidated metadata reduces network round trips when present.
        if fs.exists(root + "/.zmetadata"):
            return zarr.open_consolidated(store=store, mode="r")
    except Exception:
        pass
    return zarr.open_group(store=store, mode="r")


def find_array(group, paths: Sequence[str], label: str):
    for p in paths:
        try:
            return group[p]
        except Exception:
            pass
    available = []
    try:
        def visitor(name, obj):
            if isinstance(obj, zarr.Array):
                available.append(name)
        group.visititems(visitor)
    except Exception:
        pass
    fail(f"Could not locate {label}. Tried={list(paths)}. Example arrays={available[:50]}")


def decode_strings(a) -> np.ndarray:
    x = np.asarray(a)
    out = []
    for z in x.reshape(-1):
        if isinstance(z, (bytes, bytearray, np.bytes_)):
            out.append(z.decode("utf-8"))
        else:
            out.append(str(z))
    return np.asarray(out, dtype=object).reshape(x.shape)


def allele_vector(group, contig: str, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return REF/ALT strings without assuming a single historical Zarr schema.

    IMPORTANT: do not call ``find_array()`` while probing separate REF/ALT arrays.
    ``find_array()`` intentionally terminates on failure, which would prevent the
    combined-allele fallback used by several Ag1000G Phase-2 Zarr layouts.
    """
    refa = None
    alta = None
    for p in (f"{contig}/variants/REF", f"{contig}/variants/ref"):
        try:
            refa = group[p]
            break
        except Exception:
            pass
    for p in (f"{contig}/variants/ALT", f"{contig}/variants/alt"):
        try:
            alta = group[p]
            break
        except Exception:
            pass

    if refa is not None and alta is not None:
        refs = decode_strings(refa.get_orthogonal_selection((rows,))).reshape(-1)
        aval = decode_strings(alta.get_orthogonal_selection((rows,)))
        if aval.ndim == 1:
            alts = aval.astype(object)
        else:
            alts = np.asarray([
                ",".join([str(v) for v in row if str(v) not in {"", ".", "None", "nan"}])
                for row in aval
            ], dtype=object)
        return refs, alts

    aa = find_array(
        group,
        [f"{contig}/variants/allele", f"{contig}/variants/ALLELE", f"{contig}/variants/alleles"],
        "variant alleles",
    )
    vals = decode_strings(aa.get_orthogonal_selection((rows,)))
    if vals.ndim != 2 or vals.shape[1] < 2:
        fail(f"Unexpected combined allele array shape for {contig}: {vals.shape}")
    refs = vals[:, 0]
    alts = np.asarray([
        ",".join([str(v) for v in row[1:] if str(v) not in {"", ".", "None", "nan"}])
        for row in vals
    ], dtype=object)
    return refs, alts


def extract_rows_for_targets(pos: np.ndarray, target_rows: pd.DataFrame) -> np.ndarray:
    keep: list[np.ndarray] = []
    for s, e in zip(target_rows.genomic_start.astype(int), target_rows.genomic_end.astype(int)):
        lo = int(np.searchsorted(pos, s, side="left"))
        hi = int(np.searchsorted(pos, e, side="right"))
        if hi > lo:
            keep.append(np.arange(lo, hi, dtype=np.int64))
    return np.unique(np.concatenate(keep)) if keep else np.asarray([], dtype=np.int64)


def slice_gt(gt_arr, rows: np.ndarray, sample_idx: np.ndarray) -> np.ndarray:
    if len(rows) == 0:
        return np.empty((0, len(sample_idx), 2), dtype=np.int8)
    try:
        x = gt_arr.get_orthogonal_selection((rows, sample_idx, slice(None)))
    except Exception as exc:
        log(f"WARNING: orthogonal Zarr selection failed ({exc}); using grouped row slices.")
        pieces = []
        # Consecutive target variants are grouped to avoid one request per variant.
        starts = [0]
        for i in range(1, len(rows)):
            if rows[i] != rows[i-1] + 1:
                starts.append(i)
        starts.append(len(rows))
        for a, b in zip(starts[:-1], starts[1:]):
            lo, hi = int(rows[a]), int(rows[b-1]) + 1
            block = np.asarray(gt_arr[lo:hi, :, :])[:, sample_idx, :]
            wanted = rows[a:b] - lo
            pieces.append(block[wanted, :, :])
        x = np.concatenate(pieces, axis=0)
    x = np.asarray(x)
    if x.ndim != 3 or x.shape[0] != len(rows) or x.shape[1] != len(sample_idx):
        fail(f"Unexpected phased genotype slice shape: {x.shape}")
    if x.shape[2] != 2:
        fail(f"Expected diploid phased genotype array with ploidy 2; got {x.shape}")
    return x.astype(np.int16, copy=False)


def gt_token(pair: np.ndarray) -> str:
    vals = []
    for v in pair.tolist():
        try:
            iv = int(v)
        except Exception:
            iv = -1
        vals.append("." if iv < 0 else str(iv))
    return "|".join(vals[:2])


def write_gt(path: Path, contig: str, pos: np.ndarray, refs: np.ndarray, alts: np.ndarray, samples: np.ndarray, gt: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent)); os.close(fd)
    p = Path(tmp)
    try:
        with p.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(FIXED_GT_COLUMNS + samples.astype(str).tolist())
            for i in range(len(pos)):
                ref, alt = str(refs[i]), str(alts[i])
                if len(ref) != 1 or not alt or any(len(a) != 1 for a in alt.split(",")):
                    # Phase 2 phased scaffold is expected to be biallelic SNPs only.
                    continue
                w.writerow([f"{contig}:{int(pos[i])}:{ref}:{alt}", contig, int(pos[i]), ref, alt, "PHASE2_AR1_PHASED"] + [gt_token(gt[i,j,:]) for j in range(gt.shape[1])])
        os.replace(p, path)
    finally:
        p.unlink(missing_ok=True)


def find_accessibility_object(fs: gcsfs.GCSFileSystem, bucket: str, prefix: str) -> str:
    base = f"{bucket}/{prefix}/accessibility"
    preferred = [
        f"{base}/accessibility.h5",
        f"{base}/accessibility.hdf5",
        f"{base}/accessibility/accessibility.h5",
    ]
    for p in preferred:
        try:
            if fs.exists(p):
                return p
        except Exception:
            pass
    try:
        found = fs.find(base, maxdepth=3)
    except Exception as exc:
        fail(f"Could not list Phase-2 accessibility directory: {exc}")
    h5 = [p for p in found if p.lower().endswith((".h5", ".hdf5"))]
    h5.sort(key=lambda p: ("access" not in p.lower(), len(p)))
    if not h5:
        fail("No Phase-2 accessibility HDF5 object could be discovered in the public release.")
    return h5[0]


def download_gcs_guarded(fs, remote: str, local: Path, max_mb: float) -> None:
    info = fs.info(remote)
    size = int(info.get("size", -1))
    if size < 0:
        fail(f"Could not determine remote object size for gs://{remote}; refusing unbounded download.")
    mb = size / (1024**2)
    if mb > max_mb:
        fail(f"Accessibility fallback object is {mb:.1f} MiB, above configured ceiling {max_mb:.1f} MiB. No download performed.")
    if local.exists() and local.stat().st_size == size:
        log(f"Reusing cached accessibility object ({mb:.1f} MiB): {local}")
        return
    local.parent.mkdir(parents=True, exist_ok=True)
    log(f"Downloading Phase-2 accessibility object only ({mb:.1f} MiB; hard ceiling {max_mb:.1f} MiB)")
    fd, tmp = tempfile.mkstemp(prefix=local.name + ".", dir=str(local.parent)); os.close(fd)
    p = Path(tmp)
    try:
        with fs.open(remote, "rb") as src, p.open("wb") as dst:
            while True:
                b = src.read(4 * 1024 * 1024)
                if not b: break
                dst.write(b)
        if p.stat().st_size != size:
            fail(f"Accessibility download size mismatch: expected {size}, received {p.stat().st_size}")
        os.replace(p, local)
    finally:
        p.unlink(missing_ok=True)


def accessibility_dataset(h5, contig: str):
    candidates = [f"{contig}/is_accessible", f"{contig}/accessibility", contig]
    for p in candidates:
        try:
            obj = h5[p]
            if isinstance(obj, h5py.Dataset) and obj.ndim == 1:
                return obj
        except Exception:
            pass
    # Search only within contig group.
    try:
        grp = h5[contig]
        for name, obj in grp.items():
            if isinstance(obj, h5py.Dataset) and obj.ndim == 1 and "access" in name.lower():
                return obj
    except Exception:
        pass
    fail(f"Could not locate 1-D accessibility mask for contig {contig}. HDF5 top-level keys={list(h5.keys())}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Acquire only open Phase-2 phased data needed by the frozen CRISPR target panel.")
    ap.add_argument("--project-root", default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--metadata-only", action="store_true")
    ap.add_argument("--manual-haplotype-zarr", default=None, help="Optional local Zarr root, bypassing network source discovery.")
    args = ap.parse_args()

    root = locate_root(args.project_root)
    os.chdir(root)
    cfg = load_config(root)
    acq = cfg_get(cfg, "phase2_acquisition", default={}) or {}
    out_dir = root / str(acq.get("output_dir", "data_raw/phase2_local"))
    out_dir.mkdir(parents=True, exist_ok=True)
    frozen = root / "data_processed/00_frozen_target_sites.csv"
    if not frozen.exists():
        fail("Missing data_processed/00_frozen_target_sites.csv. Run R Step 00 first: source('R/00_preflight_and_freeze.R')")
    sites = pd.read_csv(frozen)
    required = {"population_site_id","genomic_seqid","genomic_start","genomic_end","guide_genomic_strand"}
    miss = sorted(required - set(sites.columns))
    if miss: fail(f"Frozen target table missing columns: {miss}")
    canonical = set(cfg.get("canonical_contigs", CANONICAL))
    qsites = sites[sites.genomic_seqid.astype(str).isin(canonical)].copy()
    if len(sites) != 530 or len(qsites) != 515:
        fail(f"Frozen-count guard failed: observed total={len(sites)}, canonical={len(qsites)}; expected 530/515.")

    gap_h = int(acq.get("haplotype_merge_gap_bp", 0)); gap_a = int(acq.get("accessibility_merge_gap_bp", 500)); span = int(acq.get("max_window_span_bp", 25000))
    hap_win = merge_intervals(qsites, gap_h, span); acc_win = merge_intervals(qsites, gap_a, span)
    hap_bp = sum(x.span for x in hap_win); acc_bp = sum(x.span for x in acc_win)
    if hap_bp > int(acq.get("max_haplotype_logical_bp", 250000)): fail("Haplotype logical-base safety budget exceeded.")
    if acc_bp > int(acq.get("max_accessibility_logical_bp", 2500000)): fail("Accessibility logical-base safety budget exceeded.")
    log(f"Project: {root}")
    log(f"Frozen targets: {len(sites)}; canonical targets to query: {len(qsites)}")
    log(f"Phase-2 haplotype plan: {len(hap_win)} small window(s), {hap_bp:,} logical bp")
    log(f"Phase-2 accessibility plan: {len(acc_win)} small window(s), {acc_bp:,} logical bp")
    log("Safety design: phased discovery uses public chunked Zarr; no chromosome VCF is downloaded")
    plan = pd.concat([make_plan(hap_win,"haplotype"), make_plan(acc_win,"accessibility")], ignore_index=True)
    atomic_df(plan, out_dir / "request_plan.csv")
    if args.dry_run:
        log("DRY RUN complete: no network access performed")
        return 0

    md_url = str(acq.get("metadata_url"))
    md0 = read_phase2_metadata(md_url)
    primary_taxa = {str(x).lower() for x in cfg.get("primary_taxa", ["gambiae","coluzzii"])}
    md = prepare_metadata(md0, primary_taxa)
    atomic_df(md, out_dir / "metadata" / "sample_metadata.csv")
    primary_md = md[md.primary_taxon == True].copy()
    if primary_md.empty: fail("No gambiae/coluzzii samples were recognized in Phase-2 metadata.")
    log(f"Phase-2 metadata: {len(md):,} rows; primary gambiae/coluzzii={len(primary_md):,}")
    if args.metadata_only:
        atomic_text(out_dir / "METADATA_COMPLETE.ok", SCHEMA_VERSION + "\n")
        return 0

    # Open phased haplotype hierarchy.
    if args.manual_haplotype_zarr:
        zr = Path(args.manual_haplotype_zarr).expanduser().resolve()
        if not zr.exists(): fail(f"Manual Zarr path not found: {zr}")
        group = zarr.open_group(str(zr), mode="r")
        zarr_source = str(zr)
        gfs = None
    else:
        bucket = str(acq.get("public_gcs_bucket", "ag1000g-release")); prefix = str(acq.get("public_gcs_prefix", "phase2.AR1"))
        gfs, zroot = discover_public_gcs_haplotype_root(bucket, prefix)
        group = open_zarr_group(gfs, zroot)
        zarr_source = "gs://" + zroot

    contig_done = []
    selected_axis: np.ndarray | None = None
    for contig in CANONICAL:
        ss = qsites[qsites.genomic_seqid.astype(str) == contig].copy()
        if ss.empty: continue
        pos_arr = find_array(group, [f"{contig}/variants/POS", f"{contig}/variants/position"], f"{contig} positions")
        mb = pos_arr.nbytes / 1024**2
        if mb > float(acq.get("max_position_index_mb_per_contig",150)):
            fail(f"{contig} position index would materialize {mb:.1f} MiB > configured ceiling; no read performed.")
        log(f"{contig}: loading phased-sites position index ({mb:.1f} MiB uncompressed)")
        pos_all = np.asarray(pos_arr[:], dtype=np.int64).reshape(-1)
        if len(pos_all) and np.any(np.diff(pos_all) < 0): fail(f"{contig} phased POS array is not sorted.")
        rows = extract_rows_for_targets(pos_all, ss)
        pos = pos_all[rows]
        if len(pos):
            inside = np.zeros(len(pos), dtype=bool)
            for s,e in zip(ss.genomic_start.astype(int), ss.genomic_end.astype(int)): inside |= (pos>=s)&(pos<=e)
            if not np.all(inside): fail(f"{contig}: target row selector returned positions outside frozen intervals.")

        sample_arr = find_array(group, [f"{contig}/samples", "samples"], f"{contig} sample IDs")
        zsamples = decode_strings(sample_arr[:]).reshape(-1).astype(str)
        if len(set(zsamples.tolist())) != len(zsamples): fail(f"{contig}: duplicate phased sample IDs in Zarr.")
        want = primary_md.sample_id.astype(str).tolist()
        index = {s:i for i,s in enumerate(zsamples.tolist())}
        selected = np.asarray([index[s] for s in want if s in index], dtype=np.int64)
        selected_names = np.asarray([zsamples[i] for i in selected], dtype=object)
        if len(selected_names) < 100: fail(f"{contig}: only {len(selected_names)} primary samples intersect phased panel; unexpected.")
        if selected_axis is None:
            selected_axis = selected_names
        elif not np.array_equal(selected_axis, selected_names):
            fail(f"Phased sample membership/order changed across contigs at {contig}.")

        gt_arr = find_array(group, [f"{contig}/calldata/genotype", f"{contig}/calldata/GT"], f"{contig} phased genotypes")
        if len(gt_arr.shape) != 3: fail(f"{contig}: genotype array shape is {gt_arr.shape}, expected 3-D.")
        gt = slice_gt(gt_arr, rows, selected)
        refs, alts = allele_vector(group, contig, rows) if len(rows) else (np.asarray([],dtype=object),np.asarray([],dtype=object))
        out = out_dir / "haplotypes" / f"{contig}.phased_targets.tsv"
        write_gt(out, contig, pos, refs, alts, selected_names, gt)
        done = {"contig":contig,"n_target_sites":int(len(ss)),"n_variants":int(len(pos)),"n_samples":int(len(selected_names)),"bytes":int(out.stat().st_size),"sha256":sha256_file(out)}
        atomic_json(out_dir / "haplotypes" / f"{contig}.done.json", done)
        contig_done.append(done)
        log(f"{contig}: saved {len(pos):,} phased target-overlap SNPs for {len(selected_names):,} samples")

    if selected_axis is None: fail("No phased samples acquired.")
    # Save exact sample-axis metadata in the order used by genotype files.
    phase_md = md.set_index("sample_id").loc[selected_axis.tolist()].reset_index()
    atomic_df(phase_md, out_dir / "metadata" / "phased_sample_metadata.csv")

    # Accessibility - small, guarded whole-file fallback because HDF5 is not chunk-addressable over HTTP safely.
    if args.manual_haplotype_zarr:
        log("Manual Zarr mode still uses the public Phase-2 accessibility object.")
        bucket = str(acq.get("public_gcs_bucket", "ag1000g-release")); prefix = str(acq.get("public_gcs_prefix", "phase2.AR1")); gfs = gcsfs.GCSFileSystem(token="anon")
    remote_acc = find_accessibility_object(gfs, bucket, prefix)
    local_acc = out_dir / "cache" / Path(remote_acc).name
    download_gcs_guarded(gfs, remote_acc, local_acc, float(acq.get("max_accessibility_fallback_download_mb",750)))
    acc_rows = []
    with h5py.File(local_acc, "r") as h5:
        for contig in CANONICAL:
            ss = qsites[qsites.genomic_seqid.astype(str) == contig]
            if ss.empty: continue
            arr = accessibility_dataset(h5, contig)
            for _, s in ss.iterrows():
                st, en = int(s.genomic_start), int(s.genomic_end)
                if st < 1 or en > arr.shape[0]: fail(f"Accessibility coordinate out of bounds: {s.population_site_id} {contig}:{st}-{en}")
                v = np.asarray(arr[st-1:en], dtype=bool)
                if len(v) != 23: fail(f"Expected 23 accessibility bases for {s.population_site_id}; got {len(v)}")
                pam = np.zeros(23,dtype=bool)
                if str(s.guide_genomic_strand)=="+": pam[-3:]=True
                elif str(s.guide_genomic_strand)=="-": pam[:3]=True
                else: fail(f"Invalid strand for {s.population_site_id}")
                acc_rows.append({"population_site_id":str(s.population_site_id),"contig":contig,"accessibility_records_expected":23,"accessibility_records_found":23,"accessibility_complete":True,"accessibility_missing_positions":"","accessibility_fraction_23bp":float(v.mean()),"pam_accessibility_fraction":float(v[pam].mean()),"protospacer_accessibility_fraction":float(v[~pam].mean())})
    acc = pd.DataFrame(acc_rows)
    if len(acc) != 515 or acc.population_site_id.duplicated().any(): fail(f"Accessibility handoff guard failed; rows={len(acc)}, duplicates={acc.population_site_id.duplicated().sum()}")
    atomic_df(acc, out_dir / "accessibility" / "target_accessibility.tsv", sep="\t")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "discovery_dataset": "Ag1000G Phase 2 AR1",
        "haplotype_source": zarr_source,
        "accessibility_source": "gs://" + remote_acc,
        "n_frozen_targets": int(len(sites)),
        "n_canonical_targets": int(len(qsites)),
        "n_phased_samples": int(len(selected_axis)),
        "contigs": contig_done,
        "no_chromosome_vcf_download": True,
        "validation_isolation": "Ag3 is not used in this acquisition or discovery branch",
    }
    atomic_json(out_dir / "acquisition_manifest.json", manifest)
    atomic_text(out_dir / "ACQUISITION_COMPLETE.ok", SCHEMA_VERSION + "\n")
    log("Phase 2 targeted discovery acquisition completed successfully")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
