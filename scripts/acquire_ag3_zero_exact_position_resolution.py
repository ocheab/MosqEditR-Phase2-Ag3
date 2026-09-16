#!/usr/bin/env python3
"""Position-resolved follow-up for zero-exact Ag3 CRISPR targets.

This script is deliberately separate from the completed Ag3 pilot acquisition.
It imports and reuses the exact HTTP-range/Tabix/BGZF machinery from
``scripts/acquire_ag3_external_validation.py`` and NEVER overwrites the pilot.

Scientific purpose
------------------
The completed pilot collapses the 23 coordinate-level genotype calls into
sample/site summaries (exact_23bp, pam_intact, called/nonreference allele
counts). This follow-up re-queries ONLY target/taxon combinations whose
``exact_23bp_fraction == 0`` and retains the position-level evidence.

Guide-coordinate convention
---------------------------
The convention is intentionally identical to the completed acquisition:

* guide positions 1..20 = protospacer
* guide positions 21..23 = PAM
* '+' strand: genomic_start -> guide position 1
* '-' strand: genomic_end   -> guide position 1

The script reconstructs the original sample/site states from the new
position-level rows. A completion marker is written only if those reconstructed
states and the taxon summaries exactly reproduce the completed pilot (within a
small floating-point tolerance).

Expected location
-----------------
Place this file beside ``acquire_ag3_external_validation.py`` in the project's
``scripts`` directory and run it from the project root, for example:

    python .\\scripts\\acquire_ag3_zero_exact_position_resolution.py

Outputs are written to:

    data_raw/ag3_external/pilot_position_resolution/

The completed pilot directory is read-only from this script's perspective.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests

try:
    import acquire_ag3_external_validation as base
except ImportError as exc:  # pragma: no cover - runtime placement guard
    raise SystemExit(
        "ERROR: Could not import acquire_ag3_external_validation.py. "
        "Place this script in the same scripts directory as the completed "
        "acquisition engine."
    ) from exc


SCHEMA_VERSION = "ag3-zero-exact-position-resolution-v1"
TOL = 1e-12


def log(msg: str) -> None:
    base.log(msg)


def fail(msg: str, code: int = 2) -> "NoReturn":
    base.fail(msg, code)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[])


def read_tsv_numeric(path: Path) -> pd.DataFrame:
    # For pilot summary/state files, let pandas infer numerics while retaining
    # booleans/strings in their original representation.
    return pd.read_csv(path, sep="\t")


def normalize_bool_value(x: Any) -> Any:
    if pd.isna(x):
        return np.nan
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    s = str(x).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}:
        return True
    if s in {"false", "f", "0", "no", "n"}:
        return False
    if s in {"", "na", "nan", "none", "null"}:
        return np.nan
    raise ValueError(f"cannot normalize boolean value {x!r}")


def bool_series(s: pd.Series) -> pd.Series:
    return s.map(normalize_bool_value)


def bool_equal(a: Any, b: Any) -> bool:
    if pd.isna(a) and pd.isna(b):
        return True
    if pd.isna(a) or pd.isna(b):
        return False
    return bool(a) == bool(b)


def float_equal(a: Any, b: Any, tol: float = TOL) -> bool:
    if pd.isna(a) and pd.isna(b):
        return True
    if pd.isna(a) or pd.isna(b):
        return False
    return abs(float(a) - float(b)) <= tol


def integer_equal(a: Any, b: Any) -> bool:
    if pd.isna(a) and pd.isna(b):
        return True
    if pd.isna(a) or pd.isna(b):
        return False
    return int(a) == int(b)


def reverse_complement(seq: str | None) -> str | None:
    if seq is None or pd.isna(seq):
        return None
    s = str(seq).upper()
    if not s:
        return None
    # Preserve symbolic/non-DNA alleles unchanged. This matters for VCF ALT
    # fields such as <*> or <NON_REF>.
    if not re.fullmatch(r"[ACGTN]+", s):
        return s
    table = str.maketrans("ACGTN", "TGCAN")
    return s.translate(table)[::-1]


def orient_allele(seq: str | None, strand: str) -> str | None:
    if seq is None or pd.isna(seq):
        return None
    s = str(seq)
    if strand == "+":
        return s.upper()
    if strand == "-":
        return reverse_complement(s)
    raise ValueError(f"unexpected guide strand {strand!r}")


def extract_gt_string(format_field: str, sample_field: str) -> str | None:
    fmt = str(format_field).split(":")
    vals = str(sample_field).split(":")
    try:
        idx = fmt.index("GT")
    except ValueError:
        return None
    if idx >= len(vals):
        return None
    gt = vals[idx]
    return gt if gt else None


def alt_list(alt_field: str | None) -> list[str]:
    if alt_field is None or pd.isna(alt_field):
        return []
    s = str(alt_field)
    if not s or s == ".":
        return []
    return s.split(",")


def allele_from_index(idx: int | None, ref: str | None, alt_field: str | None) -> str | None:
    if idx is None:
        return None
    if idx == 0:
        return None if ref is None else str(ref)
    alts = alt_list(alt_field)
    j = int(idx) - 1
    if j < 0 or j >= len(alts):
        return None
    return alts[j]


def genotype_class(indices: list[int | None]) -> str:
    called = [x for x in indices if x is not None]
    if not called:
        return "uncalled"
    nonref = sum(int(x) > 0 for x in called)
    if len(called) == 1:
        return "haploid_reference" if nonref == 0 else "haploid_nonreference"
    if len(called) == 2:
        if nonref == 0:
            return "hom_reference"
        if nonref == 1:
            return "heterozygous"
        return "hom_nonreference"
    if nonref == 0:
        return "polyploid_reference"
    if nonref == len(called):
        return "polyploid_nonreference"
    return "polyploid_mixed"


def expected_guide_position(start: int, end: int, genomic_pos: int, strand: str) -> int:
    if strand == "+":
        return int(genomic_pos - start + 1)
    if strand == "-":
        return int(end - genomic_pos + 1)
    raise ValueError(f"unexpected guide strand {strand!r}")


def region_for_guide_position(guide_position: int) -> tuple[str, int | None, int | None, int | None]:
    gp = int(guide_position)
    if 1 <= gp <= 20:
        # Objective PAM proximity: position 20 is immediately adjacent to PAM.
        return "protospacer", gp, None, 20 - gp
    if 21 <= gp <= 23:
        return "PAM", None, gp - 20, None
    raise ValueError(f"guide position outside 1..23: {gp}")


def load_inputs(root: Path) -> dict[str, Any]:
    cfg = base.load_cfg(root)
    ev = cfg.get("external_validation", {}) or {}
    if not bool(ev.get("enabled", True)):
        fail("external_validation.enabled is false")
    if bool(ev.get("require_no_discovery_retuning", True)) is not True:
        fail("external-validation isolation must remain TRUE")

    out_base = root / str(ev.get("output_dir", "data_raw/ag3_external"))
    pilot_dir = out_base / "pilot"
    if not (pilot_dir / "ACQUISITION_COMPLETE.ok").exists():
        fail(f"Completed pilot marker is missing: {pilot_dir / 'ACQUISITION_COMPLETE.ok'}")

    needed = {
        "pilot_states": pilot_dir / "sample_site_states.tsv",
        "pilot_taxon": pilot_dir / "site_taxon_summary.tsv",
        "pilot_targets": pilot_dir / "locked_targets.tsv",
        "pilot_manifest": pilot_dir / "run_manifest.json",
    }
    for label, path in needed.items():
        if not path.exists():
            fail(f"Missing completed-pilot input ({label}): {path}")

    lock_file = root / "data_processed" / "05b_external_validation_targets.csv"
    if not lock_file.exists():
        fail(f"Missing immutable external-validation lock: {lock_file}")

    manifest = json.loads(needed["pilot_manifest"].read_text(encoding="utf-8"))
    current_lock_hash = base.sha256_file(lock_file)
    manifest_hash = str(manifest.get("target_lock_sha256") or "")
    if manifest_hash and current_lock_hash != manifest_hash:
        fail(
            "Immutable target lock SHA-256 differs from the completed pilot manifest. "
            "Refusing position-resolution re-query."
        )

    targets = pd.read_csv(lock_file)
    pilot_targets = pd.read_csv(needed["pilot_targets"], sep="\t")
    pilot_states = read_tsv_numeric(needed["pilot_states"])
    pilot_taxon = read_tsv_numeric(needed["pilot_taxon"])

    required_target_cols = {
        "population_site_id", "gene_id", "genomic_seqid", "genomic_start",
        "genomic_end", "guide_genomic_strand", "selection_locked_before_ag3",
        "selection_reason", "protospacer_20nt", "pam",
    }
    missing = sorted(required_target_cols - set(targets.columns))
    if missing:
        fail(f"External target lock lacks required position-resolution columns: {missing}")

    if targets.population_site_id.duplicated().any():
        fail("Duplicate population_site_id in immutable target lock")
    if not targets.selection_locked_before_ag3.astype(bool).all():
        fail("External target lock is not immutable for all targets")
    if not ((targets.genomic_end.astype(int) - targets.genomic_start.astype(int) + 1) == 23).all():
        fail("Every locked external target must span exactly 23 bp")

    # Verify that the pilot's frozen copy and current immutable lock describe the
    # same coordinates/strand for every target. Sequence columns are also checked.
    compare_cols = [
        "population_site_id", "genomic_seqid", "genomic_start", "genomic_end",
        "guide_genomic_strand", "protospacer_20nt", "pam",
    ]
    a = targets[compare_cols].copy()
    b = pilot_targets[compare_cols].copy()
    merged = a.merge(b, on="population_site_id", how="outer", suffixes=("_lock", "_pilot"), indicator=True)
    if not (merged["_merge"] == "both").all():
        fail("Pilot locked_targets.tsv and immutable target lock contain different target IDs")
    for col in compare_cols[1:]:
        x = merged[f"{col}_lock"].astype(str)
        y = merged[f"{col}_pilot"].astype(str)
        if not (x == y).all():
            bad = merged.loc[x != y, "population_site_id"].head(10).tolist()
            fail(f"Pilot target lock differs from immutable lock for {col}: {bad}")

    # Derive zero-exact combinations directly from the completed pilot summary;
    # this removes any dependency on a separately exported classification CSV.
    pilot_taxon["exact_23bp_fraction"] = pd.to_numeric(
        pilot_taxon["exact_23bp_fraction"], errors="coerce"
    )
    pilot_taxon["pam_intact_fraction"] = pd.to_numeric(
        pilot_taxon["pam_intact_fraction"], errors="coerce"
    )
    pilot_taxon["callable_sample_fraction"] = pd.to_numeric(
        pilot_taxon["callable_sample_fraction"], errors="coerce"
    )
    zero = pilot_taxon.loc[pilot_taxon.exact_23bp_fraction == 0].copy()
    if zero.empty:
        fail("Completed pilot contains no target/taxon combinations with exact_23bp_fraction == 0")

    def classify_zero(r: pd.Series) -> str:
        call = float(r.callable_sample_fraction) if pd.notna(r.callable_sample_fraction) else np.nan
        pam = float(r.pam_intact_fraction) if pd.notna(r.pam_intact_fraction) else np.nan
        if pd.notna(call) and call < 0.90:
            return "LOW_CALLABILITY"
        if pd.notna(pam) and pam == 0:
            return "PAM_LOST"
        if pd.notna(pam) and pam < 0.95:
            return "PAM_POLYMORPHIC"
        if pd.notna(pam) and pam >= 0.95:
            return "PROTOSPACER_POLYMORPHIC_PAM_CONSERVED"
        return "OTHER"

    zero["conservation_class"] = zero.apply(classify_zero, axis=1)
    zero = zero[[
        "population_site_id", "gene_id", "taxon", "callable_sample_fraction",
        "exact_23bp_fraction", "pam_intact_fraction",
        "nonreference_allele_fraction_23bp", "conservation_class",
    ]].drop_duplicates()

    # Successful samples are taken from the completed state matrix, not selected
    # again from the public manifest, ensuring exact cohort identity.
    successful_samples = pilot_states[["sample_id", "taxon"]].drop_duplicates().copy()
    if successful_samples.sample_id.duplicated().any():
        d = successful_samples.loc[successful_samples.sample_id.duplicated(False), "sample_id"].tolist()
        fail(f"A completed-pilot sample is associated with multiple taxa: {d[:10]}")

    return {
        "cfg": cfg,
        "ev": ev,
        "out_base": out_base,
        "pilot_dir": pilot_dir,
        "lock_file": lock_file,
        "lock_hash": current_lock_hash,
        "manifest": manifest,
        "targets": targets,
        "pilot_states": pilot_states,
        "pilot_taxon": pilot_taxon,
        "zero": zero,
        "samples": successful_samples.sort_values(["taxon", "sample_id"]).reset_index(drop=True),
    }


def query_sample_positions(
    session: requests.Session,
    sample_id: str,
    taxon: str,
    targets: pd.DataFrame,
    lock_hash: str,
    lock_file: Path,
    merge_gap: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Re-query only the supplied targets and return one row per sample/site/base."""

    target_maps = base.build_target_maps(targets)
    vu = base.vcf_url(sample_id)
    iu = base.tbi_url(sample_id)
    counter = base.ByteCounter()
    file_size, _ = base.strict_range_probe(session, vu)

    index_blob = base.download_small_index(session, iu, counter)
    refs = base.parse_tbi(index_blob)
    missing_refs = [c for c in target_maps if c not in refs]
    if missing_refs:
        raise RuntimeError(f"TBI missing required contigs: {missing_refs}")

    chunks_by_contig: dict[str, list[tuple[int, int]]] = {}
    query_chunk_total = 0
    for contig, posmap in target_maps.items():
        all_chunks: list[tuple[int, int]] = []
        positions = sorted(posmap)
        intervals: list[tuple[int, int]] = []
        a = b = positions[0]
        for p in positions[1:]:
            if p <= b + 1:
                b = p
            else:
                intervals.append((a, b))
                a = b = p
        intervals.append((a, b))
        for start1, end1 in intervals:
            all_chunks.extend(base.chunks_for_region(refs[contig], start1 - 1, end1))
        query_chunk_total += len(all_chunks)
        chunks_by_contig[contig] = base.merge_virtual_chunks(
            all_chunks, max_compressed_gap=merge_gap
        )

    # Record-level evidence retained here instead of being collapsed immediately.
    records_by_contig: dict[str, dict[int, dict[str, Any]]] = {
        c: {} for c in target_maps
    }
    fetched_ranges = 0

    for contig, chunks in chunks_by_contig.items():
        needed_positions = target_maps[contig]
        for vbeg, vend in chunks:
            if base.sha256_file(lock_file) != lock_hash:
                raise RuntimeError("immutable external target lock changed during position resolution")
            cstart = int(vbeg >> 16)
            cend = int(vend >> 16)
            byte_end = min(file_size - 1, cend + 65535)
            payload = base.fetch_range(session, vu, cstart, byte_end, file_size, counter)
            fetched_ranges += 1
            text = base.decompress_virtual_chunk(payload, cstart, vbeg, vend)

            for raw_line in text.splitlines():
                if not raw_line or raw_line.startswith(b"#"):
                    continue
                try:
                    line = raw_line.decode("utf-8")
                except UnicodeDecodeError:
                    line = raw_line.decode("latin-1")
                fields = line.split("\t")
                if len(fields) < 10:
                    continue
                chrom = fields[0]
                if chrom != contig:
                    continue
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                if pos not in needed_positions:
                    continue

                rec = {
                    "ref": fields[3],
                    "alt": fields[4],
                    "format": fields[8],
                    "sample_field": fields[9],
                    "gt_indices": base.parse_gt(fields[8], fields[9]),
                    "gt_string": extract_gt_string(fields[8], fields[9]),
                }
                old = records_by_contig[contig].get(pos)
                if old is not None:
                    # The completed engine expects no conflicting duplicate records.
                    keys = ("ref", "alt", "gt_indices", "gt_string")
                    if any(old[k] != rec[k] for k in keys):
                        raise RuntimeError(
                            f"{sample_id} {contig}:{pos} has conflicting duplicate VCF records"
                        )
                else:
                    records_by_contig[contig][pos] = rec

    rows: list[dict[str, Any]] = []
    target_record_fracs: list[float] = []

    for r in targets.itertuples(index=False):
        contig = str(r.genomic_seqid)
        start = int(r.genomic_start)
        end = int(r.genomic_end)
        strand = str(r.guide_genomic_strand)
        site_id = str(r.population_site_id)
        gene_id = str(r.gene_id)
        protospacer = str(r.protospacer_20nt).upper()
        pam = str(r.pam).upper()
        locked_23mer = protospacer + pam
        if len(locked_23mer) != 23:
            raise RuntimeError(f"{site_id} protospacer_20nt + pam is not 23 nt")

        record_map = records_by_contig.get(contig, {})
        n_present = 0

        for genomic_pos in range(start, end + 1):
            rec = record_map.get(genomic_pos)
            record_present = rec is not None
            if record_present:
                n_present += 1
                indices = list(rec["gt_indices"])
                gt_string = rec["gt_string"]
                ref = rec["ref"]
                alt = rec["alt"]
            else:
                indices = []
                gt_string = None
                ref = None
                alt = None

            called = [a for a in indices if a is not None]
            n_called = len(called)
            n_nonref = sum(int(a) > 0 for a in called)

            gp = expected_guide_position(start, end, genomic_pos, strand)
            region, prot_pos, pam_pos, distance_to_pam = region_for_guide_position(gp)
            locked_base = locked_23mer[gp - 1]

            # VCF REF should be one AgamP4 base for these all-sites records. If it
            # is not a single base (e.g. an indel anchor), retain the raw REF but
            # do not force a misleading single-base lock comparison.
            ref_guide = None
            ref_consistent = np.nan
            if ref is not None and len(str(ref)) == 1 and re.fullmatch(r"[ACGTNacgtn]", str(ref)):
                ref_guide = orient_allele(str(ref), strand)
                ref_consistent = bool(str(ref_guide).upper() == locked_base)

            allele_strings_genomic = [allele_from_index(a, ref, alt) for a in indices]
            allele_strings_guide = [
                orient_allele(a, strand) if a is not None else None
                for a in allele_strings_genomic
            ]

            rows.append({
                "sample_id": sample_id,
                "taxon": taxon,
                "population_site_id": site_id,
                "gene_id": gene_id,
                "contig": contig,
                "genomic_position": genomic_pos,
                "guide_genomic_strand": strand,
                "guide_position": gp,
                "target_region": region,
                "protospacer_position": prot_pos,
                "pam_position": pam_pos,
                "distance_to_pam": distance_to_pam,
                "locked_23mer": locked_23mer,
                "locked_guide_reference_base": locked_base,
                "record_present": record_present,
                "vcf_ref": ref,
                "vcf_alt": alt,
                "vcf_ref_guide_oriented": ref_guide,
                "reference_consistent_with_lock": ref_consistent,
                "GT": gt_string,
                "allele_indices": "/".join("." if a is None else str(a) for a in indices) if indices else None,
                "called_alleles_genomic": "|".join("." if a is None else str(a) for a in allele_strings_genomic) if indices else None,
                "called_alleles_guide_oriented": "|".join("." if a is None else str(a) for a in allele_strings_guide) if indices else None,
                "n_called_alleles": int(n_called),
                "n_nonreference_alleles": int(n_nonref),
                "genotype_class": genotype_class(indices),
                "any_nonreference": (bool(n_nonref > 0) if n_called > 0 else np.nan),
            })

        target_record_fracs.append(n_present / 23.0)

    qc = {
        "sample_id": sample_id,
        "taxon": taxon,
        "vcf_url": vu,
        "vcf_remote_size_bytes": int(file_size),
        "tbi_bytes_downloaded": int(counter.index),
        "vcf_bytes_downloaded": int(counter.vcf),
        "http_requests": int(counter.requests),
        "tabix_candidate_chunks": int(query_chunk_total),
        "merged_vcf_ranges": int(fetched_ranges),
        "mean_zero_exact_coordinate_record_fraction": float(np.mean(target_record_fracs)),
        "min_zero_exact_target_record_fraction": float(np.min(target_record_fracs)),
        "targets_with_all_23_records": int(sum(x == 1.0 for x in target_record_fracs)),
        "targets_total": int(len(target_record_fracs)),
    }
    return rows, qc


def reconstruct_sample_site(position_states: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    keys = ["sample_id", "taxon", "population_site_id"]
    for (sample_id, taxon, site_id), d in position_states.groupby(keys, sort=False):
        if d.guide_position.nunique() != 23:
            raise RuntimeError(f"{sample_id} {site_id} does not contain exactly 23 guide positions")
        d = d.sort_values("guide_position")
        pam = d[d.target_region == "PAM"]
        callable23 = bool((d.n_called_alleles > 0).all())
        pam_callable = bool(len(pam) == 3 and (pam.n_called_alleles > 0).all())
        exact23 = bool((d.n_nonreference_alleles == 0).all()) if callable23 else np.nan
        pam_intact = bool((pam.n_nonreference_alleles == 0).all()) if pam_callable else np.nan
        rows.append({
            "sample_id": sample_id,
            "taxon": taxon,
            "population_site_id": site_id,
            "callable_23bp_recalc": callable23,
            "exact_23bp_recalc": exact23,
            "pam_callable_recalc": pam_callable,
            "pam_intact_recalc": pam_intact,
            "called_alleles_23bp_recalc": int(d.n_called_alleles.sum()),
            "nonref_alleles_23bp_recalc": int(d.n_nonreference_alleles.sum()),
            "vcf_records_present_23bp_recalc": int(d.record_present.astype(bool).sum()),
        })
    return pd.DataFrame(rows)


def compare_sample_qc(recalc: pd.DataFrame, pilot_states: pd.DataFrame, zero: pd.DataFrame) -> pd.DataFrame:
    old = pilot_states.merge(
        zero[["population_site_id", "taxon"]].drop_duplicates(),
        on=["population_site_id", "taxon"], how="inner"
    ).copy()

    for c in ["callable_23bp", "exact_23bp", "pam_callable", "pam_intact"]:
        old[c] = bool_series(old[c])
    for c in ["called_alleles_23bp", "nonref_alleles_23bp", "vcf_records_present_23bp"]:
        old[c] = pd.to_numeric(old[c], errors="coerce")

    q = old.merge(recalc, on=["sample_id", "taxon", "population_site_id"], how="outer", indicator=True)

    def row_match(r: pd.Series) -> bool:
        if r["_merge"] != "both":
            return False
        return (
            bool_equal(r.callable_23bp, r.callable_23bp_recalc)
            and bool_equal(r.exact_23bp, r.exact_23bp_recalc)
            and bool_equal(r.pam_callable, r.pam_callable_recalc)
            and bool_equal(r.pam_intact, r.pam_intact_recalc)
            and integer_equal(r.called_alleles_23bp, r.called_alleles_23bp_recalc)
            and integer_equal(r.nonref_alleles_23bp, r.nonref_alleles_23bp_recalc)
            and integer_equal(r.vcf_records_present_23bp, r.vcf_records_present_23bp_recalc)
        )

    q["qc_match"] = q.apply(row_match, axis=1)
    return q


def aggregate_recalc_taxon(recalc: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (site_id, taxon), d in recalc.groupby(["population_site_id", "taxon"], sort=False):
        called = d[d.callable_23bp_recalc == True]
        pamc = d[d.pam_callable_recalc == True]
        ca = int(d.called_alleles_23bp_recalc.sum())
        na = int(d.nonref_alleles_23bp_recalc.sum())
        rows.append({
            "population_site_id": site_id,
            "taxon": taxon,
            "n_samples_total_recalc": int(len(d)),
            "n_samples_called_23bp_recalc": int(len(called)),
            "n_samples_exact_23bp_recalc": int(pd.Series(called.exact_23bp_recalc).fillna(False).astype(bool).sum()),
            "n_samples_pam_called_recalc": int(len(pamc)),
            "n_samples_pam_intact_recalc": int(pd.Series(pamc.pam_intact_recalc).fillna(False).astype(bool).sum()),
            "callable_sample_fraction_recalc": float(len(called) / len(d)) if len(d) else np.nan,
            "exact_23bp_fraction_recalc": float(pd.to_numeric(called.exact_23bp_recalc, errors="coerce").mean()) if len(called) else np.nan,
            "pam_intact_fraction_recalc": float(pd.to_numeric(pamc.pam_intact_recalc, errors="coerce").mean()) if len(pamc) else np.nan,
            "called_alleles_23bp_recalc": ca,
            "nonreference_alleles_23bp_recalc": na,
            "nonreference_allele_fraction_23bp_recalc": float(na / ca) if ca else np.nan,
        })
    return pd.DataFrame(rows)


def compare_taxon_qc(recalc_taxon: pd.DataFrame, pilot_taxon: pd.DataFrame, zero: pd.DataFrame) -> pd.DataFrame:
    old = pilot_taxon.merge(
        zero[["population_site_id", "taxon"]].drop_duplicates(),
        on=["population_site_id", "taxon"], how="inner"
    ).copy()
    q = old.merge(recalc_taxon, on=["population_site_id", "taxon"], how="outer", indicator=True)

    int_pairs = [
        ("n_samples_total", "n_samples_total_recalc"),
        ("n_samples_called_23bp", "n_samples_called_23bp_recalc"),
        ("n_samples_exact_23bp", "n_samples_exact_23bp_recalc"),
        ("n_samples_pam_called", "n_samples_pam_called_recalc"),
        ("n_samples_pam_intact", "n_samples_pam_intact_recalc"),
        ("called_alleles_23bp", "called_alleles_23bp_recalc"),
        ("nonreference_alleles_23bp", "nonreference_alleles_23bp_recalc"),
    ]
    float_pairs = [
        ("callable_sample_fraction", "callable_sample_fraction_recalc"),
        ("exact_23bp_fraction", "exact_23bp_fraction_recalc"),
        ("pam_intact_fraction", "pam_intact_fraction_recalc"),
        ("nonreference_allele_fraction_23bp", "nonreference_allele_fraction_23bp_recalc"),
    ]

    def row_match(r: pd.Series) -> bool:
        if r["_merge"] != "both":
            return False
        return all(integer_equal(r[a], r[b]) for a, b in int_pairs) and all(
            float_equal(r[a], r[b]) for a, b in float_pairs
        )

    q["qc_match"] = q.apply(row_match, axis=1)
    for a, b in float_pairs:
        q[f"delta_{a}"] = pd.to_numeric(q[a], errors="coerce") - pd.to_numeric(q[b], errors="coerce")
    for a, b in int_pairs:
        q[f"delta_{a}"] = pd.to_numeric(q[a], errors="coerce") - pd.to_numeric(q[b], errors="coerce")
    return q


def position_summary(position_states: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "population_site_id", "gene_id", "taxon", "contig", "genomic_position",
        "guide_genomic_strand", "guide_position", "target_region",
        "protospacer_position", "pam_position", "distance_to_pam",
        "locked_guide_reference_base", "vcf_ref", "vcf_ref_guide_oriented", "vcf_alt",
    ]
    rows: list[dict[str, Any]] = []

    for keys, d in position_states.groupby(group_cols, dropna=False, sort=False):
        called = d[d.n_called_alleles > 0]
        called_alleles = int(d.n_called_alleles.sum())
        nonref = int(d.n_nonreference_alleles.sum())
        ref_checks = pd.to_numeric(d.reference_consistent_with_lock, errors="coerce")
        row = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        row.update({
            "n_samples_total": int(d.sample_id.nunique()),
            "n_samples_called": int(called.sample_id.nunique()),
            "called_alleles": called_alleles,
            "nonreference_alleles": nonref,
            "nonreference_allele_frequency": float(nonref / called_alleles) if called_alleles else np.nan,
            "n_samples_any_nonreference": int((pd.to_numeric(d.n_nonreference_alleles, errors="coerce") > 0).sum()),
            "fraction_called_samples_any_nonreference": float((called.n_nonreference_alleles > 0).mean()) if len(called) else np.nan,
            "n_hom_reference": int((d.genotype_class == "hom_reference").sum()),
            "n_heterozygous": int((d.genotype_class == "heterozygous").sum()),
            "n_hom_nonreference": int((d.genotype_class == "hom_nonreference").sum()),
            "n_haploid_reference": int((d.genotype_class == "haploid_reference").sum()),
            "n_haploid_nonreference": int((d.genotype_class == "haploid_nonreference").sum()),
            "reference_consistency_fraction": float(ref_checks.mean()) if ref_checks.notna().any() else np.nan,
        })
        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["taxon", "population_site_id", "guide_position"]
    ).reset_index(drop=True)


def site_position_summary(pos: pd.DataFrame, zero: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (site_id, gene_id, taxon), d in pos.groupby(
        ["population_site_id", "gene_id", "taxon"], sort=False
    ):
        d = d.sort_values("guide_position")
        prot = d[d.target_region == "protospacer"]
        pam = d[d.target_region == "PAM"]
        af = pd.to_numeric(d.nonreference_allele_frequency, errors="coerce")
        paf = pd.to_numeric(prot.nonreference_allele_frequency, errors="coerce")
        pamaf = pd.to_numeric(pam.nonreference_allele_frequency, errors="coerce")

        polymorphic = d.loc[af > 0, "guide_position"].astype(int).tolist()
        prot_poly = prot.loc[paf > 0, "guide_position"].astype(int).tolist()
        pam_poly = pam.loc[pamaf > 0, "guide_position"].astype(int).tolist()
        fixed = d.loc[af >= (1 - TOL), "guide_position"].astype(int).tolist()
        high = d.loc[af >= 0.95, "guide_position"].astype(int).tolist()

        finite = af[np.isfinite(af)]
        if len(finite):
            max_idx = af.idxmax()
            max_row = d.loc[max_idx]
            max_af = float(af.loc[max_idx])
            max_gp = int(max_row.guide_position)
            max_region = str(max_row.target_region)
        else:
            max_af = np.nan
            max_gp = np.nan
            max_region = None

        rows.append({
            "population_site_id": site_id,
            "gene_id": gene_id,
            "taxon": taxon,
            "n_positions_called": int((pd.to_numeric(d.called_alleles, errors="coerce") > 0).sum()),
            "n_polymorphic_positions_23bp": len(polymorphic),
            "n_polymorphic_protospacer_positions": len(prot_poly),
            "n_polymorphic_pam_positions": len(pam_poly),
            "polymorphic_guide_positions": ";".join(map(str, polymorphic)),
            "protospacer_polymorphic_guide_positions": ";".join(map(str, prot_poly)),
            "pam_polymorphic_guide_positions": ";".join(map(str, pam_poly)),
            "fixed_nonreference_guide_positions": ";".join(map(str, fixed)),
            "high_af_ge_0_95_guide_positions": ";".join(map(str, high)),
            "max_nonreference_af_23bp": float(finite.max()) if len(finite) else np.nan,
            "max_nonreference_af_protospacer": float(paf[np.isfinite(paf)].max()) if np.isfinite(paf).any() else np.nan,
            "max_nonreference_af_pam": float(pamaf[np.isfinite(pamaf)].max()) if np.isfinite(pamaf).any() else np.nan,
            "most_divergent_guide_position": max_gp,
            "most_divergent_region": max_region,
            "most_divergent_nonreference_af": max_af,
        })

    out = pd.DataFrame(rows)
    merge_cols = [
        "population_site_id", "taxon", "callable_sample_fraction",
        "exact_23bp_fraction", "pam_intact_fraction",
        "nonreference_allele_fraction_23bp", "conservation_class",
    ]
    out = out.merge(zero[merge_cols], on=["population_site_id", "taxon"], how="left")
    return out.sort_values(["taxon", "population_site_id"]).reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Position-resolved follow-up of Ag3 target/taxon combinations with "
            "exact_23bp_fraction == 0. The completed pilot is never overwritten."
        )
    )
    ap.add_argument("--project-root", default=None)
    ap.add_argument("--dry-run", action="store_true", help="validate inputs only; do not query VCFs")
    ap.add_argument("--force", action="store_true", help="replace an existing position-resolution output directory")
    args = ap.parse_args()

    root = base.locate_root(args.project_root)
    os.chdir(root)
    x = load_inputs(root)

    out_dir = x["out_base"] / "pilot_position_resolution"
    marker = out_dir / "POSITION_RESOLUTION_COMPLETE.ok"
    if marker.exists() and not args.force:
        fail(f"Position resolution is already complete: {marker}. Use --force only for an intentional rerun.")
    if out_dir.exists() and any(out_dir.iterdir()) and not args.force:
        fail(f"Output directory is non-empty: {out_dir}. Use --force to replace it.")
    if args.force and out_dir.exists():
        for p in out_dir.iterdir():
            if p.is_file() or p.is_symlink():
                p.unlink()
            elif p.is_dir():
                import shutil
                shutil.rmtree(p)
    out_dir.mkdir(parents=True, exist_ok=True)

    zero = x["zero"].copy()
    samples = x["samples"].copy()
    targets = x["targets"].copy()

    zero_ids = set(zero.population_site_id.astype(str))
    zero_targets = targets[targets.population_site_id.astype(str).isin(zero_ids)].copy()
    if len(zero_targets) != len(zero_ids):
        missing = sorted(zero_ids - set(zero_targets.population_site_id.astype(str)))
        fail(f"Zero-exact target IDs missing from immutable target lock: {missing}")

    log(f"Project: {root}")
    log(f"Completed pilot: {x['pilot_dir']}")
    log(f"Position-resolution output: {out_dir}")
    log(
        f"Zero-exact combinations: {len(zero):,}; unique targets: {zero.population_site_id.nunique():,}; "
        f"successful pilot samples: {samples.sample_id.nunique():,}"
    )
    for tax, d in zero.groupby("taxon", sort=True):
        log(f"  {tax}: zero-exact combinations={len(d):,}; successful samples={(samples.taxon == tax).sum():,}")
    for cls, d in zero.groupby("conservation_class", sort=True):
        log(f"  class {cls}: {len(d):,}")

    base.atomic_df(zero, out_dir / "zero_exact_combinations.tsv")
    base.atomic_df(zero_targets, out_dir / "zero_exact_locked_targets.tsv")
    base.atomic_df(samples, out_dir / "successful_pilot_samples.tsv")

    if args.dry_run:
        log("DRY RUN complete: inputs validated; no VCF/index genomic data opened")
        return 0

    merge_gap = int(x["ev"].get("tabix_merge_compressed_gap_bytes", 65536) or 65536)
    session = requests.Session()
    session.headers.update({"User-Agent": "MosqEditR-Ag3ZeroExactPositionResolution/1.0"})

    all_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    error_rows: list[dict[str, str]] = []

    n_samples = len(samples)
    for idx, sm in samples.iterrows():
        sample_id = str(sm.sample_id)
        taxon = str(sm.taxon)
        zero_ids_tax = set(zero.loc[zero.taxon.astype(str) == taxon, "population_site_id"].astype(str))
        tt = zero_targets[zero_targets.population_site_id.astype(str).isin(zero_ids_tax)].copy()
        if tt.empty:
            continue

        log(f"[{idx + 1}/{n_samples}] {sample_id} ({taxon}) — {len(tt)} zero-exact target(s)")
        try:
            rows, qc = query_sample_positions(
                session=session,
                sample_id=sample_id,
                taxon=taxon,
                targets=tt,
                lock_hash=x["lock_hash"],
                lock_file=x["lock_file"],
                merge_gap=merge_gap,
            )
            all_rows.extend(rows)
            qc_rows.append(qc)
            log(
                f"  recovered all 23 VCF records for {qc['targets_with_all_23_records']}/"
                f"{qc['targets_total']} targets; VCF={qc['vcf_bytes_downloaded']/1024/1024:.2f} MiB; "
                f"TBI={qc['tbi_bytes_downloaded']/1024:.1f} KiB; ranges={qc['merged_vcf_ranges']}"
            )
        except Exception as exc:
            error_rows.append({
                "sample_id": sample_id,
                "taxon": taxon,
                "vcf_url": base.vcf_url(sample_id),
                "error": repr(exc),
            })
            log(f"WARNING: position resolution failed for {sample_id}: {exc}")

    errors = pd.DataFrame(error_rows, columns=["sample_id", "taxon", "vcf_url", "error"])
    qcdf = pd.DataFrame(qc_rows)
    base.atomic_df(errors, out_dir / "position_resolution_errors.tsv")
    base.atomic_df(qcdf, out_dir / "position_resolution_transfer_qc.tsv")

    if error_rows:
        fail(
            f"{len(error_rows)} sample(s) failed position resolution. "
            f"Inspect {out_dir / 'position_resolution_errors.tsv'}. No completion marker was written."
        )

    position_states = pd.DataFrame(all_rows)
    if position_states.empty:
        fail("No position-level rows were generated")

    expected_position_rows = 0
    for _, sm in samples.iterrows():
        n_targets_tax = int((zero.taxon.astype(str) == str(sm.taxon)).sum())
        expected_position_rows += n_targets_tax * 23
    if len(position_states) != expected_position_rows:
        fail(
            f"Position matrix incomplete: rows={len(position_states):,}; expected={expected_position_rows:,}"
        )

    position_states = position_states.sort_values(
        ["taxon", "population_site_id", "sample_id", "guide_position"]
    ).reset_index(drop=True)
    base.atomic_df(position_states, out_dir / "sample_position_states_zero_exact.tsv")

    # Position-level summaries.
    pos_summary = position_summary(position_states)
    base.atomic_df(pos_summary, out_dir / "site_taxon_position_summary.tsv")

    resolved = site_position_summary(pos_summary, zero)
    base.atomic_df(resolved, out_dir / "zero_exact_position_resolved_summary.tsv")

    # Strict reconstruction against the frozen pilot outputs.
    recalc = reconstruct_sample_site(position_states)
    sample_qc = compare_sample_qc(recalc, x["pilot_states"], zero)
    base.atomic_df(sample_qc, out_dir / "position_resolution_sample_qc.tsv")

    recalc_taxon = aggregate_recalc_taxon(recalc)
    taxon_qc = compare_taxon_qc(recalc_taxon, x["pilot_taxon"], zero)
    base.atomic_df(taxon_qc, out_dir / "position_resolution_taxon_qc.tsv")

    sample_fail = int((sample_qc.qc_match != True).sum())
    taxon_fail = int((taxon_qc.qc_match != True).sum())

    # REF-vs-lock orientation QC. A non-single-base VCF REF is intentionally NA
    # and is not counted as a mismatch.
    ref_qc = position_states.copy()
    ref_qc["reference_consistent_with_lock_numeric"] = pd.to_numeric(
        ref_qc.reference_consistent_with_lock, errors="coerce"
    )
    comparable = ref_qc[ref_qc.reference_consistent_with_lock_numeric.notna()].copy()
    ref_mismatches = comparable[comparable.reference_consistent_with_lock_numeric != 1].copy()
    base.atomic_df(ref_mismatches, out_dir / "reference_lock_mismatches.tsv")

    summary = {
        "schema_version": SCHEMA_VERSION,
        "created": base.now(),
        "source_acquisition_schema": x["manifest"].get("schema_version"),
        "source_target_lock_sha256": x["lock_hash"],
        "source_pilot_dir": str(x["pilot_dir"]),
        "output_dir": str(out_dir),
        "zero_taxon_site_combinations": int(len(zero)),
        "unique_zero_exact_targets": int(zero.population_site_id.nunique()),
        "successful_pilot_samples": int(samples.sample_id.nunique()),
        "position_rows": int(len(position_states)),
        "position_summary_rows": int(len(pos_summary)),
        "sample_qc_rows": int(len(sample_qc)),
        "sample_qc_mismatches": sample_fail,
        "taxon_qc_rows": int(len(taxon_qc)),
        "taxon_qc_mismatches": taxon_fail,
        "reference_lock_comparable_rows": int(len(comparable)),
        "reference_lock_mismatches": int(len(ref_mismatches)),
        "whole_vcf_downloaded": False,
        "validation_isolation": True,
        "discovery_retuning_allowed": False,
        "vcf_bytes_downloaded_total": int(qcdf.vcf_bytes_downloaded.sum()) if len(qcdf) else 0,
        "tbi_bytes_downloaded_total": int(qcdf.tbi_bytes_downloaded.sum()) if len(qcdf) else 0,
        "http_requests_total": int(qcdf.http_requests.sum()) if len(qcdf) else 0,
    }
    base.atomic_json(out_dir / "position_resolution_manifest.json", summary)

    log("===== POSITION-RESOLUTION QC =====")
    log(f"sample/site reconstruction mismatches: {sample_fail}")
    log(f"taxon-summary reconstruction mismatches: {taxon_fail}")
    log(f"reference-lock comparable rows: {len(comparable):,}; mismatches: {len(ref_mismatches):,}")

    if sample_fail or taxon_fail:
        fail(
            "Position-level re-query did not exactly reproduce the completed pilot summaries. "
            "Inspect position_resolution_sample_qc.tsv and position_resolution_taxon_qc.tsv. "
            "No completion marker was written."
        )

    if len(ref_mismatches):
        fail(
            "VCF reference alleles do not fully agree with the locked guide-oriented "
            "reference sequence at comparable single-base records. Inspect "
            "reference_lock_mismatches.tsv before interpreting the position-level results. "
            "No completion marker was written."
        )

    marker_text = (
        "AG3 zero-exact position resolution complete\n"
        f"created={summary['created']}\n"
        f"zero_taxon_site_combinations={summary['zero_taxon_site_combinations']}\n"
        f"unique_targets={summary['unique_zero_exact_targets']}\n"
        f"position_rows={summary['position_rows']}\n"
        f"sample_qc_mismatches={sample_fail}\n"
        f"taxon_qc_mismatches={taxon_fail}\n"
        f"reference_lock_mismatches={len(ref_mismatches)}\n"
    )
    base.atomic_text(marker, marker_text)

    log(f"SUCCESS: position-resolved outputs written to {out_dir}")
    log("Key output: site_taxon_position_summary.tsv")
    log("Key output: zero_exact_position_resolved_summary.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
