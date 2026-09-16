#!/usr/bin/env python3
"""Ag3 external validation via pure-Python Tabix/BGZF HTTP range reads.

This is an EXTERNAL-VALIDATION-ONLY branch. It consumes the immutable Phase-2
target lock and never changes discovery targets, scores, classes, weights or ranks.

Why this engine exists
----------------------
The public Sanger all-sites VCFs and .tbi indexes support HTTP byte ranges, but
remote seeking through Windows Rsamtools/HTSlib produced BGZF "Invalid seek"
errors. This reader avoids HTSlib remote I/O completely:

  small .tbi index download -> Tabix chunk lookup -> exact HTTP byte ranges
  -> local BGZF decompression -> locked 23-bp genotype states.

The full ~3-GB per-sample VCF is never intentionally downloaded.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import struct
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
import requests
import yaml

SCHEMA_VERSION = "ag3-external-public-vcf-purepython-tabix-v2-sequence-aware"
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
        return yaml.safe_load(f) or {}


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


def atomic_df(df: pd.DataFrame, path: Path, sep: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    os.close(fd)
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
    b = str(url).rstrip("/").split("/")[-1]
    for suf in (".gatk.zarr.zip", ".zarr.zip", ".vcf.gz"):
        if b.endswith(suf):
            return b[: -len(suf)]
    return b.split(".")[0]


DNA_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return str(seq).upper().translate(DNA_COMP)[::-1]


def expected_genomic_23(guide_sequence_key: str, strand: str) -> str:
    g = str(guide_sequence_key).upper()
    if len(g) != 23 or re.fullmatch(r"[ACGT]{23}", g) is None:
        raise ValueError(f"invalid frozen guide_sequence_key {guide_sequence_key!r}")
    if strand == "+":
        return g
    if strand == "-":
        return revcomp(g)
    raise ValueError(f"unexpected guide strand {strand!r}")


def _get_csv(session: requests.Session, url: str, sep: str = ",") -> pd.DataFrame:
    last = None
    for attempt in range(1, 5):
        try:
            r = session.get(url, timeout=(15, 90))
            r.raise_for_status()
            return pd.read_csv(io.StringIO(r.text), sep=sep)
        except Exception as exc:
            last = exc
            if attempt < 4:
                time.sleep(min(8, 2 ** (attempt - 1)))
    raise RuntimeError(f"could not read official Ag3 metadata: {url}: {last}")


def _get_text(session: requests.Session, url: str) -> str:
    """Fetch UTF-8-ish publisher text without assuming valid CSV syntax."""
    last = None
    for attempt in range(1, 5):
        try:
            r = session.get(url, timeout=(15, 90))
            r.raise_for_status()
            # requests' apparent_encoding is useful for publisher-generated CSVs
            # containing non-ASCII country/site names.
            if not r.encoding:
                r.encoding = r.apparent_encoding or "utf-8"
            return r.text
        except Exception as exc:
            last = exc
            if attempt < 4:
                time.sleep(min(8, 2 ** (attempt - 1)))
    raise RuntimeError(f"could not fetch Ag3 metadata text: {url}: {last}")



CAPUTO2024_SUPP_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs42003-024-06809-y/MediaObjects/"
    "42003_2024_6809_MOESM4_ESM.csv"
)

# Published Ag3.0 West-African sample composition:
# An. coluzzii: 374; An. gambiae: 449; far-west intermediate taxa are excluded.
CAPUTO2024_EXPECTED_PRIMARY = {"gambiae": 449, "coluzzii": 374}


def _norm_colname(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())


def _find_sample_id_col(d: pd.DataFrame) -> str:
    aliases = {
        "sampleid", "sample", "specimenid", "specimen", "specimen_id",
        "id", "samplecode", "sample_code",
    }
    by_norm = {_norm_colname(c): c for c in d.columns}
    for a in aliases:
        if _norm_colname(a) in by_norm:
            c = by_norm[_norm_colname(a)]
            vals = d[c].astype(str).str.strip()
            hit = vals.str.fullmatch(r"[A-Za-z]{1,6}[0-9]{3,7}-C[WwXx]?").mean()
            if hit >= 0.5:
                return c

    best = None
    best_hit = -1.0
    for c in d.columns:
        vals = d[c].astype(str).str.strip()
        hit = vals.str.fullmatch(r"[A-Za-z]{1,6}[0-9]{3,7}-C[WwXx]?").mean()
        if hit > best_hit:
            best, best_hit = c, float(hit)
    if best is None or best_hit < 0.5:
        raise RuntimeError(
            "Published Ag3.0 supplementary metadata: could not identify the "
            "MalariaGEN sample-ID column."
        )
    return best


def _find_country_col(d: pd.DataFrame) -> str | None:
    by_norm = {_norm_colname(c): c for c in d.columns}
    for k in ("country", "samplingcountry", "collectioncountry", "nation"):
        if k in by_norm:
            return by_norm[k]
    for c in d.columns:
        if "country" in _norm_colname(c):
            return c
    return None


def _infer_population_code_col(d: pd.DataFrame) -> str | None:
    best = None
    best_hit = 0.0
    for c in d.columns:
        s = d[c].astype(str).str.strip()
        # Published paper population codes include BFgam/BFcol etc.
        hit = s.str.fullmatch(
            r"(?i)(?:BF|CI|GH|GN|GW|MA)(?:gam|col)"
        ).mean()
        if float(hit) > best_hit:
            best, best_hit = c, float(hit)
    return best if best_hit >= 0.05 else None


def _infer_primary_taxon_row(row: pd.Series, pop_col: str | None) -> str | None:
    if pop_col is not None:
        x = str(row.get(pop_col, "")).strip().lower()
        if re.fullmatch(r"(?:bf|ci|gh|gn|gw|ma)gam", x):
            return "gambiae"
        if re.fullmatch(r"(?:bf|ci|gh|gn|gw|ma)col", x):
            return "coluzzii"

    # Fallback to explicit species/taxon wording in any metadata field.
    vals = [str(v).strip().lower() for v in row.tolist() if pd.notna(v)]
    joined = " | ".join(vals)

    # Never coerce intermediate/cryptic/far-west categories into the primary taxa.
    if any(k in joined for k in ("intermediate", "gcx1", "gcx2", "gcx3", "bissau")):
        return None

    has_col = bool(re.search(r"\b(?:an\.?\s*)?coluzzii\b", joined))
    has_gam = bool(re.search(r"\b(?:an\.?\s*)?gambiae\b", joined))
    if has_col and not has_gam:
        return "coluzzii"
    if has_gam and not has_col:
        return "gambiae"
    return None


def _country_from_population_code(x: str) -> str | None:
    pfx = str(x).strip().upper()[:2]
    return {
        "BF": "Burkina Faso",
        "CI": "Côte d'Ivoire",
        "GH": "Ghana",
        "GN": "Guinea",
        "GW": "Guinea-Bissau",
        "MA": "Mali",
    }.get(pfx)


def _optional_col(d: pd.DataFrame, keywords: tuple[str, ...]) -> str | None:
    for c in d.columns:
        n = _norm_colname(c)
        if any(k in n for k in keywords):
            return c
    return None


def load_caputo2024_ag3_metadata(
    root: Path,
    session: requests.Session,
) -> pd.DataFrame:
    """Fallback Ag3.0 metadata from Caputo et al. Supplementary Data 1.

    Springer serves the supplementary file with .csv extension, but the file
    contains at least one malformed delimiter row for pandas' strict CSV parser.
    We therefore DO NOT infer biological metadata by skipping malformed rows.

    Instead, we extract only two fields that are explicit and independently
    auditable on every sample row:
      1. Published specimen ID (first semicolon-delimited field; e.g. AA0045-C, AV0209-CW, or VBS02049-4431STDY6672917)
      2. published population code (e.g. BFgam, BFcol)

    The population code deterministically defines both taxon and country.
    Recovery MUST match the published composition exactly:
      - 449 An. gambiae
      - 374 An. coluzzii
    Otherwise the fallback aborts.
    """
    cache = root / "metadata" / "ag3_caputo2024_westafrica_primary_taxa_metadata.csv"
    if cache.exists() and cache.stat().st_size > 1000:
        d = pd.read_csv(cache)
        req = {
            "sample_id", "country", "taxon", "population_code",
            "snp_genotypes_vcf", "metadata_source",
        }
        if req.issubset(d.columns):
            counts = d.groupby("taxon").size().to_dict()
            if (
                int(counts.get("gambiae", 0)) == 449
                and int(counts.get("coluzzii", 0)) == 374
            ):
                return d

    text = _get_text(session, CAPUTO2024_SUPP_URL)

    # Primary West-African populations reported in Table 1 of Caputo et al.
    # Far-west intermediate groups (gcx1-GM, gcx1-GW, gcx2) are deliberately
    # absent from this allowed mapping and therefore cannot enter the cohort.
    pop_info = {
        "BFgam": ("gambiae", "Burkina Faso"),
        "GHgam": ("gambiae", "Ghana"),
        "GNgam": ("gambiae", "Guinea"),
        "GWgam": ("gambiae", "Guinea-Bissau"),
        "MAgam": ("gambiae", "Mali"),
        "BFcol": ("coluzzii", "Burkina Faso"),
        "CIcol": ("coluzzii", "Côte d'Ivoire"),
        "GHcol": ("coluzzii", "Ghana"),
        "GNcol": ("coluzzii", "Guinea"),
        "MAcol": ("coluzzii", "Mali"),
    }
    pop_lookup = {k.lower(): (k, v[0], v[1]) for k, v in pop_info.items()}
    pop_pat = re.compile(
        r"(?i)\b(?:BFgam|GHgam|GNgam|GWgam|MAgam|"
        r"BFcol|CIcol|GHcol|GNcol|MAcol)\b"
    )

    records = []
    ambiguous_lines = []

    for line_no, line in enumerate(text.splitlines(), start=1):
        # The publisher file is semicolon-delimited.  We deliberately avoid
        # pandas' strict CSV parser because embedded comma-separated run lists
        # elsewhere in the file caused tokenisation failure.  The first field
        # is the explicit specimen ID in Supplementary Data 1.
        pops = sorted({m.group(0).lower() for m in pop_pat.finditer(line)})
        if not pops:
            continue

        fields = [x.strip().strip('"').strip("'") for x in line.split(";")]
        sid = fields[0].lstrip("\ufeff").strip() if fields else ""

        # Accept the published specimen ID verbatim rather than inferring one
        # particular MalariaGEN naming convention.  Examples in this table
        # include AA0045-C, AV0209-CW, and VBS02049-4431STDY6672917.
        sid_ok = bool(
            sid
            and 3 <= len(sid) <= 100
            and re.fullmatch(r"[A-Za-z0-9._-]+", sid) is not None
        )

        if len(pops) != 1 or not sid_ok:
            ambiguous_lines.append(
                (
                    line_no,
                    sid,
                    len(pops),
                    line[:180].replace("\t", " ")
                )
            )
            continue

        canonical_pop, taxon, country = pop_lookup[pops[0]]
        records.append(
            {
                "sample_id": sid,
                "sample_set": "Ag3.0-WestAfrica-Caputo2024",
                "country": country,
                "location": pd.NA,
                "year": pd.NA,
                "month": pd.NA,
                "latitude": pd.NA,
                "longitude": pd.NA,
                "sex_call": pd.NA,
                "taxon": taxon,
                "population_code": canonical_pop,
                "snp_genotypes_vcf": (
                    "https://cog.sanger.ac.uk/vo_agam_output/"
                    + sid
                    + ".vcf.gz"
                ),
                "metadata_source": (
                    "Caputo2024 Communications Biology Supplementary Data 1 "
                    "(Ag1000G Ag3.0 West Africa; semicolon-row first-field "
                    "specimen ID + explicit population-code extraction)"
                ),
            }
        )

    out = pd.DataFrame.from_records(records)
    if out.empty:
        raise RuntimeError(
            "Published Ag3.0 fallback produced zero explicit primary-taxon "
            "sample rows."
        )

    if out["sample_id"].duplicated().any():
        dup = (
            out.loc[out["sample_id"].duplicated(keep=False), "sample_id"]
            .drop_duplicates()
            .head(20)
            .tolist()
        )
        raise RuntimeError(
            "Duplicate sample IDs recovered from published Ag3.0 "
            f"supplementary metadata: {dup}"
        )

    counts = out.groupby("taxon").size().to_dict()
    got_gam = int(counts.get("gambiae", 0))
    got_col = int(counts.get("coluzzii", 0))

    # Hard validation against the article's published sample composition.
    if got_gam != CAPUTO2024_EXPECTED_PRIMARY["gambiae"] or got_col != CAPUTO2024_EXPECTED_PRIMARY["coluzzii"]:
        detail = (
            f"recovered gambiae={got_gam}, coluzzii={got_col}; "
            "expected gambiae=449, coluzzii=374"
        )
        if ambiguous_lines:
            detail += (
                f"; ambiguous population-bearing lines={len(ambiguous_lines)}; "
                f"first={ambiguous_lines[:3]}"
            )
        raise RuntimeError(
            "Published Ag3.0 fallback failed exact composition audit: " + detail
        )

    # Population-level count audit against Table 1.
    expected_pop_counts = {
        "BFcol": 135,
        "CIcol": 80,
        "GHcol": 63,
        "GNcol": 11,
        "MAcol": 85,
        "BFgam": 157,
        "GHgam": 36,
        "GNgam": 123,
        "GWgam": 8,
        "MAgam": 125,
    }
    observed_pop_counts = out["population_code"].value_counts().to_dict()
    mismatches = {
        p: (int(observed_pop_counts.get(p, 0)), int(n))
        for p, n in expected_pop_counts.items()
        if int(observed_pop_counts.get(p, 0)) != int(n)
    }
    if mismatches:
        raise RuntimeError(
            "Published Ag3.0 fallback failed population-count audit "
            f"(observed, expected): {mismatches}"
        )

    out = out.sort_values(
        ["taxon", "country", "population_code", "sample_id"],
        kind="mergesort",
    ).reset_index(drop=True)

    cache.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(cache, index=False)

    log(
        "Published Ag3.0 fallback composition audit passed: "
        f"gambiae={got_gam}, coluzzii={got_col}; "
        "all 10 population counts match Caputo et al. Table 1"
    )
    return out

def load_official_ag3_metadata(root: Path) -> pd.DataFrame:
    """Load/caches Ag3.0 primary-taxon metadata.

    Preferred source is MalariaGEN's documented Ag3.0 GCS metadata. If that
    endpoint returns an access error (e.g. HTTP 403), fall back to the published
    Ag3.0 West-African supplementary metadata of Caputo et al. (2024), which
    contains enough explicit An. gambiae and An. coluzzii samples for the locked
    250+250 external-validation design.

    The fallback changes only the external-validation sampling frame; it never
    changes Phase-2 discovery targets, scores, classes, thresholds or ranks.
    """
    cache = root / "metadata" / "ag3_official_v3_primary_taxa_metadata.csv"
    if cache.exists() and cache.stat().st_size > 1000:
        d = pd.read_csv(cache)
        req = {
            "sample_id", "sample_set", "country", "taxon",
            "snp_genotypes_vcf", "metadata_source",
        }
        if req.issubset(d.columns):
            return d

    base = "https://storage.googleapis.com/vo_agam_release_master_us_central1/v3"
    session = requests.Session()
    session.headers.update({"User-Agent": "MosqEditR-Ag3ExternalValidation/2.1"})

    try:
        manifest = _get_csv(session, f"{base}/manifest.tsv", sep="\t")
        if "sample_set" not in manifest.columns:
            raise RuntimeError("official Ag3 manifest lacks sample_set")
        sample_sets = [
            str(x) for x in manifest["sample_set"]
            if str(x) != "AG1000G-X"
        ]

        chunks = []
        for ss in sample_sets:
            general = _get_csv(
                session, f"{base}/metadata/general/{ss}/samples.meta.csv"
            )
            species = _get_csv(
                session,
                f"{base}/metadata/species_calls_aim_20220528/{ss}/"
                "samples.species_aim.csv",
            )
            catalog = _get_csv(
                session, f"{base}/metadata/general/{ss}/wgs_snp_data.csv"
            )
            if "aim_species" not in species.columns:
                raise RuntimeError(
                    f"{ss} official species metadata lacks aim_species"
                )
            if "snp_genotypes_vcf" not in catalog.columns:
                raise RuntimeError(
                    f"{ss} WGS catalog lacks snp_genotypes_vcf"
                )
            d = general.merge(
                species[["sample_id", "aim_species"]],
                on="sample_id", how="inner", validate="one_to_one",
            ).merge(
                catalog[["sample_id", "snp_genotypes_vcf"]],
                on="sample_id", how="inner", validate="one_to_one",
            )
            d["sample_set"] = ss
            d["taxon"] = d["aim_species"].astype(str).str.lower()
            chunks.append(d)

        out = pd.concat(chunks, ignore_index=True)
        out = out[out["taxon"].isin(FIXED_TAXA)].copy()
        out = out[
            out["snp_genotypes_vcf"]
            .astype(str)
            .str.contains(r"\.vcf\.gz$", regex=True)
        ].copy()
        if out["sample_id"].duplicated().any():
            dup = out.loc[
                out["sample_id"].duplicated(), "sample_id"
            ].head().tolist()
            raise RuntimeError(
                f"duplicate sample IDs in official Ag3 metadata: {dup}"
            )

        keep = [
            "sample_id", "sample_set", "country", "location", "year", "month",
            "latitude", "longitude", "sex_call", "taxon", "snp_genotypes_vcf",
        ]
        out = out[[c for c in keep if c in out.columns]].copy()
        out["metadata_source"] = (
            "MalariaGEN Ag3.0 v3 GCS metadata + AIM species calls"
        )
        cache.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(cache, index=False)
        return out

    except Exception as exc:
        log(
            "WARNING: direct MalariaGEN Ag3.0 GCS metadata access failed "
            f"({exc}). Falling back to published Ag3.0 West-African "
            "supplementary metadata."
        )
        out = load_caputo2024_ag3_metadata(root, session)

        # Also cache under the canonical path so repeated runs are reproducible
        # and do not depend on the remote supplementary file remaining reachable.
        out.to_csv(cache, index=False)
        return out

def stable_hash_key(x: str) -> str:
    return hashlib.sha256(str(x).encode("utf-8")).hexdigest()


def _proportional_alloc(counts: pd.Series, n: int) -> dict[str, int]:
    counts = counts.astype(int)
    if n >= int(counts.sum()):
        return {str(k): int(v) for k, v in counts.items()}
    ideal = counts / counts.sum() * n
    alloc = np.floor(ideal).astype(int)
    # When feasible, represent every country at least once.
    if n >= len(counts):
        for k in counts.index:
            if alloc.loc[k] == 0 and counts.loc[k] > 0:
                donors = [j for j in counts.index if alloc.loc[j] > 1]
                if donors:
                    donor = max(donors, key=lambda j: (alloc.loc[j] - ideal.loc[j], alloc.loc[j]))
                    alloc.loc[donor] -= 1
                    alloc.loc[k] += 1
    while int(alloc.sum()) < n:
        candidates = [k for k in counts.index if alloc.loc[k] < counts.loc[k]]
        if not candidates:
            break
        k = max(candidates, key=lambda j: (ideal.loc[j] - alloc.loc[j], counts.loc[j]))
        alloc.loc[k] += 1
    while int(alloc.sum()) > n:
        candidates = [k for k in counts.index if alloc.loc[k] > 0]
        k = max(candidates, key=lambda j: (alloc.loc[j] - ideal.loc[j], alloc.loc[j]))
        alloc.loc[k] -= 1
    return {str(k): int(v) for k, v in alloc.items()}


def stratified_select(md: pd.DataFrame, taxa: Sequence[str], per_taxon: int) -> pd.DataFrame:
    """Deterministic country-stratified selection using SHA-256 within strata."""
    chunks = []
    for taxon in taxa:
        d = md[md.taxon == taxon].copy()
        if d.empty:
            fail(f"No official Ag3.0 samples labelled {taxon!r}.")
        d["country_stratum"] = d["country"].fillna("UNKNOWN").astype(str)
        if per_taxon <= 0 or per_taxon >= len(d):
            take = d
        else:
            alloc = _proportional_alloc(d.groupby("country_stratum").size(), per_taxon)
            picked = []
            for country, n_take in alloc.items():
                if n_take <= 0:
                    continue
                z = d[d.country_stratum == country].copy()
                z["_hash"] = z.sample_id.map(stable_hash_key)
                z = z.sort_values(["_hash", "sample_id"]).head(n_take)
                picked.append(z)
            take = pd.concat(picked, ignore_index=True) if picked else d.iloc[0:0].copy()
        take["selection_strategy"] = "country_stratified_sha256"
        chunks.append(take)
    out = pd.concat(chunks, ignore_index=True)
    return out.drop(columns=["_hash"], errors="ignore")


def load_public_manifest(features_url: str, labels_url: str) -> pd.DataFrame:
    urls = get_lines(features_url)
    labels = get_lines(labels_url)
    if len(urls) != len(labels):
        fail(f"Ag3 public feature/label manifests have different lengths: {len(urls)} vs {len(labels)}")
    if len(urls) < 1000:
        fail(f"Ag3 public feature manifest unexpectedly contains only {len(urls)} rows.")
    out = pd.DataFrame({"source_url": urls, "taxon_raw": labels})
    out["sample_id"] = out.source_url.map(sample_id_from_url)
    out["taxon"] = out.taxon_raw.map(normalize_taxon)
    if out.sample_id.duplicated().any():
        d = out.loc[out.sample_id.duplicated(), "sample_id"].head().tolist()
        fail(f"Duplicate sample IDs in Official Ag3.0 metadata: {d}")
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


def vcf_url(sample_id: str) -> str:
    return f"https://cog.sanger.ac.uk/vo_agam_output/{sample_id}.vcf.gz"


def tbi_url(sample_id: str) -> str:
    return vcf_url(sample_id) + ".tbi"


class ByteCounter:
    def __init__(self):
        self.vcf = 0
        self.index = 0
        self.requests = 0


def request_with_retries(
    session: requests.Session,
    url: str,
    *,
    headers: dict[str, str] | None = None,
    timeout: int = 45,
    stream: bool = False,
    attempts: int = 5,
) -> requests.Response:
    last = None
    for attempt in range(1, attempts + 1):
        try:
            r = session.get(
                url,
                headers=headers,
                timeout=(10, timeout),
                stream=stream,
                allow_redirects=True,
            )
            # Retry transient HTTP failures rather than returning them to the
            # caller immediately.  Do not retry permanent 4xx responses.
            if int(r.status_code) in {408, 425, 429, 500, 502, 503, 504}:
                status = int(r.status_code)
                r.close()
                raise requests.RequestException(
                    f"transient HTTP status {status}"
                )
            return r
        except requests.RequestException as exc:
            last = exc
            if attempt < attempts:
                time.sleep(min(8, 2 ** (attempt - 1)))
    raise RuntimeError(f"HTTP request failed after {attempts} attempts: {url}: {last}")


def strict_range_probe(session: requests.Session, url: str) -> tuple[int, str]:
    r = request_with_retries(
        session, url, headers={"Range": "bytes=0-0"}, timeout=40, stream=True
    )
    try:
        status = int(r.status_code)
        cr = str(r.headers.get("Content-Range") or "")
        first = next(r.iter_content(chunk_size=2), b"")
    finally:
        r.close()
    m = re.match(r"^bytes\s+0-0/(\d+)$", cr, flags=re.I)
    if status != 206 or m is None or len(first) != 1:
        raise RuntimeError(
            f"server did not honor one-byte Range request: HTTP {status}; "
            f"Content-Range={cr!r}; bytes_read={len(first)}"
        )
    return int(m.group(1)), cr


def download_small_index(
    session: requests.Session, url: str, counter: ByteCounter, max_bytes: int = 5 * 1024 * 1024
) -> bytes:
    r = request_with_retries(session, url, timeout=60, stream=True)
    try:
        if r.status_code != 200:
            raise RuntimeError(f"index HTTP status {r.status_code}")
        cl = r.headers.get("Content-Length")
        if cl and int(cl) > max_bytes:
            raise RuntimeError(f"index Content-Length {cl} exceeds {max_bytes} bytes")
        parts = []
        total = 0
        for chunk in r.iter_content(chunk_size=64 * 1024):
            if not chunk:
                continue
            total += len(chunk)
            if total > max_bytes:
                raise RuntimeError(f"index exceeded {max_bytes} bytes while downloading")
            parts.append(chunk)
        counter.index += total
        counter.requests += 1
        return b"".join(parts)
    finally:
        r.close()


def fetch_range(
    session: requests.Session,
    url: str,
    start: int,
    end: int,
    file_size: int,
    counter: ByteCounter,
) -> bytes:
    if start < 0:
        raise ValueError("negative byte range start")
    end = min(int(end), int(file_size) - 1)
    if end < start:
        return b""
    r = request_with_retries(
        session,
        url,
        headers={"Range": f"bytes={int(start)}-{int(end)}"},
        timeout=30,
        stream=False,
    )
    counter.requests += 1
    if r.status_code != 206:
        raise RuntimeError(
            f"VCF range request was not partial: HTTP {r.status_code} for {start}-{end}"
        )
    cr = str(r.headers.get("Content-Range") or "")
    m = re.match(r"^bytes\s+(\d+)-(\d+)/(\d+)$", cr, flags=re.I)
    if m is None:
        raise RuntimeError(f"VCF range response lacks valid Content-Range: {cr!r}")
    rs, re_, total = map(int, m.groups())
    if rs != start or total != file_size:
        raise RuntimeError(
            f"VCF range response mismatch: requested {start}-{end}, got {cr}"
        )
    body = r.content
    expected = re_ - rs + 1
    if len(body) != expected:
        raise RuntimeError(
            f"VCF range body length mismatch: got {len(body)}, expected {expected}"
        )
    counter.vcf += len(body)
    return body


# ------------------------- Tabix index parsing -------------------------

def _read_i32(b: io.BytesIO) -> int:
    x = b.read(4)
    if len(x) != 4:
        raise ValueError("truncated TBI int32")
    return struct.unpack("<i", x)[0]


def _read_u32(b: io.BytesIO) -> int:
    x = b.read(4)
    if len(x) != 4:
        raise ValueError("truncated TBI uint32")
    return struct.unpack("<I", x)[0]


def _read_u64(b: io.BytesIO) -> int:
    x = b.read(8)
    if len(x) != 8:
        raise ValueError("truncated TBI uint64")
    return struct.unpack("<Q", x)[0]


def parse_tbi(blob: bytes) -> dict[str, dict[str, Any]]:
    # .tbi files are normally BGZF/gzip-compressed.
    raw = gzip.decompress(blob) if blob[:2] == b"\x1f\x8b" else blob
    b = io.BytesIO(raw)
    if b.read(4) != b"TBI\x01":
        raise ValueError("Tabix index does not begin with TBI\\1 magic")
    n_ref = _read_i32(b)
    _format = _read_i32(b)
    _col_seq = _read_i32(b)
    _col_beg = _read_i32(b)
    _col_end = _read_i32(b)
    _meta = _read_i32(b)
    _skip = _read_i32(b)
    l_nm = _read_i32(b)
    name_blob = b.read(l_nm)
    names = [x.decode("utf-8") for x in name_blob.split(b"\x00") if x]
    if len(names) != n_ref:
        raise ValueError(f"TBI reference-name count mismatch: header={n_ref}, names={len(names)}")
    refs: dict[str, dict[str, Any]] = {}
    for name in names:
        n_bin = _read_i32(b)
        bins: dict[int, list[tuple[int, int]]] = {}
        for _ in range(n_bin):
            bin_id = _read_u32(b)
            n_chunk = _read_i32(b)
            chunks = [(_read_u64(b), _read_u64(b)) for __ in range(n_chunk)]
            bins[int(bin_id)] = chunks
        n_intv = _read_i32(b)
        linear = [_read_u64(b) for _ in range(n_intv)]
        refs[name] = {"bins": bins, "linear": linear}
    return refs


def reg2bins(beg: int, end: int) -> list[int]:
    """Tabix/BAI bins overlapping 0-based half-open [beg, end)."""
    beg = max(0, int(beg))
    end = max(beg + 1, int(end))
    if end > (1 << 29):
        end = 1 << 29
    end -= 1
    bins = [0]
    bins.extend(range(1 + (beg >> 26), 1 + (end >> 26) + 1))
    bins.extend(range(9 + (beg >> 23), 9 + (end >> 23) + 1))
    bins.extend(range(73 + (beg >> 20), 73 + (end >> 20) + 1))
    bins.extend(range(585 + (beg >> 17), 585 + (end >> 17) + 1))
    bins.extend(range(4681 + (beg >> 14), 4681 + (end >> 14) + 1))
    return bins


def chunks_for_region(ref: dict[str, Any], beg: int, end: int) -> list[tuple[int, int]]:
    linear = ref["linear"]
    li = beg >> 14
    min_off = int(linear[li]) if 0 <= li < len(linear) else 0
    out = []
    for bin_id in reg2bins(beg, end):
        for cb, ce in ref["bins"].get(bin_id, []):
            if int(ce) > min_off:
                out.append((int(cb), int(ce)))
    out.sort()
    return out


def merge_virtual_chunks(
    chunks: list[tuple[int, int]], max_compressed_gap: int = 65536
) -> list[tuple[int, int]]:
    if not chunks:
        return []
    merged: list[list[int]] = [[chunks[0][0], chunks[0][1]]]
    for b, e in chunks[1:]:
        cur = merged[-1]
        cur_end_comp = cur[1] >> 16
        b_comp = b >> 16
        if b <= cur[1] or b_comp <= cur_end_comp + max_compressed_gap:
            if e > cur[1]:
                cur[1] = e
        else:
            merged.append([b, e])
    return [(int(a), int(b)) for a, b in merged]


# ----------------------------- BGZF reader -----------------------------

def bgzf_block_size(buf: bytes, offset: int) -> int:
    if offset + 12 > len(buf):
        raise ValueError("truncated BGZF header")
    if buf[offset:offset+2] != b"\x1f\x8b":
        raise ValueError("BGZF block does not begin with gzip magic")
    flg = buf[offset + 3]
    if not (flg & 4):
        raise ValueError("gzip member lacks FEXTRA; not BGZF")
    xlen = struct.unpack("<H", buf[offset+10:offset+12])[0]
    extra_start = offset + 12
    extra_end = extra_start + xlen
    if extra_end > len(buf):
        raise ValueError("truncated BGZF extra field")
    p = extra_start
    while p + 4 <= extra_end:
        si1, si2 = buf[p], buf[p+1]
        slen = struct.unpack("<H", buf[p+2:p+4])[0]
        p += 4
        if p + slen > extra_end:
            raise ValueError("malformed BGZF extra subfield")
        if si1 == 66 and si2 == 67 and slen == 2:  # 'B','C'
            bsize = struct.unpack("<H", buf[p:p+2])[0]
            return int(bsize) + 1
        p += slen
    raise ValueError("BGZF BC/BSIZE subfield not found")


def decompress_virtual_chunk(
    payload: bytes,
    payload_abs_start: int,
    vbeg: int,
    vend: int,
) -> bytes:
    cstart, ustart = int(vbeg >> 16), int(vbeg & 0xFFFF)
    cend, uend = int(vend >> 16), int(vend & 0xFFFF)
    if payload_abs_start != cstart:
        raise ValueError("payload does not start at chunk's compressed block offset")

    pieces: list[bytes] = []
    p = 0
    abs_off = cstart
    while p < len(payload):
        if abs_off > cend or (abs_off == cend and uend == 0):
            break
        total = bgzf_block_size(payload, p)
        if p + total > len(payload):
            raise ValueError(
                f"range ended before complete BGZF block at compressed offset {abs_off}"
            )
        block = payload[p:p+total]
        dec = gzip.decompress(block)

        lo = ustart if abs_off == cstart else 0
        if lo > len(dec):
            raise ValueError("virtual start offset exceeds decompressed BGZF block")

        if abs_off == cend:
            hi = uend
            if hi > len(dec):
                raise ValueError("virtual end offset exceeds decompressed BGZF block")
            pieces.append(dec[lo:hi])
            break
        pieces.append(dec[lo:])

        p += total
        abs_off += total

    return b"".join(pieces)


# --------------------------- VCF interpretation ---------------------------

def parse_gt(format_field: str, sample_field: str) -> list[int | None]:
    fmt = format_field.split(":")
    vals = sample_field.split(":")
    try:
        idx = fmt.index("GT")
    except ValueError:
        return []
    if idx >= len(vals):
        return []
    gt = vals[idx]
    if gt in ("", ".", "./.", ".|."):
        return [None, None]
    parts = re.split(r"[\/|]", gt)
    out: list[int | None] = []
    for x in parts:
        if x in ("", "."):
            out.append(None)
        else:
            try:
                out.append(int(x))
            except ValueError:
                out.append(None)
    return out


def called_bases(record: dict[str, Any] | None) -> tuple[list[str], bool]:
    if record is None:
        return [], False
    gt = record.get("gt") or []
    alleles = [str(record.get("ref", "")).upper()] + [
        str(x).upper() for x in record.get("alts", [])
    ]
    bases: list[str] = []
    fully_called = len(gt) > 0
    for a in gt:
        if a is None or a < 0 or a >= len(alleles):
            fully_called = False
            continue
        b = alleles[a]
        if len(b) != 1 or b not in "ACGT":
            fully_called = False
            continue
        bases.append(b)
    if len(bases) != len(gt):
        fully_called = False
    return bases, fully_called


def target_state_from_calls(
    calls: dict[int, dict[str, Any]],
    start: int,
    end: int,
    strand: str,
    guide_sequence_key: str,
) -> tuple[dict[str, Any], int]:
    positions = list(range(int(start), int(end) + 1))
    if len(positions) != 23:
        raise ValueError("target is not 23 bp")

    expected_genomic = expected_genomic_23(guide_sequence_key, strand)
    base_called: list[bool] = []
    mismatch_by_pos: list[int] = []
    called_by_pos: list[int] = []
    present_records = 0

    for off, pos in enumerate(positions):
        rec = calls.get(pos)
        if rec is None:
            base_called.append(False)
            mismatch_by_pos.append(0)
            called_by_pos.append(0)
            continue
        present_records += 1
        bases, fully_called = called_bases(rec)
        base_called.append(fully_called)
        called_by_pos.append(len(bases))
        mismatch_by_pos.append(sum(b != expected_genomic[off] for b in bases))

    if strand == "+":
        protospacer_idx = list(range(0, 20))
        pam_idx = [20, 21, 22]
        pam_critical_idx = [21, 22]
    elif strand == "-":
        protospacer_idx = list(range(3, 23))
        pam_idx = [0, 1, 2]
        pam_critical_idx = [0, 1]
    else:
        raise ValueError(f"unexpected guide strand {strand!r}")

    callable23 = all(base_called)
    protospacer_callable = all(base_called[i] for i in protospacer_idx)
    pam_callable = all(base_called[i] for i in pam_idx)

    strict_exact23 = bool(callable23 and sum(mismatch_by_pos) == 0)
    protospacer_exact = bool(
        protospacer_callable and sum(mismatch_by_pos[i] for i in protospacer_idx) == 0
    )
    pam_exact = bool(pam_callable and sum(mismatch_by_pos[i] for i in pam_idx) == 0)

    # Functional SpCas9 PAM = NGG. The first PAM base is intentionally allowed
    # to vary; only the two guide-oriented G positions are required to remain G.
    pam_ngg_intact = bool(
        pam_callable and all(mismatch_by_pos[i] == 0 for i in pam_critical_idx)
    )
    functional_target_intact = bool(protospacer_exact and pam_ngg_intact)

    called_alleles = int(sum(called_by_pos))
    mismatch_alleles = int(sum(mismatch_by_pos))
    state = {
        "callable_23bp": callable23,
        "exact_23bp": strict_exact23 if callable23 else np.nan,  # compatibility alias
        "strict_genotype_exact_23bp": strict_exact23 if callable23 else np.nan,
        "protospacer_callable_20bp": protospacer_callable,
        "protospacer_exact_20bp": protospacer_exact if protospacer_callable else np.nan,
        "pam_callable": pam_callable,
        "pam_intact": pam_ngg_intact if pam_callable else np.nan,  # compatibility alias
        "pam_ngg_intact": pam_ngg_intact if pam_callable else np.nan,
        "pam_exact_3bp": pam_exact if pam_callable else np.nan,
        "functional_target_intact": functional_target_intact if (protospacer_callable and pam_callable) else np.nan,
        "called_alleles_23bp": called_alleles,
        "nonref_alleles_23bp": mismatch_alleles,  # valid because reference audit is enforced
        "mismatch_alleles_23bp": mismatch_alleles,
        "called_alleles_by_position_23bp": ";".join(map(str, called_by_pos)),
        "mismatch_alleles_by_position_23bp": ";".join(map(str, mismatch_by_pos)),
        "vcf_records_present_23bp": int(present_records),
    }
    return state, present_records


def _sum_vectors(series: pd.Series, width: int = 23) -> np.ndarray:
    arr = []
    for x in series.astype(str):
        v = [int(z) for z in x.split(";")]
        if len(v) != width:
            raise ValueError(f"expected {width}-element positional vector, got {len(v)}")
        arr.append(v)
    if not arr:
        return np.zeros(width, dtype=int)
    return np.asarray(arr, dtype=int).sum(axis=0)


def _aggregate_one(d: pd.DataFrame) -> dict[str, Any]:
    called = d[d.callable_23bp == True]
    proto = d[d.protospacer_callable_20bp == True]
    pamc = d[d.pam_callable == True]
    functional_called = d[
        (d.protospacer_callable_20bp == True) & (d.pam_callable == True)
    ]
    ca_pos = _sum_vectors(d.called_alleles_by_position_23bp)
    mm_pos = _sum_vectors(d.mismatch_alleles_by_position_23bp)
    pos_af = np.divide(
        mm_pos, ca_pos, out=np.full(23, np.nan, dtype=float), where=ca_pos > 0
    )
    ca = int(ca_pos.sum())
    mm = int(mm_pos.sum())
    exact_n = int(called["strict_genotype_exact_23bp"].astype("boolean").fillna(False).sum())
    pam_ngg_n = int(pamc["pam_ngg_intact"].astype("boolean").fillna(False).sum())
    return {
        "n_samples_total": int(len(d)),
        "n_samples_called_23bp": int(len(called)),
        "n_samples_exact_23bp": exact_n,  # compatibility alias
        "n_samples_strict_genotype_exact_23bp": exact_n,
        "n_samples_protospacer_called_20bp": int(len(proto)),
        "n_samples_protospacer_exact_20bp": int(proto["protospacer_exact_20bp"].astype("boolean").fillna(False).sum()),
        "n_samples_pam_called": int(len(pamc)),
        "n_samples_pam_intact": pam_ngg_n,  # compatibility alias
        "n_samples_pam_ngg_intact": pam_ngg_n,
        "n_samples_pam_exact_3bp": int(pamc["pam_exact_3bp"].astype("boolean").fillna(False).sum()),
        "n_samples_functional_target_called": int(len(functional_called)),
        "n_samples_functional_target_intact": int(functional_called["functional_target_intact"].astype("boolean").fillna(False).sum()),
        "callable_sample_fraction": float(len(called)/len(d)) if len(d) else np.nan,
        "exact_23bp_fraction": float(called.strict_genotype_exact_23bp.astype(float).mean()) if len(called) else np.nan,
        "strict_genotype_exact_23bp_fraction": float(called.strict_genotype_exact_23bp.astype(float).mean()) if len(called) else np.nan,
        "protospacer_exact_20bp_fraction": float(proto.protospacer_exact_20bp.astype(float).mean()) if len(proto) else np.nan,
        "pam_intact_fraction": float(pamc.pam_ngg_intact.astype(float).mean()) if len(pamc) else np.nan,
        "pam_ngg_intact_fraction": float(pamc.pam_ngg_intact.astype(float).mean()) if len(pamc) else np.nan,
        "pam_exact_3bp_fraction": float(pamc.pam_exact_3bp.astype(float).mean()) if len(pamc) else np.nan,
        "functional_target_intact_fraction": float(functional_called.functional_target_intact.astype(float).mean()) if len(functional_called) else np.nan,
        "called_alleles_23bp": ca,
        "nonreference_alleles_23bp": mm,  # compatibility alias after reference audit
        "mismatch_alleles_23bp": mm,
        "nonreference_allele_fraction_23bp": float(mm/ca) if ca else np.nan,
        "max_position_nonreference_allele_fraction_23bp": float(np.nanmax(pos_af)) if np.isfinite(pos_af).any() else np.nan,
        "mean_vcf_records_present_23bp": float(d.vcf_records_present_23bp.mean()),
    }


def aggregate(states: pd.DataFrame, targets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    tax_rows = []
    for sid, d in states.groupby("population_site_id", sort=False):
        rows.append({"population_site_id": sid, **_aggregate_one(d)})
    for (sid, tax), d in states.groupby(["population_site_id", "taxon"], sort=False):
        tax_rows.append({"population_site_id": sid, "taxon": tax, **_aggregate_one(d)})
    ann = targets[["population_site_id", "gene_id", "selection_reason"]].drop_duplicates()
    return (
        pd.DataFrame(rows).merge(ann, on="population_site_id", how="left"),
        pd.DataFrame(tax_rows).merge(ann, on="population_site_id", how="left"),
    )


def build_target_maps(
    targets: pd.DataFrame,
) -> tuple[dict[str, dict[int, list[str]]], dict[str, dict[int, str]]]:
    out: dict[str, dict[int, list[str]]] = {}
    expected: dict[str, dict[int, str]] = {}
    for contig in CANONICAL:
        d = targets[targets.genomic_seqid.astype(str) == contig]
        if d.empty:
            continue
        m: dict[int, list[str]] = defaultdict(list)
        e: dict[int, str] = {}
        for r in d.itertuples(index=False):
            seq = expected_genomic_23(str(r.guide_sequence_key), str(r.guide_genomic_strand))
            for off, p in enumerate(range(int(r.genomic_start), int(r.genomic_end) + 1)):
                m[p].append(str(r.population_site_id))
                b = seq[off]
                if p in e and e[p] != b:
                    raise ValueError(
                        f"frozen targets disagree on AgamP4 reference base at {contig}:{p}: "
                        f"{e[p]} vs {b}"
                    )
                e[p] = b
        out[contig] = dict(m)
        expected[contig] = e
    return out, expected



def extract_sample_states(
    session: requests.Session,
    sample_id: str,
    taxon: str,
    targets: pd.DataFrame,
    target_maps: dict[str, dict[int, list[str]]],
    expected_refs: dict[str, dict[int, str]],
    lock_hash: str,
    lock_file: Path,
    merge_gap: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    vu = vcf_url(sample_id)
    iu = tbi_url(sample_id)
    counter = ByteCounter()
    file_size, cr = strict_range_probe(session, vu)

    index_blob = download_small_index(session, iu, counter)
    refs = parse_tbi(index_blob)
    missing_refs = [c for c in target_maps if c not in refs]
    if missing_refs:
        raise RuntimeError(f"TBI missing required contigs: {missing_refs}")

    # Gather and merge all chunks required by the locked coordinate set.
    chunks_by_contig: dict[str, list[tuple[int, int]]] = {}
    query_chunk_total = 0
    for contig, posmap in target_maps.items():
        all_chunks: list[tuple[int, int]] = []
        # Use merged locked 23-bp intervals rather than one query per position.
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
            # TBI internal coordinates are 0-based, half-open.
            all_chunks.extend(chunks_for_region(refs[contig], start1 - 1, end1))
        query_chunk_total += len(all_chunks)
        chunks_by_contig[contig] = merge_virtual_chunks(all_chunks, max_compressed_gap=merge_gap)

    calls_by_contig: dict[str, dict[int, dict[str, Any]]] = {
        c: {} for c in target_maps
    }
    fetched_ranges = 0

    for contig, chunks in chunks_by_contig.items():
        needed_positions = target_maps[contig]
        for vbeg, vend in chunks:
            if sha256_file(lock_file) != lock_hash:
                raise RuntimeError("external target-lock SHA-256 changed during Ag3 acquisition")
            cstart = int(vbeg >> 16)
            cend = int(vend >> 16)
            # One extra maximum-size BGZF block ensures the terminal block is complete.
            byte_end = min(file_size - 1, cend + 65535)
            payload = fetch_range(session, vu, cstart, byte_end, file_size, counter)
            fetched_ranges += 1
            text = decompress_virtual_chunk(payload, cstart, vbeg, vend)
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
                ref = fields[3].upper()
                alts = [] if fields[4] in ("", ".") else [x.upper() for x in fields[4].split(",")]
                expected_ref = expected_refs[contig].get(pos)
                if expected_ref is not None:
                    if len(ref) != 1 or ref not in "ACGT":
                        raise RuntimeError(
                            f"{sample_id} {contig}:{pos} has non-SNP REF={ref!r} in Ag3 all-sites SNP VCF"
                        )
                    if ref != expected_ref:
                        raise RuntimeError(
                            f"Frozen guide/reference audit failed at {contig}:{pos}: "
                            f"Ag3 AgamP4 REF={ref}, frozen expected genomic base={expected_ref}"
                        )
                rec = {"ref": ref, "alts": alts, "gt": parse_gt(fields[8], fields[9])}
                if pos in calls_by_contig[contig]:
                    if calls_by_contig[contig][pos] != rec:
                        raise RuntimeError(
                            f"{sample_id} {contig}:{pos} has conflicting duplicate VCF records"
                        )
                else:
                    calls_by_contig[contig][pos] = rec

    states: list[dict[str, Any]] = []
    record_fracs = []
    for r in targets.itertuples(index=False):
        contig = str(r.genomic_seqid)
        calls = calls_by_contig.get(contig, {})
        st, n_present = target_state_from_calls(
            calls,
            int(r.genomic_start),
            int(r.genomic_end),
            str(r.guide_genomic_strand),
            str(r.guide_sequence_key),
        )
        record_fracs.append(n_present / 23.0)
        states.append({
            "sample_id": sample_id,
            "taxon": taxon,
            "population_site_id": str(r.population_site_id),
            "contig": contig,
            **st,
        })

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
        "mean_locked_coordinate_record_fraction": float(np.mean(record_fracs)),
        "min_locked_target_record_fraction": float(np.min(record_fracs)),
        "targets_with_all_23_records": int(sum(x == 1.0 for x in record_fracs)),
        "targets_total": int(len(record_fracs)),
        "frozen_reference_audit_pass": True,
        "frozen_reference_positions_checked": int(sum(len(x) for x in expected_refs.values())),
    }
    return states, qc



def checkpoint_paths(checkpoint_dir: Path, sample_id: str) -> tuple[Path, Path]:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(sample_id))
    return (
        checkpoint_dir / f"{safe}.states.tsv",
        checkpoint_dir / f"{safe}.qc.json",
    )


def load_sample_checkpoint(
    checkpoint_dir: Path,
    sample_id: str,
    taxon: str,
    targets: pd.DataFrame,
    lock_hash: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]] | None:
    """Load one completed sample checkpoint after strict integrity checks."""
    sp, qp = checkpoint_paths(checkpoint_dir, sample_id)
    if not sp.exists() or not qp.exists():
        return None

    try:
        st = pd.read_csv(sp, sep="\t")
        qc = json.loads(qp.read_text(encoding="utf-8"))
    except Exception as exc:
        log(
            f"WARNING: ignoring unreadable checkpoint for {sample_id}: {exc}"
        )
        return None

    required_state = {
        "sample_id",
        "taxon",
        "population_site_id",
        "callable_23bp",
        "strict_genotype_exact_23bp",
        "protospacer_callable_20bp",
        "protospacer_exact_20bp",
        "pam_callable",
        "pam_ngg_intact",
        "pam_exact_3bp",
        "functional_target_intact",
        "called_alleles_by_position_23bp",
        "mismatch_alleles_by_position_23bp",
        "vcf_records_present_23bp",
    }
    if not required_state.issubset(st.columns):
        missing_cols = sorted(required_state - set(st.columns))
        log(
            f"WARNING: ignoring incomplete checkpoint for {sample_id}: "
            f"missing columns={missing_cols}"
        )
        return None

    if str(qc.get("schema_version", "")) != SCHEMA_VERSION:
        log(
            f"WARNING: ignoring checkpoint for {sample_id}: schema mismatch"
        )
        return None
    if str(qc.get("target_lock_sha256", "")) != str(lock_hash):
        log(
            f"WARNING: ignoring checkpoint for {sample_id}: target-lock hash mismatch"
        )
        return None

    if len(st) != len(targets):
        log(
            f"WARNING: ignoring checkpoint for {sample_id}: "
            f"{len(st)} rows, expected {len(targets)}"
        )
        return None
    if st["sample_id"].astype(str).nunique() != 1 or str(st["sample_id"].iloc[0]) != str(sample_id):
        log(f"WARNING: ignoring checkpoint for {sample_id}: sample-ID mismatch")
        return None
    if st["taxon"].astype(str).nunique() != 1 or str(st["taxon"].iloc[0]) != str(taxon):
        log(f"WARNING: ignoring checkpoint for {sample_id}: taxon mismatch")
        return None

    got_sites = set(st["population_site_id"].astype(str))
    exp_sites = set(targets["population_site_id"].astype(str))
    if got_sites != exp_sites:
        log(
            f"WARNING: ignoring checkpoint for {sample_id}: locked-site set mismatch"
        )
        return None

    # Every validated checkpoint must contain all 23 VCF records for every
    # locked target. This is stronger than relying on the QC JSON alone.
    vcf_present = pd.to_numeric(
        st["vcf_records_present_23bp"], errors="coerce"
    )
    if vcf_present.isna().any() or not (vcf_present.astype(int) == 23).all():
        log(
            f"WARNING: ignoring checkpoint for {sample_id}: "
            "one or more locked targets do not contain all 23 VCF records"
        )
        return None

    # Positional allele vectors are required for the later site-level
    # max-position mismatch-frequency calculation. Validate their width now.
    for col in (
        "called_alleles_by_position_23bp",
        "mismatch_alleles_by_position_23bp",
    ):
        widths = st[col].astype(str).map(lambda x: len(x.split(";")))
        if not (widths == 23).all():
            log(
                f"WARNING: ignoring checkpoint for {sample_id}: "
                f"{col} does not contain 23 positions for every locked target"
            )
            return None

    if int(qc.get("targets_with_all_23_records", -1)) != len(targets):
        log(
            f"WARNING: ignoring checkpoint for {sample_id}: "
            "not all 23-bp target records were recovered"
        )
        return None

    return st.to_dict("records"), qc


def write_sample_checkpoint(
    checkpoint_dir: Path,
    sample_id: str,
    states: list[dict[str, Any]],
    qc: dict[str, Any],
    lock_hash: str,
) -> None:
    """Atomically persist a completed sample so interrupted runs can resume."""
    sp, qp = checkpoint_paths(checkpoint_dir, sample_id)
    qcout = dict(qc)
    qcout["schema_version"] = SCHEMA_VERSION
    qcout["target_lock_sha256"] = str(lock_hash)
    qcout["checkpoint_completed_at"] = now()
    atomic_df(pd.DataFrame(states), sp)
    atomic_json(qp, qcout)



def acquire_one_sample_worker(
    *,
    idx: int,
    n_selected: int,
    sid: str,
    tax: str,
    targets: pd.DataFrame,
    target_maps: dict[str, Any],
    expected_refs: dict[str, Any],
    lock_hash: str,
    lock_file: Path,
    merge_gap: int,
    checkpoint_dir: Path,
    smoke: bool,
    min_record_fraction: float,
) -> dict[str, Any]:
    """Acquire one sample with a private requests.Session and checkpoint it."""
    cached = load_sample_checkpoint(
        checkpoint_dir, sid, tax, targets, lock_hash
    )
    if cached is not None:
        states, qc = cached
        return {
            "idx": idx,
            "sample_id": sid,
            "taxon": tax,
            "states": states,
            "qc": qc,
            "resumed": True,
            "error": None,
        }

    if sha256_file(lock_file) != lock_hash:
        return {
            "idx": idx,
            "sample_id": sid,
            "taxon": tax,
            "states": None,
            "qc": None,
            "resumed": False,
            "error": "External target-lock SHA-256 changed during Ag3 acquisition",
        }

    session = requests.Session()
    session.headers.update(
        {"User-Agent": "MosqEditR-Ag3ExternalValidation/1.2-parallel-resume"}
    )
    try:
        states, qc = extract_sample_states(
            session, sid, tax, targets, target_maps, expected_refs,
            lock_hash, lock_file, merge_gap
        )

        if (
            smoke
            and qc["mean_locked_coordinate_record_fraction"]
            < min_record_fraction
        ):
            raise RuntimeError(
                f"Smoke sample {sid} recovered only "
                f"{qc['mean_locked_coordinate_record_fraction']:.3%} "
                f"of expected locked-coordinate records "
                f"(<{min_record_fraction:.1%})."
            )

        write_sample_checkpoint(
            checkpoint_dir, sid, states, qc, lock_hash
        )
        return {
            "idx": idx,
            "sample_id": sid,
            "taxon": tax,
            "states": states,
            "qc": qc,
            "resumed": False,
            "error": None,
        }
    except Exception as exc:
        return {
            "idx": idx,
            "sample_id": sid,
            "taxon": tax,
            "states": None,
            "qc": None,
            "resumed": False,
            "error": repr(exc),
        }
    finally:
        session.close()


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Ag3 locked-site external validation by pure-Python Tabix/BGZF HTTP ranges."
    )
    ap.add_argument("--project-root", default=None)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--smoke-test", action="store_true", help="1 gambiae + 1 coluzzii")
    mode.add_argument("--pilot", action="store_true", help="configured taxon-balanced pilot")
    mode.add_argument("--final", action="store_true", help="configured final cohort")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument(
        "--workers",
        type=int,
        default=4,
        help=(
            "Concurrent sample downloads for --final (default 4; allowed 1-6). "
            "Smoke/pilot are forced to 1 worker."
        ),
    )
    args = ap.parse_args()

    root = locate_root(args.project_root)
    os.chdir(root)
    cfg = load_cfg(root)
    ev = cfg.get("external_validation", {}) or {}
    if not bool(ev.get("enabled", True)):
        fail("external_validation.enabled is false")
    if bool(ev.get("require_no_discovery_retuning", True)) is not True:
        fail("external-validation isolation must remain TRUE")

    lock_file = root / "data_processed" / "05b_external_validation_targets.csv"
    if not lock_file.exists():
        fail("Missing data_processed/05b_external_validation_targets.csv")
    targets = pd.read_csv(lock_file)
    lock_hash = sha256_file(lock_file)
    required = {
        "population_site_id", "gene_id", "genomic_seqid", "genomic_start",
        "genomic_end", "guide_genomic_strand", "protospacer_20nt", "pam",
        "guide_sequence_key", "selection_locked_before_ag3", "selection_reason",
    }
    miss = sorted(required - set(targets.columns))
    if miss:
        fail(f"External target lock missing columns: {miss}")
    if targets.population_site_id.duplicated().any():
        fail("Duplicate site IDs in external target lock")
    if not targets.selection_locked_before_ag3.astype(bool).all():
        fail("External target lock is not marked immutable")
    if not targets.genomic_seqid.astype(str).isin(CANONICAL).all():
        fail("External target lock contains noncanonical contigs")
    if not ((targets.genomic_end.astype(int) - targets.genomic_start.astype(int) + 1) == 23).all():
        fail("Every external target must be exactly 23 bp")
    targets["protospacer_20nt"] = targets["protospacer_20nt"].astype(str).str.upper()
    targets["pam"] = targets["pam"].astype(str).str.upper()
    targets["guide_sequence_key"] = targets["guide_sequence_key"].astype(str).str.upper()
    if not (targets["protospacer_20nt"].str.fullmatch(r"[ACGT]{20}")).all():
        fail("Frozen external targets contain invalid protospacer_20nt sequence(s)")
    if not (targets["pam"].str.fullmatch(r"[ACGT]GG")).all():
        fail("Frozen external targets contain non-NGG PAM(s)")
    if not (
        targets["guide_sequence_key"]
        == targets["protospacer_20nt"] + targets["pam"]
    ).all():
        fail("Frozen guide_sequence_key != protospacer_20nt + pam for one or more targets")

    md = load_official_ag3_metadata(root)
    taxa = [str(x).lower() for x in ev.get("primary_taxa", list(FIXED_TAXA))]
    md = md[md.taxon.isin(taxa)].copy()
    if md.empty:
        fail("No configured primary taxa remain in official Ag3.0 metadata")

    final = bool(args.final)
    smoke = bool(args.smoke_test)

    if args.workers < 1 or args.workers > 6:
        fail("--workers must be between 1 and 6.")
    workers = int(args.workers) if final else 1
    if final:
        nmax = int(ev.get("final_max_samples_per_taxon", 0) or 0)
        selected = stratified_select(md, taxa, nmax)
        run_mode = "final"
    elif smoke:
        selected = stratified_select(md, taxa, 1)
        run_mode = "smoke"
    else:
        # Pilot remains a technical acquisition test, but now also uses official
        # MalariaGEN taxon calls and geographically stratified deterministic sampling.
        npilot = int(ev.get("pilot_samples_per_taxon", 25) or 25)
        selected = stratified_select(md, taxa, npilot)
        run_mode = "pilot"

    out_base = root / str(ev.get("output_dir", "data_raw/ag3_external"))
    if final:
        pilot_manifest = out_base / "pilot" / "run_manifest.json"
        pilot_ok = out_base / "pilot" / "ACQUISITION_COMPLETE.ok"
        if not pilot_ok.exists() or not pilot_manifest.exists():
            fail("--final is blocked until the 50-sample Ag3 pilot completes successfully.")
        try:
            pm = json.loads(pilot_manifest.read_text(encoding="utf-8"))
        except Exception as exc:
            fail(f"Could not read completed pilot manifest: {exc}")
        if int(pm.get("n_successful_samples", 0)) < 50 or int(pm.get("n_failed_samples", 999)) != 0:
            fail("--final requires the completed 50/50 technical pilot with zero sample failures.")
        if str(pm.get("target_lock_sha256", "")) != lock_hash:
            fail("Pilot target-lock SHA-256 differs from the current immutable external target lock.")

        smoke_manifest = out_base / "smoke" / "run_manifest.json"
        smoke_ok = out_base / "smoke" / "ACQUISITION_COMPLETE.ok"
        if not smoke_ok.exists() or not smoke_manifest.exists():
            fail(
                "--final requires a fresh --smoke-test with the sequence-aware validator "
                "after applying this patch."
            )
        try:
            sm = json.loads(smoke_manifest.read_text(encoding="utf-8"))
        except Exception as exc:
            fail(f"Could not read sequence-aware smoke-test manifest: {exc}")
        if str(sm.get("schema_version", "")) != SCHEMA_VERSION:
            fail(
                "--final requires a smoke test generated by the current sequence-aware "
                f"schema {SCHEMA_VERSION!r}; rerun --smoke-test."
            )
        if int(sm.get("n_successful_samples", 0)) != 2 or int(sm.get("n_failed_samples", 999)) != 0:
            fail("--final requires a successful 2/2 sequence-aware smoke test.")
        if str(sm.get("target_lock_sha256", "")) != lock_hash:
            fail("Smoke-test target-lock SHA-256 differs from the current immutable external target lock.")
    out_dir = out_base / run_mode
    out_dir.mkdir(parents=True, exist_ok=True)
    atomic_df(selected, out_dir / "selected_samples.tsv")
    atomic_df(targets, out_dir / "locked_targets.tsv")

    log(f"Project: {root}")
    log(
        f"External-validation lock: {len(targets):,} sites / "
        f"{targets.gene_id.nunique():,} genes; SHA256={lock_hash[:16]}..."
    )
    metadata_sources = sorted(
        set(md["metadata_source"].dropna().astype(str))
    ) if "metadata_source" in md.columns else ["unspecified"]
    log(
        f"Ag3.0 validation metadata: {len(md):,} primary-taxon samples; "
        f"selected {len(selected):,} for {run_mode}"
    )
    log("Metadata source: " + " | ".join(metadata_sources))
    for tx in taxa:
        log(
            f"  {tx}: available={int((md.taxon == tx).sum()):,}; "
            f"selected={int((selected.taxon == tx).sum()):,}"
        )
    log("Engine: pure-Python Tabix index + exact HTTP BGZF byte ranges; sequence-aware frozen-guide validation")
    log(
        "Sampling: explicit Ag3.0 gambiae/coluzzii metadata; "
        "country-stratified SHA-256 selection; discovery remains frozen"
    )
    log("Missing VCF records are never interpreted as reference; they remain uncallable")

    if args.dry_run:
        log("DRY RUN complete: manifests read; no VCF/index genomic data opened")
        return 0

    target_maps, expected_refs = build_target_maps(targets)
    merge_gap = int(ev.get("tabix_merge_compressed_gap_bytes", 65536) or 65536)
    max_fail = float(ev.get("max_sample_failure_fraction", 0.10) or 0.10)
    min_record_fraction = float(ev.get("minimum_smoke_coordinate_record_fraction", 0.95) or 0.95)

    sample_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    error_rows: list[dict[str, str]] = []

    checkpoint_dir = out_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    n_selected = len(selected)
    resumed = 0
    completed = 0

    log(
        f"Acquisition concurrency: workers={workers} "
        f"({'parallel sample acquisition' if workers > 1 else 'serial'})"
    )

    jobs = []
    for idx, sm in selected.reset_index(drop=True).iterrows():
        jobs.append(
            {
                "idx": int(idx),
                "sid": str(sm.sample_id),
                "tax": str(sm.taxon),
            }
        )

    def consume_result(result: dict[str, Any]) -> None:
        nonlocal resumed, completed

        idx = int(result["idx"])
        sid = str(result["sample_id"])
        tax = str(result["taxon"])

        if result["error"] is not None:
            error_rows.append(
                {
                    "sample_id": sid,
                    "taxon": tax,
                    "vcf_url": vcf_url(sid),
                    "error": str(result["error"]),
                }
            )
            log(
                f"WARNING: [{idx+1}/{n_selected}] sample failed: "
                f"{sid} ({tax}): {result['error']}"
            )
            atomic_df(
                pd.DataFrame(
                    error_rows,
                    columns=["sample_id", "taxon", "vcf_url", "error"],
                ),
                out_dir / "sample_errors.tsv",
            )
            return

        states_one = result["states"]
        qc = result["qc"]
        sample_rows.extend(states_one)
        qc_rows.append(qc)
        completed += 1

        if bool(result["resumed"]):
            resumed += 1
            log(
                f"[{idx+1}/{n_selected}] {sid} ({tax}) already complete; "
                "loaded validated checkpoint"
            )
        else:
            log(
                f"[{idx+1}/{n_selected}] {sid}: recovered all 23 VCF records "
                f"for {qc['targets_with_all_23_records']}/{qc['targets_total']} "
                f"targets; VCF bytes="
                f"{qc['vcf_bytes_downloaded']/1024/1024:.2f} MiB; "
                f"TBI={qc['tbi_bytes_downloaded']/1024:.1f} KiB; "
                f"ranges={qc['merged_vcf_ranges']}; checkpoint saved"
            )

        atomic_json(
            out_dir / "progress.json",
            {
                "schema_version": SCHEMA_VERSION,
                "target_lock_sha256": lock_hash,
                "run_mode": run_mode,
                "workers": int(workers),
                "n_selected_samples": int(n_selected),
                "n_completed_checkpoints": int(completed),
                "n_failed_this_run": int(len(error_rows)),
                "last_completed_sample": sid,
                "updated": now(),
            },
        )

    if workers == 1:
        for job in jobs:
            idx, sid, tax = job["idx"], job["sid"], job["tax"]

            if load_sample_checkpoint(
                checkpoint_dir, sid, tax, targets, lock_hash
            ) is None:
                log(
                    f"[{idx+1}/{n_selected}] Ag3 targeted VCF: "
                    f"{sid} ({tax})"
                )

            result = acquire_one_sample_worker(
                idx=idx,
                n_selected=n_selected,
                sid=sid,
                tax=tax,
                targets=targets,
                target_maps=target_maps,
                expected_refs=expected_refs,
                lock_hash=lock_hash,
                lock_file=lock_file,
                merge_gap=merge_gap,
                checkpoint_dir=checkpoint_dir,
                smoke=smoke,
                min_record_fraction=min_record_fraction,
            )
            consume_result(result)

            if smoke and result["error"] is not None:
                fail("Smoke test requires both samples to succeed.")

            if (
                len(error_rows) / max(1, idx + 1) > max_fail
                and (idx + 1) >= 10
            ):
                fail(
                    f"Sample failure fraction exceeded {max_fail:.0%}; stopping."
                )
    else:
        with ThreadPoolExecutor(
            max_workers=workers,
            thread_name_prefix="ag3sample",
        ) as pool:
            future_map = {}

            for job in jobs:
                idx, sid, tax = job["idx"], job["sid"], job["tax"]

                cached = load_sample_checkpoint(
                    checkpoint_dir, sid, tax, targets, lock_hash
                )
                if cached is not None:
                    states_one, qc = cached
                    consume_result(
                        {
                            "idx": idx,
                            "sample_id": sid,
                            "taxon": tax,
                            "states": states_one,
                            "qc": qc,
                            "resumed": True,
                            "error": None,
                        }
                    )
                    continue

                log(
                    f"[{idx+1}/{n_selected}] queued Ag3 targeted VCF: "
                    f"{sid} ({tax})"
                )

                fut = pool.submit(
                    acquire_one_sample_worker,
                    idx=idx,
                    n_selected=n_selected,
                    sid=sid,
                    tax=tax,
                    targets=targets,
                    target_maps=target_maps,
                    expected_refs=expected_refs,
                    lock_hash=lock_hash,
                    lock_file=lock_file,
                    merge_gap=merge_gap,
                    checkpoint_dir=checkpoint_dir,
                    smoke=False,
                    min_record_fraction=min_record_fraction,
                )
                future_map[fut] = (idx, sid, tax)

            for fut in as_completed(future_map):
                consume_result(fut.result())

        if len(error_rows) / max(1, n_selected) > max_fail:
            fail(
                f"Sample failure fraction exceeded {max_fail:.0%}; stopping."
            )

    if resumed:
        log(
            f"Resume summary: reused {resumed}/{n_selected} validated "
            "sample checkpoints; no network transfer repeated for those samples"
        )

    # The manuscript-level final cohort is locked at the selected 500 samples.
    # A transient network failure must never silently reduce that cohort.
    # If any final sample failed in this invocation, stop without writing the
    # completion marker. Re-running the same command will reuse all completed
    # checkpoints and retry only the samples without a valid checkpoint.
    if final and error_rows:
        failed_ids = ", ".join(x["sample_id"] for x in error_rows[:20])
        more = "" if len(error_rows) <= 20 else f" (+{len(error_rows)-20} more)"
        fail(
            f"Final acquisition has {len(error_rows)} failed sample(s): "
            f"{failed_ids}{more}. No completion marker written. "
            "Re-run the same --final command to retry only failed/uncompleted samples."
        )

    states = pd.DataFrame(sample_rows)
    qcdf = pd.DataFrame(qc_rows)
    errors = pd.DataFrame(
        error_rows, columns=["sample_id", "taxon", "vcf_url", "error"]
    )
    atomic_df(qcdf, out_dir / "sample_transfer_qc.tsv")
    atomic_df(errors, out_dir / "sample_errors.tsv")
    if states.empty:
        fail("No external-validation states were generated.")

    good_samples = states.sample_id.nunique()
    if final:
        if good_samples != n_selected:
            fail(
                f"Final cohort incomplete: {good_samples}/{n_selected} selected "
                "samples produced valid states. Re-run --final to resume."
            )
    elif good_samples < max(2 if smoke else 10, int(0.8 * n_selected)):
        fail(f"Only {good_samples}/{n_selected} selected samples produced states")
    expected_rows = good_samples * len(targets)
    if len(states) != expected_rows:
        counts = states.groupby("sample_id").population_site_id.nunique()
        fail(
            f"External state matrix incomplete: rows={len(states)}, expected={expected_rows}; "
            f"site count range={counts.min()}-{counts.max()}"
        )

    site_summary, taxon_summary = aggregate(states, targets)
    atomic_df(states, out_dir / "sample_site_states.tsv")
    atomic_df(site_summary, out_dir / "site_summary.tsv")
    atomic_df(taxon_summary, out_dir / "site_taxon_summary.tsv")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "run_mode": run_mode,
        "external_dataset": "Ag3.0 public Sanger per-sample all-sites VCF",
        "access_method": "local .tbi parse + HTTP BGZF byte ranges + local decompression",
        "target_lock_sha256": lock_hash,
        "n_locked_targets": int(len(targets)),
        "n_locked_genes": int(targets.gene_id.nunique()),
        "n_selected_samples": int(n_selected),
        "n_successful_samples": int(good_samples),
        "n_failed_samples": int(len(errors)),
        "vcf_bytes_downloaded_total": int(qcdf.vcf_bytes_downloaded.sum()) if len(qcdf) else 0,
        "tbi_bytes_downloaded_total": int(qcdf.tbi_bytes_downloaded.sum()) if len(qcdf) else 0,
        "http_requests_total": int(qcdf.http_requests.sum()) if len(qcdf) else 0,
        "validation_isolation": True,
        "sample_metadata_source": " | ".join(
            sorted(set(md["metadata_source"].dropna().astype(str)))
        ) if "metadata_source" in md.columns else "unspecified",
        "crosses_excluded": True,
        "final_selection_strategy": "country_stratified_sha256",
        "frozen_sequence_validation": "AgamP4 REF audited against frozen guide_sequence_key after strand normalization",
        "exact_23bp_metric_semantics": "strict unphased genotype: every called allele across the 23-bp frozen target matches the frozen sequence",
        "pam_intact_metric_semantics": "functional NGG: PAM N may vary; both guide-oriented G positions must remain G in every called allele",
        "whole_vcf_downloaded": False,
        "restart_safe_per_sample_checkpoints": True,
        "sample_concurrency_workers": int(workers),
        "parallelization_level": "sample",
        "range_read_timeout_seconds": 30,
        "transient_http_retries": 5,
        "created": now(),
    }
    atomic_json(out_dir / "run_manifest.json", manifest)
    atomic_text(out_dir / "ACQUISITION_COMPLETE.ok", "OK\n")
    log(
        f"Ag3 {run_mode} external acquisition completed: "
        f"{good_samples}/{n_selected} samples; "
        f"VCF transfer={manifest['vcf_bytes_downloaded_total']/1024/1024:.2f} MiB; "
        f"TBI transfer={manifest['tbi_bytes_downloaded_total']/1024/1024:.2f} MiB"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
