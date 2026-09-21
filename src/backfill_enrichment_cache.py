"""
backfill_enrichment_cache.py
Fill the six ENA-enrichment cache columns for the 181 studies that
enrich_catalog_from_ena.py never ran on (they joined the catalog after its last run):

    n_host_species, host_species_list, body_sites, biome_labels,
    isolation_sources, sample_titles_cache

Computed from samples.tsv with the SAME definitions as enrich_catalog_from_ena.summarize()
(distinct non-empty values; sentinels like "missing" kept literally, as that script did),
so these rows are indistinguishable from ones the original script would have produced.

Why not just re-run enrich_catalog_from_ena.py: its write-back overwrites all nine of its
columns for EVERY study, including n_wgs_samples with the old formula that tests
library_strategy for "METAGENOMIC" (never matches). That would silently undo
fix_wgs_counts.py. This script touches only the six columns, only for the 181 studies,
and only cells that are blank. No ENA fetch.
"""

import shutil
import pandas as pd

R = "../results/"
CATALOG = R + "yes_catalog.tsv"
PRE_BACKFILL = R + "yes_catalog.tsv.bak2"   # state before any backfill: identifies the 181
SAMPLES = R + "samples.tsv"
BACKUP = R + "yes_catalog.tsv.bak9"

COLS = ["n_host_species", "host_species_list", "body_sites", "biome_labels",
        "isolation_sources", "sample_titles_cache"]


def distinct(sub, col):
    """Exactly enrich_catalog_from_ena.summarize().distinct()."""
    if col not in sub.columns:
        return []
    return sorted(set(str(x) for x in sub[col].dropna() if str(x).strip()))


def summarize(sub):
    hosts = distinct(sub, "host_scientific_name")
    return {
        "n_host_species": len(hosts),
        "host_species_list": "; ".join(hosts[:40]),
        "body_sites": "; ".join(distinct(sub, "host_body_site")),
        "biome_labels": "; ".join(distinct(sub, "scientific_name")[:20]),
        "isolation_sources": "; ".join(distinct(sub, "isolation_source")[:20]),
        "sample_titles_cache": "; ".join(distinct(sub, "sample_title")[:5]),
    }


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    before = pd.read_csv(PRE_BACKFILL, sep="\t", low_memory=False)
    target = set(before[before.library_strategy_mix.isna()].study_accession)
    s = pd.read_csv(SAMPLES, sep="\t", dtype=str)
    groups = dict(tuple(s[s.study_accession.isin(target)].groupby("study_accession")))
    print(f"target studies (never enriched): {len(target)}   with samples: {len(groups)}")

    shutil.copy(CATALOG, BACKUP)
    print(f"backup -> {BACKUP}")
    untouched_before = cat[~cat.study_accession.isin(target)][COLS].copy()

    filled = {c: 0 for c in COLS}
    for i, acc in cat.study_accession.items():
        if acc not in groups:
            continue
        vals = summarize(groups[acc])
        for c in COLS:
            cur = cat.at[i, c]
            blank = pd.isna(cur) or str(cur).strip() == "" or (c == "n_host_species" and cur == 0)
            if blank and vals[c] not in ("", None):
                cat.at[i, c] = vals[c]
                filled[c] += 1

    # the other 863 must be byte-for-byte unchanged in these columns
    after = cat[~cat.study_accession.isin(target)][COLS]
    assert untouched_before.fillna("<NA>").astype(str).equals(after.fillna("<NA>").astype(str)), \
        "a non-target study changed"

    cat.to_csv(CATALOG, sep="\t", index=False)
    t = cat[cat.study_accession.isin(target)]
    print("\ncells filled / still blank among the 181:")
    for c in COLS:
        still = int((t[c].isna() | (t[c].astype(str).str.strip() == "")).sum())
        print(f"  {c:22s} filled {filled[c]:3d}   still blank {still:3d}")
    print(f"\nother {len(cat) - len(target)} studies: unchanged (verified)")
    print(f"Saved -> {CATALOG}")


if __name__ == "__main__":
    run()
