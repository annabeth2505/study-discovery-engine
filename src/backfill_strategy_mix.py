"""
backfill_strategy_mix.py
Fill in library_strategy_mix for the studies that currently report nothing, and add a
CORRECTED metagenome sample count, using the sample-level catalog. No ENA refetch --
samples.tsv already holds per-sample library_strategy and library_source.

Why this is needed
------------------
234 studies report n_wgs_samples = 0, for two different reasons:

  (a) 181 have a BLANK library_strategy_mix. enrich_catalog_from_ena.py never ran on
      them (they entered the catalog afterwards), so they got 0 by default rather than
      by measurement. They hold 32,799 METAGENOMIC samples between them.

  (b) 53 have a real strategy mix with no "WGS" label (OTHER, WGA, Hi-C, WXS, Tn-Seq).
      The original count was:
          n_wgs = sum(v for k, v in strat.items() if k in ("WGS", "METAGENOMIC"))
      which tests library_STRATEGY for the value "METAGENOMIC" -- but METAGENOMIC is a
      library_SOURCE value and never appears as a strategy, so that half of the
      condition can never fire. The domain rule "library_source beats library_strategy"
      was therefore never actually implemented. 47 of these hold 1,996 METAGENOMIC
      samples.

What this writes
----------------
  library_strategy_mix   -- FILLED IN where blank. Existing values are left alone; a
                            fresh computation is compared against them and any
                            disagreement is REPORTED, not silently applied.
  library_source_mix     -- NEW. The authoritative axis per the domain rule.
  n_metagenomic_samples  -- NEW. Samples with library_source=METAGENOMIC, i.e. the count
                            n_wgs_samples was meant to be.

n_wgs_samples is deliberately NOT overwritten: it is wrong for 234 studies, and leaving
it in place beside the corrected column keeps the error auditable instead of erasing it.
Compare the two to see what changed.

Idempotent, and backs the catalog up first.
"""

import shutil
import pandas as pd

CATALOG = "../results/yes_catalog.tsv"
BACKUP = "../results/yes_catalog.tsv.bak2"
SAMPLES = "../results/samples.tsv"


def mix(series):
    vc = series.dropna().value_counts()
    return "; ".join(f"{k}:{v}" for k, v in vc.items())


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    s = pd.read_csv(SAMPLES, sep="\t", dtype=str)
    shutil.copy(CATALOG, BACKUP)
    print(f"backup -> {BACKUP}\n")

    g = s.groupby("study_accession")
    strat = g["library_strategy"].apply(mix)
    source = g["library_source"].apply(mix)
    n_meta = (s[s.library_source == "METAGENOMIC"].groupby("study_accession").size())

    fresh_strat = cat.study_accession.map(strat)
    was_blank = cat.library_strategy_mix.isna()

    # --- report disagreements on studies that ALREADY had a value (do not overwrite) ---
    both = cat.library_strategy_mix.notna() & fresh_strat.notna()
    differ = both & (cat.library_strategy_mix != fresh_strat)
    print(f"studies with an existing library_strategy_mix : {int(both.sum())}")
    print(f"  ...that disagree with a fresh computation   : {int(differ.sum())}"
          "   (left unchanged, listed below)")
    for _, r in cat[differ].head(8).iterrows():
        print(f"    {r.study_accession}: catalog={r.library_strategy_mix!r} "
              f"fresh={fresh_strat[r.name]!r}")

    # --- fill only the blanks ---
    cat.loc[was_blank, "library_strategy_mix"] = fresh_strat[was_blank]
    filled = int((was_blank & cat.library_strategy_mix.notna()).sum())

    # --- new columns ---
    cat["library_source_mix"] = cat.study_accession.map(source)
    cat["n_metagenomic_samples"] = cat.study_accession.map(n_meta).fillna(0).astype(int)

    cat.to_csv(CATALOG, sep="\t", index=False)

    still = int(cat.library_strategy_mix.isna().sum())
    print(f"\nlibrary_strategy_mix filled : {filled}")
    print(f"  still blank               : {still}")
    print(f"columns                     : {len(cat.columns)}")
    print()
    print("n_wgs_samples (original, kept) vs n_metagenomic_samples (corrected):")
    old_zero = cat[cat.n_wgs_samples == 0]
    print(f"  studies where old = 0     : {len(old_zero)}")
    print(f"  ...now counted > 0        : {int((old_zero.n_metagenomic_samples > 0).sum())}")
    print(f"  samples recovered         : {int(old_zero.n_metagenomic_samples.sum()):,}")
    print(f"\n  total n_wgs_samples       : {int(cat.n_wgs_samples.sum()):,}")
    print(f"  total n_metagenomic_samples: {int(cat.n_metagenomic_samples.sum()):,}")
    print(f"\nSaved -> {CATALOG}")
    return cat


if __name__ == "__main__":
    run()
