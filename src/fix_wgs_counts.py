"""
fix_wgs_counts.py
Make the sample-count columns mean what their names say.

The problem
-----------
n_wgs_samples was unusable for two DIFFERENT reasons that looked identical in the sheet:

  (a) 181 studies: enrich_catalog_from_ena.py never ran on them, so the column sat at
      its DEFAULT of 0. That 0 meant "never measured", not "zero WGS samples". Once
      library_strategy_mix was backfilled from samples.tsv, those rows read
      "WGS:148 ... 0", which looks like a contradiction and is really a placeholder
      sitting next to a measurement.

  (b) 53 studies: the original formula was
          n_wgs = sum(v for k, v in strat.items() if k in ("WGS", "METAGENOMIC"))
      which tests library_STRATEGY for "METAGENOMIC" -- a library_SOURCE value that
      never appears as a strategy. So the "library_source beats library_strategy" rule
      never actually fired.

After this script
-----------------
  n_wgs_samples          RECOMPUTED: samples with library_strategy == WGS. Matches its
                         name exactly -- a strategy count.
  n_metagenomic_samples  samples with library_source == METAGENOMIC. The domain-rule
                         number, and the one to use for "how many real metagenomes".
  n_wgs_samples_original  the previous column, preserved verbatim so the before/after
                         stays auditable and nothing is silently rewritten.

Backs up first. Idempotent: n_wgs_samples_original is only captured once, so re-running
never overwrites the true original.
"""

import shutil
import pandas as pd

CATALOG = "../results/yes_catalog.tsv"
BACKUP = "../results/yes_catalog.tsv.bak3"
SAMPLES = "../results/samples.tsv"


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    s = pd.read_csv(SAMPLES, sep="\t", dtype=str)
    shutil.copy(CATALOG, BACKUP)
    print(f"backup -> {BACKUP}\n")

    # preserve the original ONCE, so re-running cannot destroy it
    if "n_wgs_samples_original" not in cat.columns:
        cat["n_wgs_samples_original"] = cat["n_wgs_samples"]
        print("captured n_wgs_samples -> n_wgs_samples_original")
    else:
        print("n_wgs_samples_original already present, left as is")

    wgs = s[s.library_strategy == "WGS"].groupby("study_accession").size()
    meta = s[s.library_source == "METAGENOMIC"].groupby("study_accession").size()

    cat["n_wgs_samples"] = cat.study_accession.map(wgs).fillna(0).astype(int)
    cat["n_metagenomic_samples"] = cat.study_accession.map(meta).fillna(0).astype(int)

    cat.to_csv(CATALOG, sep="\t", index=False)

    old, new = cat.n_wgs_samples_original, cat.n_wgs_samples
    changed = (old != new)
    print(f"\nrows changed: {int(changed.sum())} of {len(cat)}")
    print(f"  was 0, now > 0 : {int(((old == 0) & (new > 0)).sum())}")
    print(f"  still 0        : {int(((old == 0) & (new == 0)).sum())}")
    print()
    print(f"totals   n_wgs_samples_original : {int(old.sum()):,}")
    print(f"         n_wgs_samples (WGS)    : {int(new.sum()):,}")
    print(f"         n_metagenomic_samples  : {int(cat.n_metagenomic_samples.sum()):,}")
    print()
    print("no row should now show WGS in the mix with n_wgs_samples == 0:")
    bad = cat[(cat.n_wgs_samples == 0) &
              cat.library_strategy_mix.str.contains("WGS", na=False)]
    print(f"  contradictory rows remaining: {len(bad)}")
    print(f"\nSaved -> {CATALOG}")
    return cat


if __name__ == "__main__":
    run()
