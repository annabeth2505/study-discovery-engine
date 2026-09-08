"""
merge_ena_field_summary.py
Merge ena_field_summary.tsv into yes_catalog.tsv.

The summary was built as a standalone file on purpose; this is the deliberate merge.
Every incoming column except the summary itself is RENAMED with an ena_ prefix, because
their generic names collide or would be ambiguous next to 79 existing columns:

    n_samples -> n_ena_samples   COLLIDES with the catalog's own n_samples
    n_fields  -> n_ena_fields    generic; there are other *_fields notions in play
    n_runs    -> n_ena_runs      generic
    status    -> ena_fetch_status  far too generic to sit in a 79-column table

Backs up first. Idempotent: re-running refreshes these columns rather than duplicating.
"""

import shutil
import pandas as pd

R = "../results/"
CATALOG = R + "yes_catalog.tsv"
BACKUP = R + "yes_catalog.tsv.bak6"
SUMMARY = R + "ena_field_summary.tsv"

RENAME = {"n_fields": "n_ena_fields", "n_samples": "n_ena_samples",
          "n_runs": "n_ena_runs", "status": "ena_fetch_status"}
ADD = ["ena_field_summary"] + list(RENAME.values())


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    summ = pd.read_csv(SUMMARY, sep="\t").rename(columns=RENAME)
    n_before, cols_before = len(cat), len(cat.columns)

    shutil.copy(CATALOG, BACKUP)
    print(f"backup -> {BACKUP}")
    print(f"renamed on the way in: {RENAME}")

    already = [c for c in ADD if c in cat.columns]
    if already:
        print(f"refreshing {len(already)} column(s) already present")
        cat = cat.drop(columns=already)

    merged = cat.merge(summ[["study_accession"] + ADD], on="study_accession", how="left")

    bad = [c for c in merged.columns if c.endswith("_x") or c.endswith("_y")]
    assert not bad, f"merge created suffixed columns: {bad}"
    assert len(merged) == n_before, f"row count changed: {n_before} -> {len(merged)}"

    merged.to_csv(CATALOG, sep="\t", index=False)

    print(f"\nrows    {n_before} -> {len(merged)}  (unchanged)")
    print(f"columns {cols_before} -> {len(merged.columns)}  (+{len(merged.columns)-cols_before})")
    print(f"\nena_field_summary populated: "
          f"{int(merged.ena_field_summary.notna().sum())}/{len(merged)}")
    print(f"  NO_ENA_ROWS (umbrella projects): "
          f"{int((merged.ena_fetch_status == 'NO_ENA_ROWS').sum())}")
    print(f"  fields per study: median {int(merged.n_ena_fields.median())}  "
          f"max {int(merged.n_ena_fields.max())}")
    print(f"\nSaved -> {CATALOG}")
    return merged


if __name__ == "__main__":
    run()
