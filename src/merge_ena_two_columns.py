"""
merge_ena_two_columns.py
Merge the two-column ENA field summary into yes_catalog.tsv, REPLACING the earlier
single-column version.

Superseded columns, dropped here (with a backup first):
    ena_field_summary   the Portal-only build. Its content is now ena_standard_fields,
                        which uses the same source but the safer " || " field separator
                        AND is paired with the custom column the Portal API could not see.
    n_ena_fields        counted fields for that superseded column.
    ena_fetch_status    status of that superseded fetch.
    n_ena_samples       duplicates n_total_samples / n_metagenomic_samples.
    n_ena_runs          duplicates the run counts already in samples.tsv.
Leaving them in place would mean two columns claiming to be "the ENA fields", one of them
stale and missing every custom tag -- exactly the ambiguity that makes an 84-column table
hard to trust.

Incoming counts are renamed with an ena_ prefix, since "n_standard_fields" alone is far
too generic to sit beside 80 other columns:
    n_standard_fields -> n_ena_standard_fields
    n_custom_fields   -> n_ena_custom_fields

Backs up first. Idempotent: re-running refreshes rather than duplicating.
"""

import shutil
import pandas as pd

R = "../results/"
CATALOG = R + "yes_catalog.tsv"
BACKUP = R + "yes_catalog.tsv.bak7"
SUMMARY = R + "ena_field_summary.tsv"

SUPERSEDED = ["ena_field_summary", "n_ena_fields", "ena_fetch_status",
              "n_ena_samples", "n_ena_runs"]
RENAME = {"n_standard_fields": "n_ena_standard_fields",
          "n_custom_fields": "n_ena_custom_fields"}
ADD = ["ena_standard_fields", "ena_custom_fields"] + list(RENAME.values())


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    summ = pd.read_csv(SUMMARY, sep="\t").rename(columns=RENAME)
    n_before, cols_before = len(cat), len(cat.columns)

    shutil.copy(CATALOG, BACKUP)
    print(f"backup -> {BACKUP}")

    drop = [c for c in SUPERSEDED if c in cat.columns]
    if drop:
        print(f"dropping superseded ({len(drop)}): {drop}")
        cat = cat.drop(columns=drop)

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
    print(f"columns {cols_before} -> {len(merged.columns)}  "
          f"({-len(drop)} superseded, +{len(ADD)} new)")
    print(f"\nena_standard_fields populated: "
          f"{int(merged.ena_standard_fields.notna().sum())}/{len(merged)}"
          f"   median {int(merged.n_ena_standard_fields.median())} fields")
    print(f"ena_custom_fields populated  : "
          f"{int((merged.n_ena_custom_fields > 0).sum())}/{len(merged)}"
          f"   max {int(merged.n_ena_custom_fields.max())} fields")
    print(f"\nSaved -> {CATALOG}")
    return merged


if __name__ == "__main__":
    run()
