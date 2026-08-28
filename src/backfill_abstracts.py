"""
backfill_abstracts.py
Fetch missing abstracts for catalog studies that HAVE a pmid but no abstract.
Pre-step before the full audit -- gives Claude more text to work with.

Reuses fetch_abstract_from_pmid from ena_fetcher. Checkpointed + resumable.
Writes the enriched abstracts back into yes_catalog.tsv (abstract column).
"""

import os
import json
import time
import pandas as pd
from src.ena_fetcher import fetch_abstract_from_pmid

CATALOG = "../results/yes_catalog.tsv"
CHECKPOINT = "../results/backfilled_abstracts.json"


def run():
    cat = pd.read_csv(CATALOG, sep="\t")

    # studies with a pmid but no abstract
    need = cat[cat["pmid"].notna() & (cat["abstract"].isna() | (cat["abstract"].astype(str).str.strip() == ""))]
    pmids = need["pmid"].apply(lambda x: str(int(float(x)))).unique().tolist()
    print(f"Studies with PMID but no abstract: {len(need)}")
    print(f"Unique PMIDs to fetch: {len(pmids)}")

    # resume
    got = {}
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            got = json.load(f)
        print(f"Resuming — {len(got)} already fetched")

    todo = [p for p in pmids if p not in got]
    print(f"Fetching {len(todo)}...\n")

    for i, pmid in enumerate(todo):
        if i % 25 == 0:
            print(f"  {i}/{len(todo)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(got, f)
        try:
            ab = fetch_abstract_from_pmid(pmid)
            got[pmid] = ab if ab else ""
        except Exception as e:
            got[pmid] = ""
        time.sleep(0.2)

    with open(CHECKPOINT, "w") as f:
        json.dump(got, f)

    # write back into catalog
    def fill(row):
        cur = row.get("abstract")
        if pd.notna(cur) and str(cur).strip():
            return cur
        if pd.notna(row.get("pmid")):
            return got.get(str(int(float(row["pmid"]))), cur)
        return cur

    cat["abstract"] = cat.apply(fill, axis=1)
    cat.to_csv(CATALOG, sep="\t", index=False)

    found = sum(1 for v in got.values() if v)
    print(f"\nAbstracts newly found: {found}/{len(pmids)}")
    print(f"Catalog abstract coverage now: {cat['abstract'].notna().sum()}/{len(cat)} "
          f"({100*cat['abstract'].notna().sum()/len(cat):.0f}%)")
    return cat


if __name__ == "__main__":
    run()