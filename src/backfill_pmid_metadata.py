"""
backfill_pmid_metadata.py
Fill paper metadata for catalog studies that HAVE a pmid but are missing
doi / title / authors / year, via one Entrez efetch pass. Also builds
fulltext_link from the best available identifier (DOI > PMID).

Writes back into yes_catalog.tsv. Checkpointed + resumable.
"""

import os
import json
import time
import pandas as pd
from Bio import Entrez

CATALOG = "../results/yes_catalog.tsv"
CHECKPOINT = "../results/pmid_metadata.json"


def fetch_pubmed_meta(pmid):
    """Return dict of doi/title/year/journal/first_author/last_author for a PMID."""
    try:
        h = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        rec = Entrez.read(h); h.close(); time.sleep(0.34)
        art = rec["PubmedArticle"][0]["MedlineCitation"]["Article"]

        # DOI can live in ELocationID or the ArticleIdList
        doi = None
        for eid in art.get("ELocationID", []):
            if eid.attributes.get("EIdType") == "doi":
                doi = str(eid)
        if not doi:
            ids = rec["PubmedArticle"][0].get("PubmedData", {}).get("ArticleIdList", [])
            for aid in ids:
                if aid.attributes.get("IdType") == "doi":
                    doi = str(aid)

        title = str(art.get("ArticleTitle", "")) or None
        try:
            year = art["Journal"]["JournalIssue"]["PubDate"].get("Year")
        except Exception:
            year = None
        journal = str(art["Journal"].get("Title", "")) or None

        authors = art.get("AuthorList", [])
        def name(a):
            return f"{a.get('LastName','')} {a.get('Initials','')}".strip()
        first = name(authors[0]) if authors else None
        last = name(authors[-1]) if len(authors) > 1 else None

        return {"doi": doi, "paper_title": title, "pub_year": year,
                "journal": journal, "first_author": first, "last_author": last}
    except Exception:
        return {}


def run():
    cat = pd.read_csv(CATALOG, sep="\t")

    # need metadata: has pmid, but missing doi OR title
    need = cat[cat["pmid"].notna() &
               (cat["doi"].isna() | cat["paper_title"].isna())]
    pmids = need["pmid"].apply(lambda x: str(int(float(x)))).unique().tolist()
    print(f"Studies with PMID but missing doi/title: {len(need)}")
    print(f"Unique PMIDs to resolve: {len(pmids)}")

    got = {}
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            got = json.load(f)
        print(f"Resuming — {len(got)} done")

    todo = [p for p in pmids if p not in got]
    print(f"Fetching {len(todo)}...\n")
    for i, pmid in enumerate(todo):
        if i % 25 == 0:
            print(f"  {i}/{len(todo)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(got, f)
        got[pmid] = fetch_pubmed_meta(pmid)
        time.sleep(0.2)
    with open(CHECKPOINT, "w") as f:
        json.dump(got, f)

    # ---- write back, only filling blanks (don't overwrite existing) ----
    def pmid_key(row):
        return str(int(float(row["pmid"]))) if pd.notna(row["pmid"]) else None

    for field in ["doi", "paper_title", "pub_year", "journal", "first_author", "last_author"]:
        def fill(row, f=field):
            cur = row.get(f)
            if pd.notna(cur) and str(cur).strip():
                return cur
            k = pmid_key(row)
            return got.get(k, {}).get(f, cur) if k else cur
        cat[field] = cat.apply(fill, axis=1)

    # ---- build fulltext_link: DOI preferred, else PubMed ----
    def link(row):
        if pd.notna(row.get("doi")) and str(row["doi"]).strip():
            return f"https://doi.org/{row['doi']}"
        if pd.notna(row.get("pmid")):
            return f"https://pubmed.ncbi.nlm.nih.gov/{str(int(float(row['pmid'])))}"
        return None
    cat["fulltext_link"] = cat.apply(link, axis=1)

    cat.to_csv(CATALOG, sep="\t", index=False)

    n = len(cat)
    print(f"\n{'='*50}\nPMID BACKFILL COMPLETE\n{'='*50}")
    for f in ["doi", "paper_title", "pub_year", "journal", "first_author", "fulltext_link"]:
        filled = cat[f].notna().sum()
        print(f"  {f:15s}: {filled:>4}/{n} ({100*filled/n:.0f}%)")
    print(f"Saved -> {CATALOG}")
    return cat


if __name__ == "__main__":
    run()