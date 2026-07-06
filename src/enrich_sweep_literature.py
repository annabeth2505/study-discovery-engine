"""
enrich_sweep_literature.py
Fill title/abstract/pmid for the sweep studies and finalize gut verdicts.

Reuses (already in src/ena_fetcher.py):
  - fetch_study_origin(accession)      -> ENA study title/description (+BioProject fallback)
  - fetch_abstract_from_pmid(pmid)     -> abstract text

Adds the ONE missing piece:
  - accession_to_pmid(accession)       -> reverse lookup via NCBI elink (bioproject->pubmed)

Adjustments baked in (from manual review of the classifier output):
  - 3 annelid studies moved DROP -> NEEDS_ABSTRACT (don't lose earthworm gut on a bare label)
  - 4 multi-host mega-surveys flagged (handle at sample level later, not as single-host gut)

Flow:
  1. load results/sweep_gut_classified.tsv
  2. apply adjustments
  3. for KEEP + NEEDS_ABSTRACT: get ENA title, then accession->PMID->abstract
  4. re-decide NEEDS_ABSTRACT using title+abstract text (title alone often settles it)
  5. write results/sweep_enriched.tsv with catalog-ready columns
"""

import time
import pandas as pd
from Bio import Entrez

# reuse existing functions + the gut regexes from the classifier
from src.ena_fetcher import fetch_study_origin, fetch_abstract_from_pmid
from src.classify_sweep_gut import GUT_RE, NOTGUT_RE

CLASSIFIED_PATH = "../results/sweep_gut_classified.tsv"
OUT_PATH = "../results/sweep_enriched.tsv"

# manual-review adjustments
ANNELID_TO_REVIEW = {"PRJNA365078", "PRJNA365079", "PRJNA365080"}
MULTIHOST_SURVEYS = {"PRJNA1200941", "PRJNA437674", "PRJEB91453", "PRJNA901878"}


# ---- accession -> PMID ----
def accession_to_pmid(accession):
    """Reverse lookup: study accession -> linked PMID via NCBI bioproject->pubmed elink."""
    try:
        h = Entrez.esearch(db="bioproject", term=accession)
        rec = Entrez.read(h); h.close()
        time.sleep(0.34)
        ids = rec.get("IdList", [])
        if not ids:
            return None
        uid = ids[0]
        h = Entrez.elink(dbfrom="bioproject", db="pubmed", id=uid)
        rec = Entrez.read(h); h.close()
        time.sleep(0.34)
        for ls in rec:
            for db in ls.get("LinkSetDb", []):
                links = db.get("Link", [])
                if links:
                    return links[0]["Id"]
        return None
    except Exception as e:
        print(f"  accession_to_pmid failed {accession}: {e}")
        return None


def gut_decision(text):
    """gut / not_gut / ambiguous from combined title+abstract text."""
    blob = str(text).lower()
    if GUT_RE.search(blob):
        return "gut"
    if NOTGUT_RE.search(blob):
        return "not_gut"
    return "ambiguous"


def run():
    df = pd.read_csv(CLASSIFIED_PATH, sep="\t")

    # --- adjustment 1: annelids DROP -> NEEDS_ABSTRACT ---
    df.loc[df["study_accession"].isin(ANNELID_TO_REVIEW), "gut_verdict"] = "NEEDS_ABSTRACT"

    # --- adjustment 2: flag mega-surveys ---
    df["is_multihost_survey"] = df["study_accession"].isin(MULTIHOST_SURVEYS)

    # only enrich studies we might keep
    to_enrich = df[df["gut_verdict"].isin(["KEEP", "NEEDS_ABSTRACT"])].copy()
    print(f"Enriching {len(to_enrich)} studies (skipping {len(df) - len(to_enrich)} DROP)...")

    titles, abstracts, pmids = {}, {}, {}
    for i, acc in enumerate(to_enrich["study_accession"]):
        if i % 10 == 0:
            print(f"  {i}/{len(to_enrich)}")
        # always-available: ENA study title
        try:
            origin = fetch_study_origin(acc)
            titles[acc] = origin.get("title") if isinstance(origin, dict) else origin
        except Exception:
            titles[acc] = None
        # reverse lookup -> abstract (many will have no paper; that's expected)
        pmid = accession_to_pmid(acc)
        pmids[acc] = pmid
        if pmid:
            try:
                abstracts[acc] = fetch_abstract_from_pmid(pmid)
            except Exception:
                abstracts[acc] = None
        time.sleep(0.2)

    df["paper_title"] = df["study_accession"].map(titles)
    df["abstract"] = df["study_accession"].map(abstracts)
    df["pmid"] = df["study_accession"].map(pmids)

    # --- re-decide NEEDS_ABSTRACT using title + abstract text ---
    def finalize(row):
        if row["gut_verdict"] == "KEEP":
            return "KEEP"
        if row["gut_verdict"] != "NEEDS_ABSTRACT":
            return row["gut_verdict"]
        text = f"{row.get('paper_title','')} {row.get('abstract','')}"
        d = gut_decision(text)
        if d == "gut":
            return "KEEP"
        if d == "not_gut":
            return "DROP"
        return "MANUAL_REVIEW"   # title+abstract still didn't settle it

    df["gut_verdict_final"] = df.apply(finalize, axis=1)

    df.to_csv(OUT_PATH, sep="\t", index=False)

    print(f"\n{'='*50}")
    print("AFTER LITERATURE ENRICHMENT")
    print(f"{'='*50}")
    print(df["gut_verdict_final"].value_counts().to_string())

    keep = df[df["gut_verdict_final"] == "KEEP"]
    print(f"\nCONFIRMED GUT (final): {len(keep)} studies")
    print(f"  single-host: {(~keep['is_multihost_survey']).sum()}")
    print(f"  multi-host surveys (handle at sample level): {keep['is_multihost_survey'].sum()}")
    print(f"\nConfirmed-gut by clade:")
    print(keep[~keep["is_multihost_survey"]]["clade"].value_counts().to_string())

    manual = df[df["gut_verdict_final"] == "MANUAL_REVIEW"]
    if len(manual):
        print(f"\nSTILL NEEDS MANUAL REVIEW ({len(manual)}) — no paper or unclear text:")
        print(manual[["study_accession", "clade", "paper_title"]].to_string(index=False))

    print(f"\nAbstract coverage on enriched: {df['abstract'].notna().sum()}/{len(to_enrich)}")
    print(f"(Low coverage is expected — most ENA studies are data-only deposits.)")
    print(f"\nSaved -> {OUT_PATH}")
    return df


if __name__ == "__main__":
    run()