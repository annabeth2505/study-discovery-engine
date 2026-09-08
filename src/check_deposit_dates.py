"""
check_deposit_dates.py
Classify suspect PMIDs as REUSE/meta-analysis vs PRIMARY study.

Discriminator: deposit date vs paper year. If a paper's linked ENA accessions went
public YEARS BEFORE the paper, the paper reused existing data and the catalog link is a
misattribution. If the deposits are contemporaneous with the paper, it is a primary
study and the link is correct.

yes_catalog.tsv has no deposit-date column, so first_public is fetched from ENA. The
fetch is BULK (one OR'd query for all accessions, same approach as fetch_samples.py) --
result=study returns one row per study with first_public, so no per-study loop.

Verdicts, computed from the gap between paper year and deposit year:
  REUSE   -- EVERY accession predates the paper by >= GAP_YEARS (default 3). Even the
             newest deposit is old, so no data was generated for this paper.
  MIXED   -- some accessions are old, some contemporaneous. A primary study that ALSO
             reused older data -- reported separately rather than forced into a binary,
             because the catalog link is right for some accessions and wrong for others.
  PRIMARY -- deposits are contemporaneous with the paper.
  UNKNOWN -- ENA returned no date, or the catalog has no pub_year. Kept visible rather
             than counted as passing.

Output: ../results/deposit_date_check.tsv  (one row per PMID)

Usage:
  python check_deposit_dates.py
  python check_deposit_dates.py --pmids 38882494,39947133 --gap 3
"""

import io
import argparse
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
CATALOG = "../results/yes_catalog.tsv"
OUT = "../results/deposit_date_check.tsv"

GAP_YEARS = 3
BATCH = 100

SUSPECT_PMIDS = [38882494, 39947133, 37398234, 42282812, 33868800, 27388460,
                 35440042, 32473013, 39505912, 40149936, 33144315, 37208617,
                 39902954, 41432437]


def make_session():
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=0.5,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET", "POST"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


SESSION = make_session()


def fetch_first_public(accessions):
    """Bulk: one query per batch of accessions -> {accession: first_public}."""
    out = {}
    for i in range(0, len(accessions), BATCH):
        batch = accessions[i:i + BATCH]
        query = " OR ".join(f'study_accession="{a}"' for a in batch)
        r = SESSION.post(ENA_SEARCH,
                         data={"result": "study", "query": query,
                               "fields": "study_accession,first_public,last_updated",
                               "format": "tsv", "limit": 0},
                         timeout=180)
        if r.status_code == 200 and r.text.strip():
            df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
            for _, row in df.iterrows():
                acc = row["study_accession"]
                # an umbrella accession also returns its CHILD studies; keep only what
                # was asked for, so the dated/requested counts stay honest
                if acc in batch:
                    out[acc] = row.get("first_public")
    return out


def year_of(d):
    try:
        return int(str(d)[:4])
    except (TypeError, ValueError):
        return None


def run(pmids=None, gap_years=GAP_YEARS, out=OUT):
    pmids = pmids or SUSPECT_PMIDS
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    cat["pmid_i"] = pd.to_numeric(cat["pmid"], errors="coerce")
    sub = cat[cat.pmid_i.isin(pmids)]

    accs = sorted(sub.study_accession.dropna().astype(str).unique())
    print(f"{len(pmids)} PMIDs, {len(accs)} accessions -- fetching first_public from ENA...")
    dates = fetch_first_public(accs)
    print(f"  dated: {len(dates)}/{len(accs)}\n")

    rows = []
    for pmid in pmids:
        s = sub[sub.pmid_i == pmid]
        if s.empty:
            rows.append({"pmid": pmid, "verdict": "UNKNOWN",
                         "note": "no accession in yes_catalog"})
            continue

        paper_year = year_of(s.pub_year.dropna().iloc[0]) if s.pub_year.notna().any() else None
        acc_list = sorted(s.study_accession.dropna().astype(str).unique())
        years = {a: year_of(dates.get(a)) for a in acc_list}
        dated = {a: y for a, y in years.items() if y}

        if not dated or paper_year is None:
            verdict, gmin, gmax, ymin, ymax = "UNKNOWN", None, None, None, None
        else:
            ymin, ymax = min(dated.values()), max(dated.values())
            gmax = paper_year - ymin          # gap from the OLDEST deposit
            gmin = paper_year - ymax          # gap from the NEWEST deposit
            if gmin >= gap_years:
                verdict = "REUSE"
            elif gmax >= gap_years:
                verdict = "MIXED"
            else:
                verdict = "PRIMARY"

        rows.append({
            "pmid": pmid,
            "paper_year": paper_year,
            "paper_title": (s.paper_title.dropna().iloc[0][:70]
                            if s.paper_title.notna().any() else ""),
            "n_accessions": len(acc_list),
            "n_dated": len(dated),
            "deposit_year_min": ymin,
            "deposit_year_max": ymax,
            "gap_from_oldest": gmax,
            "gap_from_newest": gmin,
            "verdict": verdict,
            "accession_years": "; ".join(
                f"{a}:{years[a] or 'no-date'}" for a in acc_list),
        })

    res = pd.DataFrame(rows)
    res.to_csv(out, sep="\t", index=False)

    # ------------------------------ report ------------------------------
    W = 100
    print("=" * W)
    print(f"DEPOSIT DATE vs PAPER YEAR   (REUSE if every deposit predates paper by >={gap_years}y)")
    print("=" * W)
    order = {"REUSE": 0, "MIXED": 1, "PRIMARY": 2, "UNKNOWN": 3}
    for _, r in res.sort_values("verdict", key=lambda c: c.map(order)).iterrows():
        print(f"\nPMID {int(r.pmid)}  [{r.verdict}]   paper {r.paper_year}")
        print(f"  {r.paper_title}")
        if r.verdict == "UNKNOWN":
            print(f"  {r.get('note', 'no deposit date returned by ENA')}")
            continue
        print(f"  deposits {r.deposit_year_min}-{r.deposit_year_max}  "
              f"({r.n_dated}/{r.n_accessions} dated)   "
              f"gap: {r.gap_from_newest}y from newest, {r.gap_from_oldest}y from oldest")
        print(f"  {r.accession_years[:150]}")

    print("\n" + "=" * W)
    for k, v in res.verdict.value_counts().items():
        print(f"  {k:9s} {v}")
    print(f"\nSaved -> {out}")
    return res


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pmids", type=str, default=None)
    p.add_argument("--gap", type=int, default=GAP_YEARS)
    p.add_argument("--out", type=str, default=OUT)
    a = p.parse_args()
    run(pmids=[int(x) for x in a.pmids.split(",")] if a.pmids else None,
        gap_years=a.gap, out=a.out)
