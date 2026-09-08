"""
fix_paper_host_layers.py
Correct the paper-derived host layers for reuse-linkage contamination, then recompute
multi-host and human-by-biome-label on the studies whose paper link we can trust.

The contamination
-----------------
Some accessions are linked to a paper that REUSED their data rather than to their own
primary paper. Everything extracted from that paper -- host_title_value,
host_abstract_value, is_model_organism -- then describes the WRONG study.

Proven case, PMID 27388460 (S24-7, "homeothermic animals"), fan-out of 5:
    PRJEB6456    real title: Dynamics ... Human Gut Microbiome during the First Year of Life
    PRJEB7759    real title: A Catalogue of the Mouse Gut Metagenome
    PRJNA278393  real title: Metagenome sequencing of the Hadza hunter-gatherer gut microbiota
    PRJEB1220    real title: A method for identifying metagenomic species ... co-abundance
    PRJNA313232  real title: Bacteroidales group S24-7 population genome recovery  <- the
                 paper's OWN deposit; the only correct link of the five
All five carry "human; mouse; koala; guinea pig" from the S24-7 paper.

Why the date rule alone does not catch it
-----------------------------------------
The specified rule is fanout >= 3 AND paper_deposit_gap >= 3. The S24-7 deposits are
from 2014-2016 against a 2016 paper, so the gaps are 2,1,1,1,0 -- all below 3. The date
discriminator only detects OLD reuse; this is RECENT reuse and slips through. The five
are therefore added by an explicit override, and paper_link_suspect_reason records which
mechanism flagged each study so the two are never confused.

No LLM re-run: everything here is recomputed from existing columns plus ENA first_public.

New columns: pmid_fanout, paper_deposit_gap, paper_link_suspect,
             paper_link_suspect_reason, paper_host_reliable, is_multi_host_paper
Nothing is deleted or re-linked, and host_source_summary is only touched for the
reuse-suspect relabel in step 4.
"""

import io
import shutil
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

R = "../results/"
CATALOG = R + "yes_catalog.tsv"
BACKUP = R + "yes_catalog.tsv.bak5"
ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"

FANOUT_MIN = 3
GAP_MIN = 3
BATCH = 100

# confirmed reuse the date rule cannot see (recent reuse, gap < 3)
KNOWN_REUSE = {
    "PRJEB1220": 27388460, "PRJEB6456": 27388460, "PRJEB7759": 27388460,
    "PRJNA278393": 27388460, "PRJNA313232": 27388460,
}
# PRJNA313232 is the S24-7 paper's OWN deposit -- correctly linked, so not a suspect
KNOWN_REUSE_EXCEPT = {"PRJNA313232"}

NEW_COLS = ["pmid_fanout", "paper_deposit_gap", "paper_link_suspect",
            "paper_link_suspect_reason", "paper_host_reliable", "is_multi_host_paper"]


def session():
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=5, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"])))
    return s


SESSION = session()


def fetch_first_public(accs):
    out = {}
    for i in range(0, len(accs), BATCH):
        batch = accs[i:i + BATCH]
        q = " OR ".join(f'study_accession="{a}"' for a in batch)
        r = SESSION.post(ENA_SEARCH, data={
            "result": "study", "query": q,
            "fields": "study_accession,first_public", "format": "tsv", "limit": 0},
            timeout=180)
        if r.status_code == 200 and r.text.strip():
            df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
            for _, row in df.iterrows():
                if row["study_accession"] in batch:   # umbrellas also return children
                    out[row["study_accession"]] = row.get("first_public")
    return out


def split_hosts(v):
    """Distinct hosts in a ';'-joined value, collapsing general-vs-specific.

    'human; human patients' is ONE host at two specificities, per the domain rule that
    a general term and a specific one are agreement, never two organisms.
    """
    if not isinstance(v, str) or not v.strip():
        return []
    parts = [p.strip().lower() for p in v.split(";") if p.strip()]
    parts = [p for p in parts if p not in ("", "n/a", "none", "nan")]
    kept = []
    for p in sorted(set(parts), key=len):        # shortest first = most general
        if not any(p in k or k in p for k in kept):
            kept.append(p)
    return kept


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    shutil.copy(CATALOG, BACKUP)
    print(f"backup -> {BACKUP}\n")
    cat = cat.drop(columns=[c for c in NEW_COLS if c in cat.columns])

    cat["pmid_i"] = pd.to_numeric(cat.pmid, errors="coerce")

    # ---------------- step 1: reuse linkage ----------------
    fan = cat.pmid_i.value_counts()
    cat["pmid_fanout"] = cat.pmid_i.map(fan).fillna(0).astype(int)

    need = sorted(cat.loc[cat.pmid_fanout >= 2, "study_accession"].astype(str).unique())
    print(f"STEP 1  fetching first_public for {len(need)} accessions (fanout >= 2)...")
    dates = fetch_first_public(need)
    print(f"        dated {len(dates)}/{len(need)}")

    dep_year = cat.study_accession.map(
        {a: (int(d[:4]) if isinstance(d, str) and d[:4].isdigit() else None)
         for a, d in dates.items()})
    pub_year = pd.to_numeric(cat.pub_year, errors="coerce")
    cat["paper_deposit_gap"] = (pub_year - dep_year)

    by_rule = (cat.pmid_fanout >= FANOUT_MIN) & (cat.paper_deposit_gap >= GAP_MIN)
    by_override = (cat.study_accession.isin(KNOWN_REUSE) &
                   ~cat.study_accession.isin(KNOWN_REUSE_EXCEPT))

    cat["paper_link_suspect"] = by_rule | by_override
    cat["paper_link_suspect_reason"] = [
        "rule:fanout+gap" if r else ("override:confirmed-reuse" if o else "")
        for r, o in zip(by_rule, by_override)]

    print(f"        flagged by rule     : {int(by_rule.sum())}")
    print(f"        flagged by override : {int(by_override.sum())}")
    print(f"        paper_link_suspect  : {int(cat.paper_link_suspect.sum())}")
    s247 = cat[cat.study_accession.isin(KNOWN_REUSE)]
    print("\n        S24-7 fan-out (PMID 27388460) -- rule vs override:")
    for _, r in s247.iterrows():
        print(f"          {r.study_accession:14s} fanout={r.pmid_fanout} "
              f"gap={r.paper_deposit_gap}  rule={bool(by_rule[r.name])}  "
              f"suspect={r.paper_link_suspect}  {r.paper_link_suspect_reason}")

    # ---------------- step 2: reliability ----------------
    cat["paper_host_reliable"] = ~cat.paper_link_suspect
    print(f"\nSTEP 2  paper_host_reliable = True for "
          f"{int(cat.paper_host_reliable.sum())} studies")

    # ---------------- step 3: multi-host ----------------
    # Evaluate each field SEPARATELY. Unioning them manufactures false multi-host:
    # PRJEB55711 title "patients with systemic lupus erythematosus" + abstract "human"
    # is ONE host phrased two ways, not two hosts.
    n_abs = cat.host_abstract_value.apply(lambda v: len(split_hosts(v)))
    n_tit = cat.host_title_value.apply(lambda v: len(split_hosts(v)))
    cat["_n_paper_hosts"] = pd.concat([n_abs, n_tit], axis=1).max(axis=1)

    # repair the model-organism flag where the extracted host text plainly says so
    says_model = (cat.host_abstract_value.fillna("").str.contains(r"model of|\(model", case=False, regex=True) |
                  cat.host_title_value.fillna("").str.contains(r"model of|\(model", case=False, regex=True))
    misfired = says_model & ~cat.is_model_organism.fillna(False).astype(bool)
    if misfired.any():
        print(f"\n        REPAIRED is_model_organism on {int(misfired.sum())} studies "
              "whose host text says 'model of':")
        for _, r in cat[misfired].iterrows():
            print(f"          {r.study_accession:14s} {str(r.host_abstract_value)[:58]}")
        cat.loc[misfired, "is_model_organism"] = True

    model = cat.is_model_organism.fillna(False).astype(bool)
    multi_raw = cat._n_paper_hosts >= 2
    cat["is_multi_host_paper"] = multi_raw & ~model & cat.paper_host_reliable

    excluded_reuse = cat[multi_raw & ~model & ~cat.paper_host_reliable]
    print(f"\nSTEP 3  >=2 paper hosts                    : {int(multi_raw.sum())}")
    print(f"        ...minus model-organism studies    : {int((multi_raw & ~model).sum())}")
    print(f"        ...minus reuse-suspects            : {int(cat.is_multi_host_paper.sum())}"
          "   <- is_multi_host_paper")
    print(f"\n        EXCLUDED as reuse-suspect (not silent): {len(excluded_reuse)}")
    for _, r in excluded_reuse.iterrows():
        print(f"          {r.study_accession:14s} {str(r.host_abstract_value)[:58]}")

    for acc in ("PRJNA940065", "PRJNA1148165"):
        r = cat[cat.study_accession == acc]
        if len(r):
            r = r.iloc[0]
            print(f"\n        CHECK {acc}: is_model_organism={r.is_model_organism} "
                  f"disease={r.model_organism_disease!r}")
            print(f"          abstract host: {str(r.host_abstract_value)[:80]}")
            print(f"          -> is_multi_host_paper={r.is_multi_host_paper}")

    # ---------------- step 4: human-by-biome ----------------
    hb = cat.human_by_biome_label.fillna(False).astype(bool)
    impure = hb & (cat._n_paper_hosts >= 2)
    print(f"\nSTEP 4  human_by_biome_label before        : {int(hb.sum())}")
    print(f"        ...naming >=2 hosts (not human alone): {int(impure.sum())}")
    for _, r in cat[impure].iterrows():
        print(f"          {r.study_accession:14s} reliable={r.paper_host_reliable}  "
              f"{str(r.host_abstract_value)[:52]}")

    cat.loc[impure, "human_by_biome_label"] = False
    # reliable + multi-host -> genuinely multi-host
    promote = impure & cat.paper_host_reliable & ~model
    cat.loc[promote, "is_multi_host_paper"] = True
    # unreliable -> the recovered host cannot be trusted at all
    relabel = impure & ~cat.paper_host_reliable
    cat.loc[relabel, "host_source_summary"] = "paper-link-suspect"

    print(f"\n        promoted to is_multi_host_paper    : {int(promote.sum())}")
    print(f"        relabelled 'paper-link-suspect'    : {int(relabel.sum())}")

    hb2 = cat.human_by_biome_label.fillna(False).astype(bool)
    samples = pd.to_numeric(cat.n_total_samples, errors="coerce").fillna(0)
    print(f"\n        human_by_biome_label after         : {int(hb2.sum())}")
    print(f"        samples over the CLEAN set         : {int(samples[hb2].sum()):,}"
          f"   (was {int(samples[hb].sum()):,} over the old {int(hb.sum())})")

    # ---------------- step 5: corrected picture ----------------
    print(f"\n{'='*66}\nSTEP 5  CORRECTED PICTURE\n{'='*66}")
    sus = cat[cat.paper_link_suspect]
    print(f"paper_link_suspect: {len(sus)} studies")
    print("  driven by these PMIDs:")
    for p, n in sus.pmid_i.value_counts().items():
        rows = sus[sus.pmid_i == p]
        print(f"    PMID {int(p)}  {n} accessions  fanout={rows.pmid_fanout.iloc[0]}  "
              f"reason={rows.paper_link_suspect_reason.iloc[0]}")

    mh = cat.is_multi_host_paper.fillna(False).astype(bool)
    div = cat.sample_host_pattern.isin(["DIVERSE", "DIVERSE_AND_SPARSE"])
    print(f"\nmulti-host by PAPER   (is_multi_host_paper)      : {int(mh.sum())}")
    print(f"multi-host by SAMPLES (DIVERSE*)                 : {int(div.sum())}")
    print(f"  overlap (both methods agree)                   : {int((mh & div).sum())}")
    print(f"  UNION = true multi-host population             : {int((mh | div).sum())}")
    print(f"\nhuman_by_biome_label (clean): {int(hb2.sum())} studies, "
          f"{int(samples[hb2].sum()):,} samples")

    cat = cat.drop(columns=["pmid_i", "_n_paper_hosts"])
    bad = [c for c in cat.columns if c.endswith("_x") or c.endswith("_y")]
    assert not bad, f"suffixed columns: {bad}"
    cat.to_csv(CATALOG, sep="\t", index=False)
    print(f"\ncolumns: {len(cat.columns)}   Saved -> {CATALOG}")
    return cat


if __name__ == "__main__":
    run()
