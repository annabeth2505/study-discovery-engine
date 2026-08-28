"""
detect_discrepancies.py
Cross-check the PAPER (PubMed title + abstract) against the DEPOSIT (ENA host,
samples, body sites, library strategies) to find REPORTING ERRORS.
Runs on paper-linked studies only (both sources present).

Checks (with false-positive guards baked in after hand-validation):
  1. host_mismatch   -- paper's SAMPLED organism genuinely differs from deposit host
  2. host_missing    -- paper names a host but deposit field is blank/generic
  3. gut_sample_issue-- deposit sample-level fields reveal non-gut samples
  4. multihost_issue -- deposit sample data shows MULTIPLE distinct host species
                        the paper presents as one (must be explicit in the data)

Adjudication: body-site/sample-type -> deposit sample metadata is ground truth;
host identity -> paper is ground truth.

Output:
  ../results/discrepancies.tsv          -- ALL checked studies (for catalog merge)
  ../results/discrepancies_flagged.tsv  -- flagged subset 
Checkpointed + resumable.
"""

import os
import json
import time
import pandas as pd
import anthropic

MODEL = "claude-haiku-4-5-20251001"
CATALOG = "../results/yes_catalog.tsv"
CHECKPOINT = "../results/discrepancy_verdicts.json"
OUT_ALL = "../results/discrepancies.tsv"
OUT_FLAGGED = "../results/discrepancies_flagged.tsv"

LIMIT = None   # 30 for validation first

client = anthropic.Anthropic()

SYSTEM = (
    "You are a metagenomics data-quality auditor. For each study you get TWO INDEPENDENT "
    "sources: (A) the PAPER (PubMed title + abstract) and (B) the DEPOSIT (ENA host field, "
    "sample body-sites, isolation sources, sample titles, library-strategy mix). Find "
    "genuine REPORTING ERRORS by comparing them. Adjudication: for BODY-SITE / SAMPLE-TYPE "
    "the DEPOSIT sample-level metadata is ground truth; for HOST IDENTITY the PAPER is "
    "ground truth. "
    "CRITICAL RULES to avoid false positives — apply them strictly: "
    "(1) MODEL ORGANISMS: gut-microbiome papers about a HUMAN disease are very often done "
    "in an ANIMAL MODEL. If the abstract mentions or implies an animal model (mouse model, "
    "mice, murine, rats, germ-free/gnotobiotic animals, 'animal model of'), then a "
    "mouse/rat host in the deposit is CORRECT — do NOT flag host_mismatch. Only flag "
    "host_mismatch when the paper describes samples taken DIRECTLY from an organism the "
    "deposit clearly contradicts. "
    "(2) GENERAL vs SPECIFIC is AGREEMENT, never a conflict: a paper umbrella term "
    "('ruminant','mammal','bird','salmonid','rodent','fish') matching a specific deposit "
    "species (buffalo, mouse, tit, salmon) is AGREEMENT. Common name vs Latin name (mouse "
    "vs Mus musculus, buffalo vs Bubalus bubalis) is AGREEMENT. NEVER infer that a general "
    "term implies multiple species. "
    "(3) EMPTY is not WRONG: if the deposit host field is blank or generic ('metagenome', "
    "'gut metagenome', none) while the paper names a host, that is host_missing (a metadata "
    "gap), NOT host_mismatch. "
    "(4) MULTI-HOST requires EXPLICIT EVIDENCE: only flag multihost_issue when the DEPOSIT's "
    "own sample-level data (host_species_list, biome labels, sample titles) explicitly names "
    "MULTIPLE DISTINCT host species while the paper presents a single one. NEVER flag it "
    "because a paper uses a general term, and NEVER flag it merely because the library "
    "strategy is mixed (WGS + 16S + RNA-Seq together is normal). "
    "Return ONLY valid JSON."
)

PROMPT = """Study {acc}.

=== SOURCE A: THE PAPER (PubMed) ===
Title: {title}
Abstract: {abstract}

=== SOURCE B: THE DEPOSIT (ENA) ===
Host field: {host}
Host species seen across samples ({n_host} distinct): {host_list}
Sample body sites: {body_sites}
Isolation sources: {isolation}
Sample titles: {sample_titles}
Biome labels: {biomes}
Library strategy mix: {strat_mix}

Check FOUR things. Apply the false-positive rules from the system prompt strictly —
when in doubt, do NOT flag.

1. host_mismatch: the organism the paper SAMPLES genuinely differs from the deposit host.
   Exclude animal-model studies of human disease. Exclude general-vs-specific and
   common-vs-Latin. Flag true only for a real, direct conflict.
   IMPORTANT: if the mouse/rat host is APPROPRIATE (animal model, or any reason the deposit
   host is actually fine), set flagged=FALSE. Do NOT set flagged=true and then say in the
   note that it is appropriate — the flag must match your judgment. If your note would say
   "appropriate / not an error / expected / correct", then flagged MUST be false.
2. host_missing: the paper names a host but the deposit host field is blank or generic.
3. gut_sample_issue: the deposit sample-level fields reveal NON-gut samples (surfaces,
   soil, oral, skin, water, blood) in a study treated as gut. Deposit metadata is ground
   truth here.
4. multihost_issue: the deposit's sample-level data explicitly names MULTIPLE DISTINCT
   host species that the paper presents as a single host. Must be explicit in the data,
   not inferred from a general term. Mixed library strategy alone is NOT this.

Return ONLY:
{{"host_mismatch":{{"flagged":false,"paper_value":"","deposit_value":"","note":""}},
"host_missing":{{"flagged":false,"paper_host":"","note":""}},
"gut_sample_issue":{{"flagged":false,"evidence":"","note":""}},
"multihost_issue":{{"flagged":false,"evidence":"","note":""}},
"any_flag":false,
"summary":"one sentence overall"}}"""


def call_llm(prompt, retries=6):
    for a in range(retries):
        try:
            m = client.messages.create(model=MODEL, max_tokens=700, system=SYSTEM,
                                       messages=[{"role": "user", "content": prompt}])
            t = m.content[0].text.strip().replace("```json", "").replace("```", "").strip()
            return json.loads(t)
        except anthropic.RateLimitError:
            time.sleep(13 * (a + 1))
        except json.JSONDecodeError:
            return {"any_flag": None, "summary": "PARSE_ERROR"}
        except Exception as e:
            return {"any_flag": None, "summary": f"ERROR: {str(e)[:120]}"}
    return {"any_flag": None, "summary": "ERROR: max retries"}


def trunc(x, n):
    s = str(x) if pd.notna(x) else "(none)"
    return s[:n] + "..." if len(s) > n else s


def run():
    cat = pd.read_csv(CATALOG, sep="\t")
    checkable = cat[cat["pmid"].notna() &
                    cat["abstract"].notna() &
                    (cat["abstract"].astype(str).str.strip() != "")].copy()
    if LIMIT:
        checkable = checkable.head(LIMIT)
        print(f"** VALIDATION RUN — first {LIMIT} **")
    print(f"Cross-checkable studies: {len(checkable)}")

    verdicts = {}
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            verdicts = json.load(f)
        print(f"Resuming — {len(verdicts)} done")

    todo = checkable[~checkable["study_accession"].isin(verdicts.keys())]
    print(f"Checking {len(todo)}...\n")

    for i, r in enumerate(todo.itertuples()):
        acc = r.study_accession
        if i % 10 == 0:
            print(f"  {i}/{len(todo)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(verdicts, f)
        prompt = PROMPT.format(
            acc=acc,
            title=trunc(getattr(r, "paper_title", ""), 300),
            abstract=trunc(getattr(r, "abstract", ""), 1200),
            host=trunc(getattr(r, "host_species", ""), 200),
            n_host=getattr(r, "n_host_species", "?"),
            host_list=trunc(getattr(r, "host_species_list", ""), 400),
            body_sites=trunc(getattr(r, "body_sites", ""), 200),
            isolation=trunc(getattr(r, "isolation_sources", ""), 200),
            sample_titles=trunc(getattr(r, "sample_titles_cache", ""), 250),
            biomes=trunc(getattr(r, "biome_labels", ""), 200),
            strat_mix=getattr(r, "library_strategy_mix", "(none)"),
        )
        verdicts[acc] = call_llm(prompt)
        time.sleep(0.2)

    with open(CHECKPOINT, "w") as f:
        json.dump(verdicts, f)

    rows = []
    for acc, v in verdicts.items():
        def g(k, f): return (v.get(k) or {}).get(f, "") if isinstance(v.get(k), dict) else ""
        hm = g("host_mismatch", "flagged") is True
        # guard: if the model's OWN note says the host is fine, trust the reasoning
        # over the mechanical boolean (catches self-contradicting flags)
        note_hm = str(g("host_mismatch", "note")).lower()
        summary = str(v.get("summary", "")).lower()
        exonerating = ["appropriate", "not a discrepancy", "no error", "no mismatch",
                       "is correct", "animal model", "mouse model", "expected", "no genuine"]
        if hm and any(p in note_hm or p in summary for p in exonerating):
            hm = False
        hmiss = g("host_missing", "flagged") is True
        gs = g("gut_sample_issue", "flagged") is True
        mh = g("multihost_issue", "flagged") is True
        flags = [name for name, on in
                 [("host_mismatch", hm), ("host_missing", hmiss),
                  ("gut_sample_issue", gs), ("multihost_issue", mh)] if on]
        rows.append({
            "study_accession": acc,
            "discrepancy_flags": "; ".join(flags),
            "discrepancy_host": (f"paper: {g('host_mismatch','paper_value')} vs "
                                 f"deposit: {g('host_mismatch','deposit_value')}") if hm else "",
            "discrepancy_host_missing": (f"paper names: {g('host_missing','paper_host')}; "
                                         f"deposit blank/generic") if hmiss else "",
            "discrepancy_gut": g("gut_sample_issue", "evidence") if gs else "",
            "discrepancy_multihost": g("multihost_issue", "evidence") if mh else "",
            "discrepancy_notes": v.get("summary", ""),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_ALL, sep="\t", index=False)

    flagged = out[out["discrepancy_flags"].astype(str).str.len() > 0].copy()
    flagged = flagged.merge(cat[["study_accession", "paper_title", "host_species"]],
                            on="study_accession", how="left")
    flagged.to_csv(OUT_FLAGGED, sep="\t", index=False)

    parse_err = sum(1 for v in verdicts.values() if v.get("summary", "").startswith(("PARSE_ERROR", "ERROR")))
    print(f"\n{'='*54}\nDISCREPANCY DETECTION COMPLETE\n{'='*54}")
    print(f"Studies checked: {len(out)}")
    print(f"Flagged with >=1 discrepancy: {len(flagged)}")
    print(f"\nBy type:")
    print(f"  host_mismatch (real conflict): {out['discrepancy_flags'].str.contains('host_mismatch').sum()}")
    print(f"  host_missing (blank/generic):  {out['discrepancy_flags'].str.contains('host_missing').sum()}")
    print(f"  gut_sample_issue:              {out['discrepancy_flags'].str.contains('gut_sample_issue').sum()}")
    print(f"  multihost_issue:               {out['discrepancy_flags'].str.contains('multihost_issue').sum()}")
    print(f"\nParse errors: {parse_err}")
    print(f"Saved -> {OUT_ALL}  (all checked, for catalog merge)")
    print(f"Saved -> {OUT_FLAGGED}  (flagged only)")
    return out, flagged


if __name__ == "__main__":
    run()