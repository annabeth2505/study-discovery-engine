"""
parse_sample_patterns.py
Characterize each study by the PATTERN of its sample-level host annotations --
the "uniform vs diverse" question.

Input : ../results/samples.tsv            (one row per biosample)
Output: ../results/study_sample_patterns.tsv  (one row per study)

Never writes into yes_catalog.tsv -- that file is heavily layered and merge-collision
prone. This is a standalone file to be merged later, deliberately.

TWO MEASURES, KEPT STRICTLY SEPARATE:
  1. COMPLETENESS -- what fraction of a study's samples have a real resolved host.
     A blank is host_resolved null / UNRESOLVED / a sentinel string.
  2. DIVERSITY -- among ONLY the populated samples, how many DISTINCT hosts.
     A blank is "unknown", never a distinct host, and never counts toward diversity.

Diversity is counted on host_tax_id, falling back to host_resolved only where tax_id
is missing. Free text has variants (RAT / mice / C57BL/6J / Mus musculus are all one
species) that would otherwise inflate the count into false diversity; tax_id collapses
them. Within a study, a name seen elsewhere WITH a tax_id is mapped onto that tax_id,
so "Mus musculus" (no tax_id) and 10090 don't split into two phantom hosts.

  2x2 label (sample_host_pattern):
                    1 distinct host      >=2 distinct hosts
    complete >0.8   UNIFORM              DIVERSE
    complete <=0.8  UNIFORM_WITH_GAPS    DIVERSE_AND_SPARSE
  Plus NO_HOST_DATA when a study has ZERO populated hosts -- the 2x2 has no cell for
  it (0 distinct hosts is neither 1 nor >=2), and calling it UNIFORM_WITH_GAPS would
  assert a uniformity the data never showed. A gap stays visible as a gap.

STAGE 2 -- LLM adjudication, DIVERSE / DIVERSE_AND_SPARSE studies only. The pandas
pass calls anything with >=2 distinct populated hosts "diverse", but some of that is
general-vs-specific labelling of ONE animal (donkey + "herbivore" = not multi-host).
Claude judges which. UNIFORM / UNIFORM_WITH_GAPS / NO_HOST_DATA are never sent.
Checkpointed + resumable.

Usage:
  python parse_sample_patterns.py --stage pandas --limit 20      # test the pass
  python parse_sample_patterns.py --stage pandas --studies A,B   # specific studies
  python parse_sample_patterns.py --stage pandas                 # full pandas pass
  python parse_sample_patterns.py --stage llm                    # adjudicate diverse
"""

import os
import re
import json
import time
import argparse
import pandas as pd

SAMPLES = "../results/samples.tsv"
OUT = "../results/study_sample_patterns.tsv"
CHECKPOINT = "../results/host_pattern_verdicts.json"

MODEL = "claude-haiku-4-5-20251001"
COMPLETE_THRESHOLD = 0.8

SENTINELS = {"missing", "na", "n/a", "not applicable", "none", "null", "unknown",
             "not collected", "not provided", "uncalculated", "unresolved", ""}

COLS = ["study_accession", "n_samples", "n_host_populated", "pct_host_complete",
        "n_distinct_hosts", "distinct_host_list", "sample_host_pattern",
        "seq_type_mix", "n_distinct_seq_types", "has_genomic", "n_genomic",
        "host_pattern_confirmed", "host_pattern_llm_note"]


def _clean(v):
    """Usable string, or None for blank / sentinel / UNRESOLVED."""
    if v is None or (isinstance(v, float) and pd.isna(v)) or pd.isna(v):
        return None
    s = str(v).strip()
    return None if not s or s.lower() in SENTINELS else s


def _norm(name):
    """Normalize a host name for grouping: case and whitespace only, nothing clever."""
    return re.sub(r"\s+", " ", str(name).strip().lower())


def study_pattern(sub):
    """All per-study measures for one study's sample rows."""
    n_samples = len(sub)

    # ---- 1. COMPLETENESS: host_resolved is the resolved field ----
    populated = sub[[_clean(v) is not None for v in sub["host_resolved"]]]
    n_pop = len(populated)
    pct = round(n_pop / n_samples, 3) if n_samples else 0.0

    # ---- 2. DIVERSITY: populated samples only, keyed on tax_id ----
    # map name -> tax_id from rows carrying both, so a name-only row joins its species
    name_to_tax = {}
    for nm, tx in zip(populated["host_resolved"], populated["host_tax_id"]):
        nm, tx = _clean(nm), _clean(tx)
        if nm and tx:
            name_to_tax.setdefault(_norm(nm), tx)

    keys = {}   # diversity key -> representative display name
    for nm, tx in zip(populated["host_resolved"], populated["host_tax_id"]):
        nm, tx = _clean(nm), _clean(tx)
        if nm is None:
            continue
        key = tx or name_to_tax.get(_norm(nm)) or _norm(nm)
        keys.setdefault(key, nm)

    n_distinct = len(keys)
    host_list = "; ".join(sorted(keys.values()))

    # ---- 2x2 label ----
    if n_pop == 0:
        pattern = "NO_HOST_DATA"
    elif n_distinct == 1:
        pattern = "UNIFORM" if pct > COMPLETE_THRESHOLD else "UNIFORM_WITH_GAPS"
    else:
        pattern = "DIVERSE" if pct > COMPLETE_THRESHOLD else "DIVERSE_AND_SPARSE"

    # ---- sequencing types ----
    src = sub["library_source"].dropna().value_counts()
    seq_mix = "; ".join(f"{k}:{v}" for k, v in src.items())
    n_genomic = int(src.get("GENOMIC", 0))

    return {
        "n_samples": n_samples,
        "n_host_populated": n_pop,
        "pct_host_complete": pct,
        "n_distinct_hosts": n_distinct,
        "distinct_host_list": host_list,
        "sample_host_pattern": pattern,
        "seq_type_mix": seq_mix,
        "n_distinct_seq_types": int(len(src)),
        "has_genomic": bool(n_genomic > 0),
        "n_genomic": n_genomic,
    }


def run_pandas(limit=None, studies=None, out=OUT):
    d = pd.read_csv(SAMPLES, sep="\t", dtype=str)

    if studies:
        d = d[d["study_accession"].isin(studies)]
    accs = sorted(d["study_accession"].dropna().unique())
    if limit:
        accs = accs[:limit]
        d = d[d["study_accession"].isin(accs)]

    rows = []
    for acc, sub in d.groupby("study_accession"):
        rows.append({"study_accession": acc, **study_pattern(sub)})

    prof = pd.DataFrame(rows)
    # stage-2 columns exist from the start; the LLM fills only the diverse subset
    prof["host_pattern_confirmed"] = prof["sample_host_pattern"]
    prof["host_pattern_llm_note"] = ""
    prof = prof[COLS].sort_values("study_accession").reset_index(drop=True)
    prof.to_csv(out, sep="\t", index=False)

    print(f"{'='*58}\nSAMPLE-PATTERN PASS (pandas)\n{'='*58}")
    print(f"Studies: {len(prof)}   Samples: {int(prof.n_samples.sum())}")
    print()
    print("pattern distribution:")
    for k, v in prof.sample_host_pattern.value_counts().items():
        n = int(prof.loc[prof.sample_host_pattern == k, "n_samples"].sum())
        print(f"  {k:20s} {v:5d} studies  ({n:7d} samples)")
    print()
    n_div = int(prof.sample_host_pattern.isin(["DIVERSE", "DIVERSE_AND_SPARSE"]).sum())
    print(f"-> {n_div} studies queued for LLM adjudication")
    print(f"multi-seq-type studies: {int((prof.n_distinct_seq_types > 1).sum())}")
    print(f"studies containing GENOMIC samples: {int(prof.has_genomic.sum())}")
    print(f"\nSaved -> {out}")
    return prof


# ----------------------------- stage 2: LLM ------------------------------------

SYSTEM = (
    "You are a metagenomics metadata auditor. A study's samples carry host labels written "
    "by the submitter. You are given the DISTINCT host labels found across one study's "
    "samples, and must judge whether they represent GENUINELY DIFFERENT ANIMALS or are "
    "just different ways of labelling ONE animal. "
    "CRITICAL RULES: "
    "(1) GENERAL vs SPECIFIC is the SAME animal, not diversity: 'herbivore' + 'Equus "
    "asinus' = one donkey study labelled inconsistently. Same for 'mammal', 'rodent', "
    "'ruminant', 'fish', 'bird' sitting alongside a specific species THEY CONTAIN. "
    "(2) Common name vs Latin name is the SAME animal: mouse / Mus musculus, "
    "chicken / Gallus gallus, buffalo / Bubalus bubalis. "
    "(3) Strain and breed names are the SAME species: C57BL/6J, BALB/c, ICR = Mus "
    "musculus. Subspecies (Gallus gallus subsp domesticus) = the species. "
    "(4) A general term that does NOT contain the specific one is genuine diversity: "
    "'bird' alongside 'Equus asinus' are different animals -- a bird is not a donkey. "
    "(5) Two distinct named species (Mus musculus + Rattus norvegicus) = GENUINE, even "
    "when both are model organisms. A human + animal-model mix is also GENUINE. "
    "Return ONLY valid JSON."
)

PROMPT = """Study {acc}.

Distinct host labels across its {n_pop} host-annotated samples ({n_distinct} distinct):
{host_list}

Are these GENUINELY different animals, or variants/general-vs-specific of ONE animal?

- "GENUINE"     -> at least two truly different animals (real multi-host study, or a
                   submitter error putting a wrong animal in the deposit)
- "SAME_ANIMAL" -> all labels describe ONE animal at different levels of precision,
                   or as name variants / strains

Return ONLY (list AT MOST 5 examples in distinct_animals, never the whole list --
a long echo truncates the response and the verdict is lost):
{{"verdict":"GENUINE","distinct_animals":["..."],"note":"one short sentence"}}"""


def call_llm(client, prompt, retries=6):
    import anthropic
    for a in range(retries):
        try:
            m = client.messages.create(model=MODEL, max_tokens=1000, system=SYSTEM,
                                       messages=[{"role": "user", "content": prompt}])
            t = m.content[0].text.strip().replace("```json", "").replace("```", "").strip()
            return json.loads(t)
        except anthropic.RateLimitError:
            time.sleep(13 * (a + 1))
        except json.JSONDecodeError:
            return {"verdict": "PARSE_ERROR", "note": "unparseable model output"}
        except Exception as e:
            return {"verdict": "ERROR", "note": str(e)[:120]}
    return {"verdict": "ERROR", "note": "max retries"}


# a downgraded label keeps the study's completeness half intact
DOWNGRADE = {"DIVERSE": "UNIFORM", "DIVERSE_AND_SPARSE": "UNIFORM_WITH_GAPS"}


def run_llm(out=OUT, checkpoint=CHECKPOINT, limit=None):
    import anthropic
    from pathlib import Path
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).resolve().parent.parent / ".env")   # same as host_enrichment.py
    client = anthropic.Anthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))

    prof = pd.read_csv(out, sep="\t")
    diverse = prof[prof.sample_host_pattern.isin(["DIVERSE", "DIVERSE_AND_SPARSE"])]
    if limit:
        diverse = diverse.head(limit)

    done = {}
    if os.path.exists(checkpoint):
        with open(checkpoint) as f:
            done = json.load(f)
        print(f"Resuming -- {len(done)} studies already adjudicated")

    todo = [r for _, r in diverse.iterrows() if r.study_accession not in done]
    print(f"Adjudicating {len(todo)} diverse studies with {MODEL}...")

    for i, r in enumerate(todo):
        if i % 10 == 0:
            print(f"  {i}/{len(todo)}")
            with open(checkpoint, "w") as f:
                json.dump(done, f)
        done[r.study_accession] = call_llm(client, PROMPT.format(
            acc=r.study_accession, n_pop=r.n_host_populated,
            n_distinct=r.n_distinct_hosts, host_list=r.distinct_host_list))
        time.sleep(0.3)

    with open(checkpoint, "w") as f:
        json.dump(done, f)

    # ---- fold verdicts back in; failures stay visible as NOT_CHECKED ----
    def confirmed(row):
        v = done.get(row.study_accession)
        if v is None:
            return row.sample_host_pattern          # never sent to the LLM
        verdict = v.get("verdict")
        if verdict == "GENUINE":
            return row.sample_host_pattern
        if verdict == "SAME_ANIMAL":
            return DOWNGRADE.get(row.sample_host_pattern, row.sample_host_pattern)
        return "NOT_CHECKED"                        # PARSE_ERROR / ERROR

    prof["host_pattern_confirmed"] = prof.apply(confirmed, axis=1)
    prof["host_pattern_llm_note"] = prof.study_accession.map(
        lambda a: (done.get(a) or {}).get("note", ""))
    prof.to_csv(out, sep="\t", index=False)

    adj = prof[prof.study_accession.isin(done)]
    print(f"\n{'='*58}\nLLM ADJUDICATION COMPLETE\n{'='*58}")
    print(f"Studies adjudicated: {len(adj)}")
    print(f"  confirmed GENUINE diversity : {int((adj.host_pattern_confirmed.isin(['DIVERSE','DIVERSE_AND_SPARSE'])).sum())}")
    print(f"  downgraded to uniform       : {int((adj.host_pattern_confirmed.isin(['UNIFORM','UNIFORM_WITH_GAPS'])).sum())}")
    print(f"  NOT_CHECKED (LLM failure)   : {int((adj.host_pattern_confirmed == 'NOT_CHECKED').sum())}")
    print("\nfinal pattern distribution:")
    print(prof.host_pattern_confirmed.value_counts().to_string())
    print(f"\nSaved -> {out}")
    return prof


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--stage", choices=["pandas", "llm", "all"], default="pandas")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--studies", type=str, default=None,
                   help="comma-separated study accessions")
    p.add_argument("--out", type=str, default=OUT)
    a = p.parse_args()
    studies = [s.strip() for s in a.studies.split(",")] if a.studies else None

    if a.stage in ("pandas", "all"):
        run_pandas(limit=a.limit, studies=studies, out=a.out)
    if a.stage in ("llm", "all"):
        run_llm(out=a.out, limit=a.limit)
