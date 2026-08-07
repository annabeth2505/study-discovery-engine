"""
audit_catalog_v3.py
Four-criteria LLM audit of the enriched catalog. Reads already-populated columns
(no ENA re-fetch).

KEY PRINCIPLE (corrected): the audit judges VALIDITY (is this a real animal gut WGS
study?), NOT SCOPE or COMPLETENESS. It never excludes a study for being human or for
having an unresolved host -- those are decisions for the curator, not the classifier.

EXCLUDE only for genuine disqualifiers:
  - gut = no (skin/oral/soil/etc., or not a microbiome study)
  - no WGS content (n_wgs_samples == 0 AND library_source != METAGENOMIC)
  - body_site = non_gut
  - host = incidental (trace host in a survey about something else)
REVIEW (flag, don't exclude): gut uncertain, low confidence, or method mismatch.
KEEP everything else -- INCLUDING human, unclear-host, genus-only, and genuine
multi-host comparative studies. Human is flagged (is_human) for the separate scope call.

Adds llm_notes column: the audit's full reasoning inline for every study.
Output: ../results/catalog_audit.tsv
"""

import os
import json
import time
import pandas as pd
import anthropic

MODEL = "claude-haiku-4-5-20251001"
CATALOG = "../results/yes_catalog.tsv"
CHECKPOINT = "../results/audit_verdicts.json"
OUT_PATH = "../results/catalog_audit.tsv"

LIMIT = None   # 30 for validation; None = full catalog

client = anthropic.Anthropic()

SYSTEM = (
    "You are a rigorous metagenomics data curator auditing a catalog of animal gut "
    "shotgun (WGS) metagenomics studies. Assess four criteria per study, based ONLY on the "
    "evidence provided. Distinguish STATED facts from INFERENCE; when inferring, say so and "
    "lower confidence. If evidence is insufficient, return 'uncertain' rather than guessing. "
    "Watch for: (a) multi-host SURVEYS where the tagged animal is only an incidental/trace "
    "sample in a study about something else -> mark host 'incidental'; but a genuine "
    "comparative study OF several animals' guts is NOT incidental -> mark 'multi_host'. "
    "(b) a stated host contradicting the paper title -> flag mismatch. (c) a study that only "
    "MENTIONS the gut but sampled skin/oral -> it is a skin/oral study. "
    "CRITICAL: humans ARE animals. A human gut microbiome study is a valid animal gut study "
    "-> true_animal_gut = 'yes'. NEVER return 'no' just because the host is human. The human "
    "vs non-human distinction is captured ONLY in the is_human flag, never in the gut verdict. "
    "Return ONLY valid JSON."
    
)

PROMPT = """AUDIT study {acc}.

EVIDENCE
- Paper title: {title}
- Abstract: {abstract}
- Host field: {host}
- Host species in samples ({n_host} distinct): {host_list}
- library_strategy_mix (per-sample): {strat_mix}
- WGS/shotgun sample count: {n_wgs}
- library_source: {lib_source}   (METAGENOMIC=community DNA; GENOMIC=single organism)
- Biome label(s): {biomes}
- Body site(s): {body_sites}
- Isolation source(s): {isolation}
- Sample titles: {sample_titles}

Assess FOUR criteria. Each: verdict, confidence (high/medium/low), evidence quote.

1. true_animal_gut: genuine ANIMAL GUT microbiome study? (yes / no / uncertain)
   NOTE: humans count as animals here -- a human gut study is 'yes'. Only mark 'no' for
   non-gut (skin/oral/soil), non-microbiome (host genome), or not-a-real-study cases.
2. body_site: are samples GUT vs non-gut?  (gut_explicit / gut_inferred / non_gut / ambiguous)
3. host_id: identify the host.  (single_clear / multi_host / incidental / genus_only / unclear)
   set "host" to the host, and "is_human" true/false.
4. method_note: does the described method agree with the ENA WGS fields?  (agree / mismatch / unclear)

Return ONLY:
{{"true_animal_gut":{{"verdict":"...","confidence":"...","evidence":"..."}},
"body_site":{{"verdict":"...","confidence":"...","evidence":"..."}},
"host_id":{{"verdict":"...","confidence":"...","host":"...","is_human":false,"evidence":"..."}},
"method_note":{{"verdict":"...","evidence":"..."}},
"rationale":"one to two sentences explaining the overall picture"}}"""


def call_llm(prompt, retries=6):
    for a in range(retries):
        try:
            m = client.messages.create(model=MODEL, max_tokens=600, system=SYSTEM,
                                       messages=[{"role": "user", "content": prompt}])
            t = m.content[0].text.strip().replace("```json", "").replace("```", "").strip()
            return json.loads(t)
        except anthropic.RateLimitError:
            time.sleep(13 * (a + 1))
        except json.JSONDecodeError:
            return {"rationale": "PARSE_ERROR"}
        except Exception as e:
            return {"rationale": f"ERROR: {str(e)[:120]}"}
    return {"rationale": "ERROR: max retries"}


def trunc(x, n):
    s = str(x) if pd.notna(x) else "(none)"
    return s[:n] + "..." if len(s) > n else s


def run():
    cat = pd.read_csv(CATALOG, sep="\t")
    if LIMIT:
        cat = cat.head(LIMIT)
        print(f"** VALIDATION RUN — first {LIMIT} **")

    verdicts = {}
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            verdicts = json.load(f)
        print(f"Resuming — {len(verdicts)} done")

    todo = cat[~cat["study_accession"].isin(verdicts.keys())]
    print(f"Auditing {len(todo)} studies...\n")

    for i, r in enumerate(todo.itertuples()):
        acc = r.study_accession
        if i % 10 == 0:
            print(f"  {i}/{len(todo)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(verdicts, f)
        prompt = PROMPT.format(
            acc=acc,
            title=trunc(getattr(r, "paper_title", ""), 300),
            abstract=trunc(getattr(r, "abstract", ""), 800),
            host=trunc(getattr(r, "host_species", ""), 200),
            n_host=getattr(r, "n_host_species", "?"),
            host_list=trunc(getattr(r, "host_species_list", ""), 300),
            strat_mix=getattr(r, "library_strategy_mix", "(none)"),
            n_wgs=getattr(r, "n_wgs_samples", "?"),
            lib_source=getattr(r, "library_source", "(none)"),
            biomes=trunc(getattr(r, "biome_labels", ""), 200),
            body_sites=trunc(getattr(r, "body_sites", ""), 150),
            isolation=trunc(getattr(r, "isolation_sources", ""), 150),
            sample_titles=trunc(getattr(r, "sample_titles_cache", ""), 200),
        )
        verdicts[acc] = call_llm(prompt)
        time.sleep(0.2)

    with open(CHECKPOINT, "w") as f:
        json.dump(verdicts, f)

    wgs_by_acc = dict(zip(cat["study_accession"], cat["n_wgs_samples"]))
    src_by_acc = dict(zip(cat["study_accession"], cat["library_source"]))

    rows = []
    for acc, v in verdicts.items():
        def g(crit, f): return (v.get(crit) or {}).get(f, "") if isinstance(v.get(crit), dict) else ""
        gut  = g("true_animal_gut", "verdict")
        gconf = g("true_animal_gut", "confidence")
        site = g("body_site", "verdict")
        host = g("host_id", "verdict")
        method = g("method_note", "verdict")
        n_wgs = wgs_by_acc.get(acc, 0) or 0
        src = src_by_acc.get(acc, "")
        is_human = str(g("host_id", "is_human")).lower() == "true"

        # ---- keep-rule: validity only, never scope/completeness ----
        has_wgs = (n_wgs > 0) or (src == "METAGENOMIC")
        parse_failed = (not gut) or str(gut).strip() in ("", "nan") or \
                       str(v.get("rationale","")).startswith(("PARSE_ERROR", "ERROR"))
        if parse_failed:
            keep = "REVIEW"          # couldn't audit -> human review, never silent exclude
        elif gut == "no" or site == "non_gut" or host == "incidental" or not has_wgs:
            keep = "EXCLUDE"
        elif gut == "uncertain" or gconf == "low" or method == "mismatch":
            keep = "REVIEW"
        else:
            keep = "KEEP"

        # ---- llm_notes: full reasoning inline ----
        notes = (
            f"GUT: {gut} ({gconf}) — {g('true_animal_gut','evidence')} || "
            f"BODY_SITE: {site} — {g('body_site','evidence')} || "
            f"HOST: {g('host_id','host')} [{host}]"
            f"{' (HUMAN)' if is_human else ''} — {g('host_id','evidence')} || "
            f"METHOD: {method} — {g('method_note','evidence')} || "
            f"SUMMARY: {v.get('rationale','')}"
        )

        rows.append({
            "study_accession": acc, "keep": keep,
            "gut": gut, "gut_conf": gconf,
            "body_site": site, "host_verdict": host, "host": g("host_id", "host"),
            "is_human": is_human, "method_note": method,
            "n_wgs_samples": n_wgs, "library_source": src,
            "llm_notes": notes,
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_PATH, sep="\t", index=False)

    print(f"\n{'='*52}\nAUDIT COMPLETE\n{'='*52}")
    print("Keep decision:");         print(out["keep"].value_counts().to_string())
    print("\nCriterion 1 gut:");      print(out["gut"].value_counts().to_string())
    print("\nCriterion 2 body_site:"); print(out["body_site"].value_counts().to_string())
    print("\nCriterion 3 host:");     print(out["host_verdict"].value_counts().to_string())
    print("\nmethod agreement:");     print(out["method_note"].value_counts().to_string())
    print(f"\nFlagged human (kept, for scope decision): {out['is_human'].sum()}")
    print(f"EXCLUDE: {(out['keep']=='EXCLUDE').sum()}  REVIEW: {(out['keep']=='REVIEW').sum()}  KEEP: {(out['keep']=='KEEP').sum()}")
    print(f"Saved -> {OUT_PATH}")
    return out


if __name__ == "__main__":
    run()