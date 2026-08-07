"""
audit_catalog_v2.py
Four-criteria LLM audit of the ENRICHED catalog. Reads the columns already
populated by enrich_catalog_from_ena.py + backfill passes -- does NOT re-fetch ENA.

Criteria (Claude scores each; keep-decision computed in CODE, not the prompt):
  1. true_animal_gut : genuine animal gut microbiome study?
  2. has_wgs         : does it contain WGS/shotgun samples? (from n_wgs_samples;
                       16S samples are fine as long as WGS samples exist)
  3. body_site       : gut vs non-gut (oral/skin/soil/water)?
  4. host_id         : identifiable single host? (flags incidental/multi-host & human)

Keep-rule (transparent, code-side):
  KEEP     if gut=yes AND n_wgs_samples>0 AND body_site in (gut_explicit,gut_inferred)
              AND host not incidental
  EXCLUDE  if gut=no OR n_wgs_samples==0 (unless library_source=METAGENOMIC salvage)
              OR body_site=non_gut OR host incidental
  REVIEW   otherwise (uncertain / low-confidence / mismatch)

Output: ../results/catalog_audit.tsv  (per-study verdicts + evidence + rationale)
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

LIMIT = None   # set to 30 for a validation run first; None = full catalog

client = anthropic.Anthropic()

SYSTEM = (
    "You are a rigorous metagenomics data curator auditing a catalog of animal gut "
    "shotgun (WGS) metagenomics studies. Assess four criteria per study and base every "
    "judgment ONLY on the evidence provided. Distinguish what is STATED from what you "
    "INFER; when inferring, say so and lower confidence. If evidence is insufficient, "
    "return 'uncertain' rather than guessing. Watch for: (a) multi-host surveys where the "
    "tagged animal is only an incidental/trace sample -> not a genuine study of it; "
    "(b) a stated host that contradicts the paper title's subject -> flag mismatch; "
    "(c) a study that merely MENTIONS the gut but sampled skin/oral -> it is a skin/oral "
    "study. Note explicitly when the host is human. Return ONLY valid JSON."
)

PROMPT = """AUDIT study {acc}.

EVIDENCE
- Paper title: {title}
- Abstract: {abstract}
- Host field: {host}
- Host species seen in samples ({n_host} distinct): {host_list}
- ENA library_strategy_mix (per-sample): {strat_mix}
- WGS/shotgun sample count: {n_wgs}
- ENA library_source: {lib_source}   (METAGENOMIC=community DNA; GENOMIC=single organism)
- Biome label(s): {biomes}
- Body site(s): {body_sites}
- Isolation source(s): {isolation}
- Sample titles: {sample_titles}

Assess FOUR criteria. Each: verdict, confidence (high/medium/low), evidence quote.

1. true_animal_gut: genuine study of an ANIMAL GUT microbiome?  (yes / no / uncertain)
2. body_site: are the samples GUT vs non-gut?  (gut_explicit / gut_inferred / non_gut / ambiguous)
3. host_id: can we identify the host?  (single_clear / multi_host / incidental / genus_only / unclear)
   -- 'incidental' = tagged animal is a trace host in a study about something else.
   -- set "host" to the identified host, and "is_human" true/false.
4. method_note: does the paper's described method agree with the ENA fields (WGS/shotgun)?
   (agree / mismatch / unclear) -- note if the paper says 16S-only but ENA shows WGS, etc.

Return ONLY:
{{"true_animal_gut":{{"verdict":"...","confidence":"...","evidence":"..."}},
"body_site":{{"verdict":"...","confidence":"...","evidence":"..."}},
"host_id":{{"verdict":"...","confidence":"...","host":"...","is_human":false,"evidence":"..."}},
"method_note":{{"verdict":"...","evidence":"..."}},
"rationale":"one to two sentences on the overall picture"}}"""


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
    print(f"Auditing {len(todo)} studies with {MODEL}...\n")

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

    # ---- flatten + compute keep-decision IN CODE ----
    wgs_by_acc = dict(zip(cat["study_accession"], cat["n_wgs_samples"]))
    src_by_acc = dict(zip(cat["study_accession"], cat["library_source"]))

    rows = []
    for acc, v in verdicts.items():
        def g(crit, f): return (v.get(crit) or {}).get(f, "") if isinstance(v.get(crit), dict) else ""
        gut = g("true_animal_gut", "verdict")
        site = g("body_site", "verdict")
        host = g("host_id", "verdict")
        n_wgs = wgs_by_acc.get(acc, 0) or 0
        src = src_by_acc.get(acc, "")

        # transparent keep-rule
        has_wgs = (n_wgs > 0) or (src == "METAGENOMIC")   # salvage odd-strategy metagenomes
        if gut == "no" or site == "non_gut" or host == "incidental" or not has_wgs:
            keep = "EXCLUDE"
        elif (gut == "yes" and site in ("gut_explicit", "gut_inferred")
              and host in ("single_clear", "genus_only")):
            keep = "KEEP"
        else:
            keep = "REVIEW"

        rows.append({
            "study_accession": acc, "keep": keep,
            "gut": gut, "gut_conf": g("true_animal_gut", "confidence"),
            "body_site": site, "site_conf": g("body_site", "confidence"),
            "host_verdict": host, "host": g("host_id", "host"),
            "is_human": g("host_id", "is_human"),
            "method_note": g("method_note", "verdict"),
            "n_wgs_samples": n_wgs, "library_source": src,
            "gut_evidence": g("true_animal_gut", "evidence"),
            "site_evidence": g("body_site", "evidence"),
            "host_evidence": g("host_id", "evidence"),
            "rationale": v.get("rationale", ""),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_PATH, sep="\t", index=False)

    print(f"\n{'='*52}\nAUDIT COMPLETE\n{'='*52}")
    print("Keep decision:");        print(out["keep"].value_counts().to_string())
    print("\nCriterion 1 gut:");     print(out["gut"].value_counts().to_string())
    print("\nCriterion 2 body_site:"); print(out["body_site"].value_counts().to_string())
    print("\nCriterion 3 host:");    print(out["host_verdict"].value_counts().to_string())
    print("\nmethod agreement:");    print(out["method_note"].value_counts().to_string())
    print(f"\nFlagged human: {(out['is_human'].astype(str).str.lower()=='true').sum()}")
    print(f"EXCLUDE or REVIEW: {(out['keep']!='KEEP').sum()}")
    print(f"Saved -> {OUT_PATH}")
    return out


if __name__ == "__main__":
    run()