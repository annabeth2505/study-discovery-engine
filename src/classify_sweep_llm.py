"""
classify_sweep_llm.py
Re-classify the 64 host-sweep studies through the SAME LLM classifier used for
the PubMed batches, so every study in the catalog is judged by one standard.

Why: the sweep's original verdicts (biome-label + body-site logic) let multi-host
mislabels through -- e.g. a "Turtles" study that is actually the American Gut
Project (human), or a "Lizards" study that is actually about butterflies. The LLM
reads the paper title + full host list together and catches these.

Evidence comes straight from sweep_enriched.tsv (already has host, biome_labels,
body_sites, paper_title, abstract) -- no re-fetch needed for classification.
(host_tax_id for taxonomy is fetched later, only for the confirmed ones.)

Output: ../results/sweep_llm_classified.tsv
"""

import os
import json
import time
import pandas as pd
import anthropic

MODEL = "claude-haiku-4-5-20251001"
SWEEP_PATH = "../results/sweep_enriched.tsv"
CHECKPOINT = "../results/sweep_llm_verdicts.json"
OUT_PATH = "../results/sweep_llm_classified.tsv"

client = anthropic.Anthropic()

SYSTEM = (
    "You are a metagenomics data curator. You decide whether a study is a SHOTGUN/WGS "
    "METAGENOMIC study of an ANIMAL's GUT microbiome (gut, intestinal, fecal/faecal, "
    "cecal, rumen, hindgut, digesta; cloacal/rectal count as gut-proxy). "
    "IMPORTANT: many of these studies list MANY hosts. Decide based on what the study is "
    "PRIMARILY about, judged from the paper title and the dominant host(s). If the study "
    "is a large multi-host survey where the target animal is incidental (e.g. a human gut "
    "project that happens to include one reptile sample), that is NOT a gut study of that "
    "animal -- return NOT_GUT or NOT_TARGET. Base decisions ONLY on the evidence. If body "
    "site/method is inferred rather than stated, set confidence='low'. Return ONLY JSON."
)

PROMPT = """Study {acc} (originally tagged as clade: {clade}):
- Paper title: {paper_title}
- Abstract: {abstract}
- Host organism(s) listed: {host}
- Biome label(s): {biome_labels}
- Body site(s): {body_sites}
- Flagged as multi-host survey: {multihost}

Classify whether this is a shotgun/WGS metagenomic study of an ANIMAL GUT microbiome
for the tagged clade. Watch for multi-host studies where the tagged animal is incidental.

- verdict: GUT_WGS (gut metagenomics, tagged animal is a genuine subject)
           NOT_GUT (metagenomics but non-gut site)
           NOT_TARGET (gut study but the tagged animal is incidental / a trace host in a
                       study primarily about something else, e.g. human gut project)
           NOT_METAGENOMIC (genome, transcriptome, or 16S/amplicon)
           UNCERTAIN
- confidence: high | medium | low
- body_site_signal: explicit | inferred | none
- evidence: quote the specific evidence that drove the decision
- reasoning: one sentence

Return ONLY:
{{"verdict":"...","confidence":"...","body_site_signal":"...","evidence":"...","reasoning":"..."}}"""


def call_llm(prompt, max_retries=6):
    for attempt in range(max_retries):
        try:
            msg = client.messages.create(
                model=MODEL, max_tokens=300, system=SYSTEM,
                messages=[{"role": "user", "content": prompt}])
            text = msg.content[0].text.strip().replace("```json", "").replace("```", "").strip()
            return json.loads(text)
        except anthropic.RateLimitError:
            time.sleep(13 * (attempt + 1))
        except json.JSONDecodeError:
            return {"verdict": "PARSE_ERROR", "confidence": "", "body_site_signal": "",
                    "evidence": "", "reasoning": "non-JSON"}
        except Exception as e:
            return {"verdict": "ERROR", "confidence": "", "body_site_signal": "",
                    "evidence": "", "reasoning": str(e)[:120]}
    return {"verdict": "ERROR", "confidence": "", "body_site_signal": "",
            "evidence": "", "reasoning": "max retries"}


def run():
    sweep = pd.read_csv(SWEEP_PATH, sep="\t")
    print(f"Sweep studies to re-classify: {len(sweep)}")

    verdicts = {}
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            verdicts = json.load(f)
        print(f"Resuming — {len(verdicts)} done")

    todo = sweep[~sweep["study_accession"].isin(verdicts.keys())]
    print(f"Classifying {len(todo)}...\n")

    for i, r in enumerate(todo.itertuples()):
        if i % 10 == 0:
            print(f"  {i}/{len(todo)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(verdicts, f)

        def trunc(x, n):  # keep prompts bounded — host lists can be huge
            s = str(x) if pd.notna(x) else "(none)"
            return s[:n] + "..." if len(s) > n else s

        prompt = PROMPT.format(
            acc=r.study_accession,
            clade=r.clade,
            paper_title=trunc(r.paper_title, 300),
            abstract=trunc(r.abstract, 600),
            host=trunc(r.host, 500),
            biome_labels=trunc(r.biome_labels, 300),
            body_sites=trunc(r.body_sites, 200),
            multihost=getattr(r, "is_multihost_survey", False),
        )
        v = call_llm(prompt)
        v["clade"] = r.clade
        v["paper_title"] = r.paper_title
        v["old_verdict"] = r.gut_verdict_final
        verdicts[r.study_accession] = v
        time.sleep(0.3)

    with open(CHECKPOINT, "w") as f:
        json.dump(verdicts, f)

    out = pd.DataFrame([{"study_accession": a, **v} for a, v in verdicts.items()])

    def tier(r):
        if r["verdict"] == "GUT_WGS" and r["confidence"] in ("high", "medium"):
            return "CONFIRMED_GUT"
        if r["verdict"] == "GUT_WGS":
            return "LIKELY_GUT_REVIEW"
        if r["verdict"] == "UNCERTAIN":
            return "REVIEW"
        return "EXCLUDED"

    out["tier"] = out.apply(tier, axis=1)
    out.to_csv(OUT_PATH, sep="\t", index=False)

    print(f"\n{'='*52}\nSWEEP RE-CLASSIFICATION (LLM)\n{'='*52}")
    print(out["tier"].value_counts().to_string())
    print(f"\nverdict x confidence:")
    print(pd.crosstab(out["verdict"], out["confidence"]).to_string())

    # show how the new verdict compares to the old one
    print(f"\nOLD verdict -> NEW tier (did the LLM change calls?):")
    print(pd.crosstab(out["old_verdict"], out["tier"]).to_string())

    conf = out[out["tier"] == "CONFIRMED_GUT"]
    print(f"\nCONFIRMED GUT after LLM: {len(conf)}")
    print(conf["clade"].value_counts().to_string())

    print(f"\nNewly EXCLUDED that old logic had KEPT (caught mislabels):")
    caught = out[(out["tier"] == "EXCLUDED") & (out["old_verdict"] == "KEEP")]
    for _, r in caught.iterrows():
        print(f"  {r['study_accession']} [{r['clade']}] -> {r['verdict']}: {r['reasoning']}")

    print(f"\nSaved -> {OUT_PATH}")
    return out


if __name__ == "__main__":
    run()