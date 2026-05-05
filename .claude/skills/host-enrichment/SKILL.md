---
name: host-enrichment
description: Derive the host species (single scientific binomial) for a list of ENA/NCBI study accession codes by fetching each study's title and PubMed abstract, then applying a regex + Claude Haiku hybrid classifier. Use when the user provides one or more study accessions (PRJNA*, PRJEB*, PRJDB*, SRP*, ERP*, DRP*) and asks to derive, enrich, or identify their host species.
---

This skill enriches a list of study accession codes with derived host species information by reusing the project's pipeline at `src/host_enrichment.py`.

## Inputs

The user provides one or more study accession codes. Accept any of these shapes:
- Comma-separated: `PRJNA801645, PRJEB47613`
- Space-separated: `PRJNA801645 PRJEB47613`
- One per line
- Embedded in prose ("can you enrich PRJNA801645 and PRJEB47613?")

Recognised accession prefixes: `PRJNA`, `PRJEB`, `PRJDB`, `SRP`, `ERP`, `DRP`.

If the user references a file (e.g. "the accessions in foo.txt"), read the file and extract accessions from it.

## How to run

1. Extract accession codes from the user's message into a list. Deduplicate, preserve order.
2. Pick an output path:
   - If the user specifies one, use it.
   - Otherwise default to `results/host_enriched_<N>studies.csv` where `<N>` is the count.
3. Invoke the enrichment from the project root:

   ```bash
   python -m src.host_enrichment --output <output_path> <acc1> <acc2> ...
   ```

   The CLI also accepts a single comma-separated string, e.g. `python -m src.host_enrichment --output out.csv "PRJNA801645,PRJEB47613"`.

   Resume is automatic: if the output CSV exists, already-processed accessions are skipped.

4. After it finishes, read the output CSV and report:
   - The path written
   - A breakdown of `derivation_source` counts (how many came from regex vs LLM)
   - Any rows where `derived_host_species` is empty (`llm_unknown` or `llm_error`) — these warrant manual review
   - The `derived_host_species` column for each accession

## Pipeline summary (for context only — do not re-implement)

For each accession the pipeline:
1. Calls `fetch_study_origin(acc)` → study title (ENA, NCBI BioProject fallback).
2. Calls `fetch_pubmed_abstract(acc)` → PubMed abstract (or `None` if no linked paper).
3. **Regex layer**: scans both texts for ~27 curated common-name + binomial patterns mapped to canonical scientific names. Confirms when title ∩ abstract = exactly one species, or when abstract is missing and title has exactly one species.
4. **LLM fallback** (`claude-haiku-4-5`, cached system prompt): if regex can't confirm, asks the model for the single primary host binomial or `UNKNOWN`.

## Output schema

| column | meaning |
|---|---|
| `study_accession` | input accession |
| `title` | fetched study title |
| `abstract` | fetched PubMed abstract (may be empty) |
| `derived_host_species` | single scientific binomial, or empty if unidentifiable |
| `derivation_source` | one of: `regex_confirmed`, `regex_title_only_no_abstract`, `llm`, `llm_unknown`, `llm_error` |

## Requirements

- Run from the project root (`/Users/ananyakharya/Documents/study-discovery-engine`) so `src/` is importable.
- `.env` must contain `NCBI_EMAIL`, `NCBI_API_KEY`, and `ANTHROPIC_API_KEY`.
- NCBI calls are rate-limited; expect ~0.5–1s per accession.

## Examples

**User says:** "enrich PRJNA801645 and PRJEB47613"
→ Run: `python -m src.host_enrichment --output results/host_enriched_2studies.csv PRJNA801645 PRJEB47613`

**User says:** "derive host species for these: PRJNA801645, PRJDB10675, SRP123456, save to /tmp/out.csv"
→ Run: `python -m src.host_enrichment --output /tmp/out.csv PRJNA801645 PRJDB10675 SRP123456`

**User says:** "the accessions are in mylist.txt"
→ Read `mylist.txt`, extract accessions, then run the command above with the parsed list.
