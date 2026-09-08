# study-discovery-engine

Automated, AI-enriched catalog of **animal gut WGS (shotgun) metagenomics studies**
across the animal kingdom. Undergraduate research project in the Knight Lab
(mentor: Sam). The catalog is the input to Sam's tree-of-life resubmission.

Benchmark being matched/beaten: **AnimalMetagenome DB** (PMC9203544) — 1,044 studies,
10,672 WGS samples, manual, human-excluded, through May 2021.

## Environment
- Conda env: **study-triage** (activate before running anything)
- Notebooks in `notebooks/` (main: `06_build_catalog.ipynb`), scripts in `src/`,
  results in `../results/` (relative to notebooks/ or src/)
- Fresh kernel needs: `import sys, os; sys.path.insert(0, os.path.abspath(".."))`
- NCBI Entrez: configure via `src/fetcher.py` `configure_entrez()`, email
  akharya@ucsd.edu, needs ~0.34s sleep between calls
- Anthropic API: key in `ANTHROPIC_API_KEY` env var. Model used for
  classification/judgment: **claude-haiku-4-5-20251001** ($1/$5 per M tokens).
  Free-tier rate limit ~5/min → sequential calls + backoff. Personal API
  credit (lab account still pending — ask Sam).

## Hard-won gotchas (do NOT relearn these the hard way)
- **Confirm files are actually SAVED into `src/` before running** — many failures
  traced to editing a file that was never saved.
- `importlib.reload` is UNRELIABLE (stale-module cache). To force-refresh:
  `del sys.modules["src.module_name"]` then re-import, or restart kernel.
- Save DataFrames to disk immediately — vars are lost on kernel restart.
- All long fetch/LLM passes must be **checkpointed + resumable** (write a JSON
  checkpoint every N items, skip already-done on rerun).
- Watch for **column-name collisions** on merge — catalog has grown many columns
  across layers (original pipeline + enrichment + audit + discrepancy). Drop
  overlapping cols before merge, or they become `_x`/`_y` and break selection.
  Known trio: `body_site` (original), `body_sites` (enrichment), `llm_body_site_signal`.
- Empty-string truthiness: Python `and`/`or` return operands, not bools. Wrap
  masks in `bool()` when building boolean columns via `.apply`.

## Data: results/yes_catalog.tsv (1,044 studies, ~107k samples)
One row per study. Column layers (each merged in by study_accession):
1. Original pipeline: study_accession, data_source, host_species, body_site, pmid...
2. Enrichment (`enrich_catalog_from_ena.py`): library_strategy_mix, n_wgs_samples,
   n_host_species, host_species_list, body_sites, biome_labels, isolation_sources,
   sample_titles_cache, n_total_samples
3. Paper metadata (`backfill_pmid_metadata.py`, `backfill_abstracts.py`):
   paper_title, abstract, doi, pub_year, journal, first_author, last_author.
   Coverage = 561/1044 (PMID ceiling); abstract coverage 561/561 of those.
4. Audit (`audit_catalog_v3.py` → catalog_audit.tsv): keep (KEEP/REVIEW/EXCLUDE),
   is_human, gut, body_site, host_verdict, method_note, llm_notes
5. Discrepancy (`detect_discrepancies.py` → discrepancies.tsv): discrepancy_flags,
   discrepancy_host, discrepancy_host_missing, discrepancy_gut,
   discrepancy_multihost, discrepancy_notes

## Domain rules (these encode real biology — apply consistently)
- **library_source beats library_strategy.** Studies with odd strategy labels
  (OTHER/WGA/Hi-C) but library_source=METAGENOMIC are legit metagenomes.
- **Humans ARE animals: keep + flag** (`is_human`), do not silently exclude.
  Benchmark excludes humans, so report both total AND animal-only numbers.
- **Unclear/multi/genus-only host = KEEP** (a metadata gap is not invalidity).
- **Model organisms:** gut-microbiome papers about a HUMAN disease are often done
  in a MOUSE/RAT model. Mouse host + human-disease paper is NOT a host mismatch.
  (Known remaining false-positive: mouse+human still trips the multihost flag.)
- **General vs specific = agreement**, never a conflict ("ruminant" vs buffalo;
  mouse vs Mus musculus).
- **Failures must be visible, never silent-clean.** A parse-error/never-checked
  study must be marked distinctly (e.g. NOT_CHECKED), not counted as passing.

## Audit result: KEEP 741, REVIEW 272, EXCLUDE 31, human-flagged 207.
REVIEW = "needs a glance" (WGA-method, thin metadata), not "problem".

## Discrepancy result (cross-check paper vs ENA deposit, 560 paper-linked studies):
170 flagged — host_mismatch 33, host_missing 86, gut_sample_issue 24, multihost 60.
Two independent sources = PubMed title+abstract (paper) vs ENA metadata (deposit).
Adjudication: deposit sample-metadata wins body-site; paper wins host identity.

## Open threads / next steps
- **Sample-level catalog (Sam's current ask):** pull ALL samples per study from
  ENA, then have Claude parse per-study sample patterns — are samples all same
  host (uniform, maybe wrong) vs all different (diverse, maybe meticulous, e.g.
  PRJNA678871), find multi-seq-type studies. Two-phase: (1) real ENA bulk fetch →
  samples.tsv (one row/sample) + pandas pattern stats, no LLM; (2) Claude API
  interprets ambiguous patterns. Fetch must be a BULK query (not per-study loop —
  those are slow). Claude cannot fetch data itself — never let it hallucinate
  sample rows; the fetch runs against real ENA.
- Merge discrepancy columns into yes_catalog.tsv (distinct names, watch collisions)
- 483 deposit-only studies: deposit-internal consistency check (awaiting Sam)
- Human scope decision (awaiting Sam)
- Taxonomy cleanup: standardize Actinopteri/Actinopterygii; fix plant FPs
- Self-updating pipeline: SQLite → incremental date-filtered refresh → diff