"""
build_catalog_xlsx.py
Turn yes_catalog.tsv (1,044 studies, 55 columns across all enrichment layers) into a
readable Excel workbook.

55 columns is too many to scan, so the main tab leads with the columns people actually
use and keeps the rest to the right, in layer order. Nothing is dropped.

Tabs:
  README          -- column glossary BY LAYER, plus the caveats that change the counts
  Catalog         -- all 1,044 studies, all 55 columns
  Host sources    -- the host-provenance view: where each study's host was found
  Model organisms -- the 122 lab studies modelling a condition in an animal
  No host in ENA  -- studies with no sample-level ENA host, rescued vs not

Presentation layer: reformats the TSV, changes no values.
"""

import pandas as pd
from datetime import date
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

R = "../results/"
CATALOG = R + "yes_catalog.tsv"
OUT = R + "yes_catalog.xlsx"

HEAD_FILL = PatternFill("solid", fgColor="2F5496")
HEAD_FONT = Font(name="Arial", bold=True, color="FFFFFF")
BODY_FONT = Font(name="Arial")
TITLE_FONT = Font(name="Arial", bold=True, size=13)
NOTE_FONT = Font(name="Arial", italic=True, color="666666")

# the columns worth seeing first, in this order; everything else follows
LEAD = ["study_accession", "host_species", "host_class", "body_site", "n_samples",
        "n_total_samples", "n_metagenomic_samples",
        "library_strategy", "data_source", "pub_year",
        "paper_title", "pmid",
        "host_source_summary", "host_ena_value", "host_title_value",
        "host_abstract_value", "is_model_organism", "model_organism_disease",
        "human_by_biome_label", "sample_host_pattern",
        "keep", "is_human", "gut", "discrepancy_flags",
        "n_ena_standard_fields", "n_ena_custom_fields"]

WIDE = {"ena_standard_fields": 110, "ena_custom_fields": 90,
        "n_ena_standard_fields": 20, "n_ena_custom_fields": 19, "llm_notes": 70, "method_note": 46, "discrepancy_notes": 60,
        "discrepancy_host": 52, "discrepancy_host_missing": 46,
        "discrepancy_gut": 56, "discrepancy_multihost": 56,
        "host_source_note": 66, "paper_link_suspect_reason": 24,
        "discrepancy_flags": 24, "host_verdict": 18, "audit_body_site": 18,
        "audit_host": 24, "keep": 10, "gut": 9, "n_metagenomic_samples": 23, "library_source_mix": 30, "paper_title": 62, "abstract": 70, "host_source_note": 66, "llm_evidence": 50,
        "host_species_list": 44, "sample_titles_cache": 44, "isolation_sources": 38,
        "biome_labels": 32, "body_sites": 30, "library_strategy_mix": 26,
        "host_abstract_value": 40, "host_title_value": 34, "host_ena_value": 26,
        "model_organism_disease": 34, "host_source_summary": 22,
        "sample_host_pattern": 20, "fulltext_link": 34, "doi": 30, "journal": 30}

LAYERS = [
    ("1. Original pipeline", "study_accession, host_species, host_tax_id, body_site, "
     "country, n_samples, library_strategy, library_source, data_source, pmid, "
     "classification, needs_review"),
    ("2. ENA enrichment", "library_strategy_mix, n_wgs_samples, n_host_species, "
     "host_species_list, body_sites, biome_labels, isolation_sources, "
     "sample_titles_cache, n_total_samples"),
    ("3. Paper metadata", "paper_title, abstract, doi, pub_year, journal, first_author, "
     "last_author, fulltext_link  (561/1044 -- the PMID ceiling)"),
    ("4. Taxonomy", "host_kingdom, host_phylum, host_class, host_order, host_family, "
     "host_genus, host_species_taxonomy"),
    ("5. LLM classification", "classification_final, llm_confidence, llm_body_site_signal, "
     "llm_evidence"),
    ("6. Audit (NEW)", "keep (KEEP/REVIEW/EXCLUDE), gut, gut_conf, audit_body_site, "
     "host_verdict, audit_host, is_human, method_note, llm_notes -- all 1,044 studies"),
    ("7. Discrepancy (NEW)", "discrepancy_flags, discrepancy_host, "
     "discrepancy_host_missing, discrepancy_gut, discrepancy_multihost, "
     "discrepancy_notes -- only the 560 paper-linked studies could be checked"),
    ("8. Sample-derived (NEW)", "library_source_mix, n_metagenomic_samples -- "
     "recomputed from samples.tsv; library_strategy_mix backfilled where it was blank"),
    ("9. ENA fields (NEW)", "ena_standard_fields + ena_custom_fields (with their counts) "
     "-- TOGETHER these are every ENA field the study used. Standard = ENA's normalized "
     "vocabulary from the Portal API; custom = submitter tags with no standard twin, read "
     "verbatim from the sample XML. Clean partition: no field appears in both."),
    ("10. Host sources (NEW)", "has_paper, host_ena_referenced/_value, "
     "host_title_referenced/_value, host_abstract_referenced/_value, is_model_organism, "
     "model_organism_disease, human_by_biome_label, host_source_summary, "
     "host_source_note, sample_host_pattern, n_host_populated"),
]

GLOSSARY = [
    ("keep", "Audit verdict: KEEP 741 / REVIEW 272 / EXCLUDE 31. REVIEW means 'needs a "
     "glance' (thin metadata, WGA method), NOT 'problem'."),
    ("is_human", "211 studies. Humans ARE animals here: kept and flagged, never silently "
     "dropped. The benchmark excludes them, so report both totals."),
    ("gut / gut_conf", "Is this gut-derived, and how confident: yes 834 / uncertain 175 / no 27."),
    ("llm_notes", "The audit's per-study evidence, for ALL 1,044 studies. The llm_evidence "
     "column covers only the 177 sweep-classified studies -- different pass, not a gap."),
    ("discrepancy_flags", "NOT_CHECKED (484, no paper to compare) | none (390, checked and "
     "clean) | host_mismatch / host_missing / gut_sample_issue / multihost_issue (170)."),
    ("ena_standard_fields", "Every STANDARDIZED ENA field the study populated, across "
     "sample/run/experiment/study level. Median 47 fields."),
    ("ena_custom_fields", "Every SUBMITTER-SPECIFIC tag with no standardized twin, read "
     "verbatim from the sample XML (ref_biomaterial, source_material_id, replicate, "
     "breed, ngdc_sample_id...). 676 of 1,044 studies use at least one."),
    ("GRAMMAR (both columns)", "' || ' between fields; ' = ' key-value (split on the "
     "FIRST occurrence); '; ' between distinct values; '[N distinct]' above 8 values; "
     "'...[truncated]' for a single value over 200 chars. Sentinels ('missing', 'not "
     "applicable') shown literally; unused fields omitted. See the 'ENA fields' tab."),
    ("host_source_summary", "Where the host was found: ENA + paper | paper only | "
     "ENA only | human-by-biome-label | unknown."),
    ("host_ena_value", "Host from the ENA deposit. Biome labels ('gut metagenome') do "
     "NOT count as a host."),
    ("host_title_value / host_abstract_value", "Host named in the paper title and in the "
     "abstract, extracted INDEPENDENTLY -- they can disagree, and that is signal."),
    ("is_model_organism", "TRUE when the SEQUENCED SAMPLES came from an animal modelling "
     "a condition -- not merely because the paper mentions mice."),
    ("model_organism_disease", "The condition modelled (colitis, Alzheimer's, ...)."),
    ("human_by_biome_label", "TRUE = no ENA host, but the biome label / abstract shows "
     "the host is human. Turns 'no host' into 'host is human, unstated'."),
    ("sample_host_pattern", "From the sample-level catalog: UNIFORM / UNIFORM_WITH_GAPS "
     "/ DIVERSE / DIVERSE_AND_SPARSE / NO_HOST_DATA."),
    ("n_host_populated", "How many of the study's biosamples carry a real ENA host."),
    ("library_source_mix", "Per-sample library_source counts. The AUTHORITATIVE axis: "
     "library_source beats library_strategy."),
    ("n_metagenomic_samples", "Samples with library_source=METAGENOMIC -- what "
     "n_wgs_samples was meant to count. See the caveat below."),
]

CAVEATS = [
    "host_species (layer 1) is the STUDY-level host and is often a mode across samples, "
    "or a biome label like 'feces metagenome'. host_ena_value excludes biome labels, so "
    "the two legitimately disagree.",

    "48% of studies with no sample-level ENA host (118 of 246) are still recoverable: "
    "64 from the paper, 54 human-by-biome-label. See the 'No host in ENA' tab.",

    "A model-organism study is never human-by-biome-label -- its host is the model "
    "animal, not the human whose disease is modelled.",

    "The host is the organism the SEQUENCED SAMPLES came from. Papers that sequence a "
    "human cohort and then validate in mice are HUMAN studies; host_source_note says which.",

    "Humans ARE animals here: kept and flagged, never silently excluded. The "
    "AnimalMetagenome DB benchmark excludes them, so report both totals.",

    "Host values keep the paper's own specificity: 'ruminant' stays 'ruminant'. No "
    "species is invented and a general term is never expanded into a list.",

    "n_wgs_samples is WRONG for 234 studies and is kept only for audit. Its formula "
    "tested library_STRATEGY for the value 'METAGENOMIC', which is a library_SOURCE "
    "value and never appears as a strategy, so the source-beats-strategy rule never "
    "fired. Use n_metagenomic_samples: 136,058 samples vs the old 101,334.",

    "6 studies have no METAGENOMIC samples at all (OTHER/WGA/Hi-C, library_source "
    "GENOMIC or OTHER, 1-6 samples each). Their zero is real, not a counting bug.",

    "13 studies have an existing library_strategy_mix that disagrees with a fresh "
    "count from samples.tsv -- the deposits grew after the original enrichment ran. "
    "Existing values were left as they were, not overwritten.",

    "The field separator is ' || ', NOT a single '|'. ENA uses a bare '|' INSIDE values "
    "as its own multi-value separator, so a single pipe is unsafe to split on.",

    "Values over 200 chars are truncated FOR DISPLAY. Two values differing only past the "
    "cutoff render identically, but the distinct-value count reflects the FULL values -- "
    "that is not a duplicate.",

    "ena_custom_fields is empty for 368 studies: they used only standardized fields. It is "
    "also empty for 16 studies whose samples ENA has no XML for (recent NCBI BioSamples "
    "not yet mirrored) -- see results/ena_refetch.tsv, where a bounded 40-sample probe per "
    "study recorded whether the check was exhaustive.",

    "Custom tag spelling variants are NOT merged: 'host_subject_id' (55 studies) and "
    "'host subject id' (32) stay distinct, faithful to what submitters wrote.",

    "discrepancy_flags = NOT_CHECKED means the study has no linked paper, so paper-vs-"
    "deposit could not be compared. It does NOT mean the study is clean -- do not count "
    "those 484 as passing.",

    "llm_confidence / llm_body_site_signal / llm_evidence are populated for only 177 "
    "studies, and that is by design: they come from the PubMed-sweep classifier, which "
    "ran on a disjoint set from the 867 carrying a 'classification' value. For per-study "
    "evidence covering everything, use llm_notes from the audit layer.",

    "biome_labels is the catalog's original column; the host-source merge deliberately "
    "did not overwrite it.",
]


def style_row(ws, r, ncols, font=None, fill=None):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=r, column=c)
        cell.font = font or BODY_FONT
        if fill:
            cell.fill = fill


def add_table(wb, title, df):
    ws = wb.create_sheet(title)
    cols = list(df.columns)
    ws.append(cols)
    style_row(ws, 1, len(cols), HEAD_FONT, HEAD_FILL)
    for row in df.itertuples(index=False, name=None):
        ws.append(["" if pd.isna(v) else v for v in row])
    ws.freeze_panes = "B2"          # keep study_accession visible while scrolling right
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{len(df) + 1}"
    for i, name in enumerate(cols, start=1):
        ws.column_dimensions[get_column_letter(i)].width = WIDE.get(name, 17)
    return ws


# columns whose content is prose, not a value -- wrap them instead of letting them
# spill across their neighbours
WRAP = {"ena_standard_fields", "ena_custom_fields", "discrepancy_host", "discrepancy_host_missing", "discrepancy_gut",
        "discrepancy_multihost", "discrepancy_notes", "llm_notes", "method_note",
        "host_source_note", "paper_title"}


def add_wrapped_table(wb, title, df):
    """Like add_table, but prose columns wrap and rows top-align, so a 400-character
    note reads as a paragraph in its own cell rather than a line across the sheet."""
    ws = add_table(wb, title, df)
    cols = list(df.columns)
    wrap_idx = [i for i, c in enumerate(cols, start=1) if c in WRAP]
    for r in range(2, len(df) + 2):
        for i in wrap_idx:
            ws.cell(row=r, column=i).alignment = Alignment(wrap_text=True,
                                                           vertical="top")
    return ws


def build_readme(wb, cat):
    ws = wb.create_sheet("README")
    ws.column_dimensions["A"].width = 36
    ws.column_dimensions["B"].width = 102

    def line(a="", b="", font=BODY_FONT, bfont=None):
        ws.append([a, b])
        r = ws.max_row
        ws.cell(row=r, column=1).font = font
        ws.cell(row=r, column=2).font = bfont or BODY_FONT
        ws.cell(row=r, column=2).alignment = Alignment(wrap_text=True, vertical="top")

    line("ANIMAL GUT WGS METAGENOMICS CATALOG", "", TITLE_FONT)
    line("", f"Generated {date.today().isoformat()} from results/yes_catalog.tsv")
    line("", f"{len(cat):,} studies x {len(cat.columns)} columns")
    line()
    line("WHERE THE HOST COMES FROM", "", TITLE_FONT)
    for k, v in cat.host_source_summary.value_counts().items():
        line(f"  {k}", f"{v} studies")
    line()
    line("SAMPLE-LEVEL HOST PATTERN", "", TITLE_FONT)
    for k, v in cat.sample_host_pattern.value_counts().items():
        line(f"  {k}", f"{v} studies")
    line()
    line("COLUMN LAYERS", "", TITLE_FONT)
    for name, cols in LAYERS:
        line(f"  {name}", cols)
    line()
    line("KEY COLUMNS", "", TITLE_FONT)
    for k, v in GLOSSARY:
        line(f"  {k}", v)
    line()
    line("READ THIS BEFORE COUNTING", "", TITLE_FONT)
    for c in CAVEATS:
        line("  -", c, BODY_FONT, NOTE_FONT)
    return ws


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    ordered = [c for c in LEAD if c in cat.columns] + \
              [c for c in cat.columns if c not in LEAD]
    # keep each corrected column NEXT TO the one it corrects -- a stale n_wgs_samples
    # sitting 15 columns away from its fix reads as a bug rather than as history
    for corrected, anchor in [("library_source_mix", "library_strategy_mix"),
                              ("n_metagenomic_samples_dup", "n_wgs_samples")]:
        if corrected in ordered and anchor in ordered:
            ordered.remove(corrected)
            ordered.insert(ordered.index(anchor) + 1, corrected)
    cat = cat[[c for c in ordered if c in cat.columns]]

    hs = cat[["study_accession", "paper_title", "host_source_summary",
              "host_ena_value", "host_title_value", "host_abstract_value",
              "is_model_organism", "model_organism_disease", "human_by_biome_label",
              "sample_host_pattern", "n_host_populated", "host_source_note"]]

    mo = cat[cat.is_model_organism == True][
        ["study_accession", "model_organism_disease", "host_abstract_value",
         "host_ena_value", "host_species", "paper_title", "n_samples",
         "host_source_note"]].sort_values("model_organism_disease")

    nh = cat[cat.sample_host_pattern == "NO_HOST_DATA"].copy()
    nh["recovery"] = nh.host_source_summary.map({
        "paper only": "rescued: paper names host",
        "human-by-biome-label": "rescued: human, implicit",
    }).fillna("not recoverable")
    nh = nh.sort_values(["recovery", "n_total_samples"], ascending=[True, False])[
        ["study_accession", "recovery", "host_abstract_value", "host_title_value",
         "biome_labels", "n_total_samples", "paper_title", "host_source_note"]]

    wb = Workbook()
    wb.remove(wb.active)
    build_readme(wb, cat)
    add_table(wb, "Catalog", cat)
    add_table(wb, "Host sources", hs)
    add_table(wb, "Model organisms", mo)
    add_table(wb, "No host in ENA", nh)

    # one row per FLAGGED study, evidence side by side and readable
    dis = cat[~cat.discrepancy_flags.isin(["none", "NOT_CHECKED"])][
        ["study_accession", "discrepancy_flags", "discrepancy_host",
         "discrepancy_host_missing", "discrepancy_gut", "discrepancy_multihost",
         "discrepancy_notes", "paper_title"]].sort_values("discrepancy_flags")
    add_wrapped_table(wb, "Discrepancies", dis)
    print(f"Discrepancies  : {len(dis)} flagged studies")

    # the summary is ~1,900 chars on average -- it needs a wrapped tab of its own
    ena = cat[["study_accession", "n_ena_standard_fields", "n_ena_custom_fields",
               "ena_custom_fields", "ena_standard_fields"]].sort_values(
                   "n_ena_custom_fields", ascending=False)
    add_wrapped_table(wb, "ENA fields", ena)
    print(f"ENA fields     : {len(ena)} studies")
    wb.save(OUT)

    print(f"Catalog        : {len(cat)} studies x {len(cat.columns)} columns")
    print(f"Host sources   : {len(hs)}")
    print(f"Model organisms: {len(mo)}")
    print(f"No host in ENA : {len(nh)}  (rescued "
          f"{int((nh.recovery != 'not recoverable').sum())})")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    run()
