"""
build_host_sources_xlsx.py
Turn the host-source analysis into a readable Excel workbook.

Tabs:
  README            -- what each column means and the caveats that change how to read it
  Host sources      -- all 1,044 studies: where the host was found, and what each source named
  Model organisms   -- the 122 lab studies modelling a condition in an animal, with the
                       disease named and the model animal derived
  No host in ENA    -- the 246 NO_HOST_DATA studies, split into rescued vs unknown
  Conflicts         -- deposit says one thing, paper says another (worth a human glance)

Presentation layer: reformats results/host_sources.tsv, changes no values.
"""

import re
import pandas as pd
from datetime import date
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

R = "../results/"
SOURCES = R + "host_sources.tsv"
PATTERNS = R + "study_sample_patterns.tsv"
CATALOG = R + "yes_catalog.tsv"
OUT = R + "host_sources.xlsx"

HEAD_FILL = PatternFill("solid", fgColor="2F5496")
HEAD_FONT = Font(name="Arial", bold=True, color="FFFFFF")
BODY_FONT = Font(name="Arial")
TITLE_FONT = Font(name="Arial", bold=True, size=13)
NOTE_FONT = Font(name="Arial", italic=True, color="666666")

WIDTHS = {
    "study_accession": 16, "has_paper": 11, "paper_title": 60,
    "host_ena_referenced": 19, "host_ena_value": 26,
    "host_title_referenced": 20, "host_title_value": 34,
    "host_abstract_referenced": 23, "host_abstract_value": 40,
    "is_model_organism": 18, "model_organism_disease": 34, "model_animal": 14,
    "human_by_biome_label": 21, "host_source_summary": 22, "host_source_note": 68,
    "sample_host_pattern": 20, "n_host_populated": 17, "n_samples": 11,
    "biome_labels": 30, "recovery": 22,
}

ANIMALS = ["mouse", "mice", "murine", "rat", "zebrafish", "piglet", "pig",
           "chicken", "macaque", "monkey", "ferret", "dog", "rabbit"]
CANON = {"mice": "mouse", "murine": "mouse", "piglet": "pig", "macaque": "monkey"}

GLOSSARY = [
    ("study_accession", "The study. Join key to yes_catalog.tsv and samples.tsv."),
    ("has_paper", "FALSE = deposit-only study. Title/abstract columns read N/A, not blank."),
    ("host_ena_referenced", "Does the ENA deposit give a real host? Biome labels ('gut metagenome') do NOT count."),
    ("host_ena_value", "The ENA host. Falls back to the sample-level host list when host_species is a biome label."),
    ("host_title_referenced / _value", "Host named in the paper TITLE, at the specificity the title used."),
    ("host_abstract_referenced / _value", "Host named in the ABSTRACT. Independent of the title -- the two can disagree."),
    ("is_model_organism", "TRUE when the SEQUENCED SAMPLES came from an animal modelling a condition."),
    ("model_organism_disease", "The condition being modelled (colitis, Alzheimer's...). Blank if unnamed."),
    ("model_animal", "DERIVED here from the host text, for filtering. Not from the model."),
    ("human_by_biome_label", "TRUE = no ENA host, but biome label / abstract shows the host is human."),
    ("host_source_summary", "ENA + paper | paper only | ENA only | human-by-biome-label | unknown."),
    ("host_source_note", "Where the sequenced samples came from, in the model's words."),
    ("sample_host_pattern", "From study_sample_patterns.tsv: UNIFORM / DIVERSE / NO_HOST_DATA etc."),
]

CAVEATS = [
    "'paper only' and 'human-by-biome-label' are the RESCUE cases: no host in the ENA "
    "deposit, but the host is recoverable. 118 of the 246 NO_HOST_DATA studies (48%) "
    "fall here -- see the 'No host in ENA' tab.",

    "A model-organism study is never labelled human-by-biome-label: its host is the "
    "model animal, not the human whose disease is modelled.",

    "The host is the organism the SEQUENCED SAMPLES came from -- not every animal the "
    "paper mentions. Papers that sequence a human cohort and then validate a mechanism "
    "in mice are HUMAN studies; the mice are downstream. host_source_note records which.",

    "Title and abstract are extracted INDEPENDENTLY and can disagree. That is signal, "
    "not noise -- see the 'Conflicts' tab.",

    "Values are recorded at the specificity the paper used. 'ruminant' stays 'ruminant'; "
    "no species is invented, and a general term is never expanded into a list.",

    "model_animal is derived by string-matching the host text in THIS script, for "
    "filtering convenience. model_organism_disease comes from the model.",

    "Nothing here is written back into yes_catalog.tsv, and no sample host is ever "
    "backfilled from a study-level value.",
]


def derive_animal(v):
    s = str(v).lower()
    for a in ANIMALS:
        if re.search(rf"\b{a}", s):
            return CANON.get(a, a)
    return "other/unclear"


def style_row(ws, r, ncols, font=None, fill=None, wrap_col=None):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=r, column=c)
        cell.font = font or BODY_FONT
        if fill:
            cell.fill = fill
        if wrap_col and c == wrap_col:
            cell.alignment = Alignment(wrap_text=True, vertical="top")


def add_table(wb, title, df):
    ws = wb.create_sheet(title)
    cols = list(df.columns)
    ws.append(cols)
    style_row(ws, 1, len(cols), HEAD_FONT, HEAD_FILL)
    for row in df.itertuples(index=False, name=None):
        ws.append(["" if pd.isna(v) else v for v in row])
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{len(df) + 1}"
    for i, name in enumerate(cols, start=1):
        ws.column_dimensions[get_column_letter(i)].width = WIDTHS.get(name, 16)
    return ws


def build_readme(wb, src, nh):
    ws = wb.create_sheet("README")
    ws.column_dimensions["A"].width = 34
    ws.column_dimensions["B"].width = 104

    def line(a="", b="", font=BODY_FONT, bfont=None):
        ws.append([a, b])
        r = ws.max_row
        ws.cell(row=r, column=1).font = font
        ws.cell(row=r, column=2).font = bfont or BODY_FONT
        ws.cell(row=r, column=2).alignment = Alignment(wrap_text=True, vertical="top")

    line("HOST SOURCES", "", TITLE_FONT)
    line("", f"Generated {date.today().isoformat()} by src/extract_host_sources.py")
    line("", "Where each study's host can be found: the ENA deposit, the paper title, "
             "the paper abstract -- and which host each source names.")
    line()
    line("WHERE THE HOST COMES FROM", "", TITLE_FONT)
    for k, v in src.host_source_summary.value_counts().items():
        line(f"  {k}", f"{v} studies")
    line()
    line("THE 246 STUDIES WITH NO ENA HOST", "", TITLE_FONT)
    for k, v in nh.recovery.value_counts().items():
        line(f"  {k}", f"{v} studies")
    line("  ->", f"{int((nh.recovery != 'not recoverable').sum())} of {len(nh)} rescued "
                 f"({100*(nh.recovery != 'not recoverable').mean():.0f}%)")
    line()
    line("TABS", "", TITLE_FONT)
    line("  Host sources", "All 1,044 studies, every column.")
    line("  Model organisms", "The 122 lab studies modelling a condition in an animal.")
    line("  No host in ENA", "The 246 NO_HOST_DATA studies, rescued vs not.")
    line("  Conflicts", "Deposit and paper disagree, or title and abstract disagree.")
    line()
    line("COLUMNS", "", TITLE_FONT)
    for k, v in GLOSSARY:
        line(f"  {k}", v)
    line()
    line("READ THIS BEFORE COUNTING", "", TITLE_FONT)
    for c in CAVEATS:
        line("  -", c, BODY_FONT, NOTE_FONT)
    return ws


def run():
    src = pd.read_csv(SOURCES, sep="\t")
    pat = pd.read_csv(PATTERNS, sep="\t")[["study_accession", "n_samples"]]
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)[["study_accession", "paper_title"]]

    src = src.merge(cat, on="study_accession", how="left").merge(pat, on="study_accession", how="left")
    src["model_animal"] = src.apply(
        lambda r: derive_animal(r.host_abstract_value if pd.notna(r.host_abstract_value)
                                else r.host_title_value) if r.is_model_organism else "", axis=1)

    main_cols = ["study_accession", "paper_title", "has_paper",
                 "host_source_summary", "host_ena_referenced", "host_ena_value",
                 "host_title_referenced", "host_title_value",
                 "host_abstract_referenced", "host_abstract_value",
                 "is_model_organism", "model_organism_disease", "model_animal",
                 "human_by_biome_label", "sample_host_pattern", "n_samples",
                 "host_source_note", "biome_labels"]

    # ---- model organisms ----
    mo = src[src.is_model_organism].sort_values(
        ["model_animal", "model_organism_disease"])[
        ["study_accession", "model_organism_disease", "model_animal",
         "host_abstract_value", "host_ena_value", "paper_title",
         "n_samples", "host_source_note"]]

    # ---- no ENA host ----
    nh = src[src.sample_host_pattern == "NO_HOST_DATA"].copy()
    nh["recovery"] = nh.host_source_summary.map({
        "paper only": "rescued: paper names host",
        "human-by-biome-label": "rescued: human, implicit",
    }).fillna("not recoverable")
    nh = nh.sort_values(["recovery", "n_samples"], ascending=[True, False])[
        ["study_accession", "recovery", "host_abstract_value", "host_title_value",
         "biome_labels", "n_samples", "paper_title", "host_source_note"]]

    # ---- conflicts ----
    def disagree(r):
        t, a = str(r.host_title_value).lower(), str(r.host_abstract_value).lower()
        both = r.host_title_referenced and r.host_abstract_referenced
        return bool(both and t not in a and a not in t)

    conf = src[(src.is_model_organism & (src.host_ena_value == "Homo sapiens")) |
               src.apply(disagree, axis=1)].copy()
    conf["conflict_type"] = conf.apply(
        lambda r: "ENA=human, paper=animal model"
        if (r.is_model_organism and r.host_ena_value == "Homo sapiens")
        else "title vs abstract differ", axis=1)
    conf = conf.sort_values("conflict_type")[
        ["study_accession", "conflict_type", "host_ena_value", "host_title_value",
         "host_abstract_value", "model_organism_disease", "paper_title"]]

    wb = Workbook()
    wb.remove(wb.active)
    build_readme(wb, src, nh)
    add_table(wb, "Host sources", src[main_cols])
    add_table(wb, "Model organisms", mo)
    add_table(wb, "No host in ENA", nh)
    add_table(wb, "Conflicts", conf)
    wb.save(OUT)

    print(f"Host sources : {len(src)}")
    print(f"Model organisms: {len(mo)}")
    print(f"No host in ENA : {len(nh)}  (rescued {int((nh.recovery != 'not recoverable').sum())})")
    print(f"Conflicts      : {len(conf)}")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    run()
