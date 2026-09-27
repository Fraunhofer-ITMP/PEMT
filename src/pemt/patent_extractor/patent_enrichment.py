# -*- coding: utf-8 -*-

"""Script for extracting patent literature from the SureChEMBL bulk data."""

import json
import logging
import os
from typing import Iterable, Optional

import pandas as pd

from pemt.constants import PATENT_DIR, VALID_CODES
from pemt.surechembl import SECTIONS, SureChEMBLBulk, resolve_source

logger = logging.getLogger(__name__)

os.makedirs(PATENT_DIR, exist_ok=True)

#: Columns of ``<analysis_name>_patent_data.tsv``. The first six match earlier PEMT versions.
PATENT_COLUMNS = [
    "chembl",
    "surechembl",
    "patent_id",
    "date",
    "ipc",
    "assignee",
    "family_id",
    "sections",
]


def _join_unique(values, order=None) -> str:
    """Join the distinct ``"; "``-separated items of ``values``, optionally in a given order."""
    items = {item for value in values if value for item in str(value).split("; ")}
    key = (lambda x: order.get(x, len(order))) if order else None
    return "; ".join(sorted(items, key=key))


def merge_duplicate_records(patents: pd.DataFrame) -> pd.DataFrame:
    """Keep one row per chemical and patent.

    SureChEMBL sometimes holds several compound records with the same InChIKey, so one
    ChEMBL chemical can reach the same patent through more than one SureChEMBL id. Such
    rows are merged: ``surechembl`` and ``sections`` list all values, the patent columns
    are identical anyway. Only the (few) duplicated rows are grouped, to stay fast.
    """
    key = ["chembl", "patent_id"]
    duplicated = patents.duplicated(key, keep=False)
    if not duplicated.any():
        return patents

    merged = (
        patents[duplicated]
        .groupby(key, sort=False)
        .agg(
            surechembl=("surechembl", _join_unique),
            date=("date", "first"),
            ipc=("ipc", "first"),
            assignee=("assignee", "first"),
            family_id=("family_id", "first"),
            sections=("sections", lambda v: _join_unique(v, SECTIONS)),
        )
        .reset_index()
    )
    return pd.concat([patents[~duplicated], merged[PATENT_COLUMNS]], ignore_index=True)


def extract_patent(
    analysis_name: str,
    patent_year: int = 2000,
    source: Optional[str] = None,
    sections: Optional[Iterable[str]] = None,
    ipc_codes: Optional[Iterable[str]] = VALID_CODES,
    bulk: Optional[SureChEMBLBulk] = None,
) -> pd.DataFrame:
    """Extract and store the patents mentioning the chemicals of an analysis.

    Reads ``<analysis_name>_chemicals.tsv`` (written by
    :func:`pemt.patent_extractor.patent_chemical_harmonizer.harmonize_chemicals`) and looks
    up all its SureChEMBL compounds in one query. Results go to
    ``<analysis_name>_patent_data.tsv``, with the settings and the compounds already queried
    in ``<analysis_name>_patent_data.json`` so that re-runs with the same settings only
    query new compounds.

    :param analysis_name: Name of the analysis.
    :param patent_year: Keep patents published in or after this year.
    :param source: SureChEMBL bulk data release, URL or local folder (see :func:`pemt.surechembl.resolve_source`).
    :param sections: Only count a chemical found in these patent sections, e.g. ``["claims"]``
        (see :data:`pemt.surechembl.SECTIONS`). ``None`` counts any section.
    :param ipc_codes: Keep patents with an IPC code starting with one of these; ``None`` keeps all.
    :param bulk: An open :class:`pemt.surechembl.SureChEMBLBulk` to reuse instead of ``source``.
    :returns: One row per chemical and patent with the columns in :data:`PATENT_COLUMNS`.
    """
    chemical_file = f"{PATENT_DIR}/{analysis_name}_chemicals.tsv"
    if not os.path.exists(chemical_file):
        logger.warning(f"{chemical_file} not found; run the chemical harmonizer first")
        return pd.DataFrame(columns=PATENT_COLUMNS)

    chemicals = pd.read_csv(
        chemical_file,
        sep="\t",
        dtype=str,
        usecols=["chembl", "schembl_id"],
    ).dropna()
    if chemicals.empty:
        return pd.DataFrame(columns=PATENT_COLUMNS)

    patent_file = f"{PATENT_DIR}/{analysis_name}_patent_data.tsv"
    state_file = f"{PATENT_DIR}/{analysis_name}_patent_data.json"
    settings = {
        "source": bulk.source if bulk is not None else resolve_source(source),
        "year": int(patent_year),
        "sections": sorted(set(sections)) if sections is not None else None,
        "ipc_codes": sorted(set(ipc_codes)) if ipc_codes is not None else None,
    }

    # Reuse earlier results only if they were made with the same settings.
    patents = pd.DataFrame(columns=PATENT_COLUMNS)
    done = set()
    if os.path.exists(patent_file) and os.path.exists(state_file):
        with open(state_file) as f:
            state = json.load(f)
        if state.get("settings") == settings:
            patents = pd.read_csv(patent_file, sep="\t", dtype=str).reindex(columns=PATENT_COLUMNS)
            done = set(state.get("compounds_done", []))
            # Files from before duplicate records were merged are fixed up in place.
            merged = merge_duplicate_records(patents)
            if len(merged) != len(patents):
                patents = merged.sort_values(["chembl", "date", "patent_id"], ignore_index=True)
                patents.to_csv(patent_file, sep="\t", index=False)
        else:
            logger.info("Patent settings changed since the last run; querying all compounds again")

    todo = sorted(set(chemicals["schembl_id"]) - done)
    logger.info(f"{chemicals['schembl_id'].nunique()} SureChEMBL compounds, {len(todo)} to query")

    if todo:
        own_connection = bulk is None
        bulk = bulk or SureChEMBLBulk(source)
        try:
            hits = bulk.patents_for_compounds(
                todo, ipc_prefixes=settings["ipc_codes"], min_year=patent_year, sections=sections
            )
        finally:
            if own_connection:
                bulk.close()

        hits = hits.assign(
            surechembl="SCHEMBL" + hits["compound_id"].astype(str),
            patent_id=hits["patent_number"],
            date=pd.to_datetime(hits["publication_date"]).dt.strftime("%Y-%m-%d"),
        )
        new = chemicals.rename(columns={"schembl_id": "surechembl"}).merge(hits, on="surechembl")
        new = new[PATENT_COLUMNS].fillna("").astype(str)

        patents = pd.concat([patents, new], ignore_index=True) if not patents.empty else new
        patents = merge_duplicate_records(patents.drop_duplicates())
        patents = patents.sort_values(["chembl", "date", "patent_id"], ignore_index=True)
        patents.to_csv(patent_file, sep="\t", index=False)
        with open(state_file, "w") as f:
            json.dump({"settings": settings, "compounds_done": sorted(done | set(todo))}, f, indent=2)

    current = patents[patents["chembl"].isin(set(chemicals["chembl"]))]
    logger.info(
        f"{current['patent_id'].nunique()} patents for {current['chembl'].nunique()} of "
        f"{chemicals['chembl'].nunique()} chemicals"
    )
    return current.reset_index(drop=True)
