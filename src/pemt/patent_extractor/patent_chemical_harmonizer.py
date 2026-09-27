# -*- coding: utf-8 -*-

"""Script for harmonizing the ChEMBL chemicals with patent chemicals."""

import json
import logging
import os
from typing import Optional

import pandas as pd

from pemt.constants import MAPPER_DIR, PATENT_DIR
from pemt.surechembl import SureChEMBLBulk
from pemt.utils import get_chembl_structures

logger = logging.getLogger(__name__)

os.makedirs(MAPPER_DIR, exist_ok=True)
os.makedirs(PATENT_DIR, exist_ok=True)

CHEMICAL_COLUMNS = ["chembl", "schembl_id", "name", "inchi_key"]


def _chembl_ids_for_run(analysis_name: str, from_genes: bool) -> list:
    """Collect the ChEMBL compound ids an analysis needs."""
    if from_genes:
        gene_file = f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json"
        if not os.path.exists(gene_file):
            raise FileNotFoundError(
                "Please ensure that you run the experimental data extractor file first."
            )
        with open(gene_file) as f:
            gene_chemical_dict = json.load(f)
        genes_skipped = sum(1 for chemicals in gene_chemical_dict.values() if not chemicals)
        if genes_skipped:
            logger.info(f"{genes_skipped} genes have no active chemicals and are skipped")
        return sorted({c for chemicals in gene_chemical_dict.values() for c in chemicals})

    chemical_file = f"{PATENT_DIR}/{analysis_name}_chemicals.tsv"
    df = pd.read_csv(chemical_file, sep="\t", dtype=str)
    if "chembl" not in df.columns:
        raise ValueError(f'{chemical_file} needs a "chembl" column with ChEMBL compound ids')
    return sorted(set(df["chembl"].dropna()))


def harmonize_chemicals(
    analysis_name: str,
    from_genes: bool = True,
    source: Optional[str] = None,
    bulk: Optional[SureChEMBLBulk] = None,
) -> pd.DataFrame:
    """Map ChEMBL compounds to SureChEMBL identifiers through their InChIKeys.

    Names and InChIKeys come from ChEMBL; the InChIKeys are then looked up in the
    SureChEMBL bulk data. Results are cached in ``<analysis_name>_chemicals.tsv``
    (including compounds without a SureChEMBL match), so re-runs only look up new ids.

    :param analysis_name: The name of the analysis you want to run. This name would be used to save the resultant file.
    :param from_genes: Take the chemicals from the chemical extractor output (True) or from a user provided
        ``<analysis_name>_chemicals.tsv`` with a ``chembl`` column (False).
    :param source: SureChEMBL bulk data release, URL or local folder (see :func:`pemt.surechembl.resolve_source`).
    :param bulk: An open :class:`pemt.surechembl.SureChEMBLBulk` to reuse instead of ``source``.
    :returns: One row per ChEMBL compound of this run with the columns in :data:`CHEMICAL_COLUMNS`.
    """
    chemical_file = f"{PATENT_DIR}/{analysis_name}_chemicals.tsv"
    chembl_ids = _chembl_ids_for_run(analysis_name, from_genes)

    # Reuse earlier results; a user-provided input file has no inchi_key column yet.
    cached = pd.DataFrame(columns=CHEMICAL_COLUMNS)
    previous_columns = []
    if os.path.exists(chemical_file):
        previous = pd.read_csv(chemical_file, sep="\t", dtype=str)
        previous_columns = list(previous.columns)
        if "inchi_key" in previous.columns:
            cached = previous.reindex(columns=CHEMICAL_COLUMNS)

    todo = sorted(set(chembl_ids) - set(cached["chembl"]))
    logger.info(f"{len(chembl_ids)} chemicals, {len(todo)} not yet mapped to SureChEMBL")

    if todo:
        structures = get_chembl_structures(todo)
        keys = structures["inchi_key"].dropna().tolist()

        if keys:
            own_connection = bulk is None
            bulk = bulk or SureChEMBLBulk(source)
            try:
                mapping = bulk.compounds_for_inchikeys(keys)
            finally:
                if own_connection:
                    bulk.close()
            mapping["schembl_id"] = "SCHEMBL" + mapping.pop("surechembl_id").astype(str)
        else:
            mapping = pd.DataFrame(columns=["inchi_key", "schembl_id"])

        new = structures.merge(mapping, on="inchi_key", how="left")[CHEMICAL_COLUMNS]
        logger.info(
            f"Mapped {new.dropna(subset=['schembl_id'])['chembl'].nunique()} of {len(todo)} "
            f"new chemicals to SureChEMBL"
        )
        cached = pd.concat([cached, new], ignore_index=True) if not cached.empty else new
        cached = cached.sort_values(["chembl", "schembl_id"], ignore_index=True)

    # Always write the file (even when empty) so the patent step finds it.
    if todo or not os.path.exists(chemical_file) or "inchi_key" not in previous_columns:
        cached.to_csv(chemical_file, sep="\t", index=False)

    return cached[cached["chembl"].isin(chembl_ids)].reset_index(drop=True)
