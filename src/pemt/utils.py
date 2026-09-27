# -*- coding: utf-8 -*-

import logging
from typing import Dict, Iterable, Optional

import pandas as pd

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)

"""Protein mapper functions"""


def get_hgnc_id() -> Dict[str, str]:
    """Mapping dictionary for HGNC symbol to HGNC identifiers"""
    protein_mapping = pd.read_csv(
        f"https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_status&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit",
        sep="\t",
        index_col="Approved symbol",
    ).to_dict()["HGNC ID"]

    return protein_mapping


def hgnc_to_chembl(
    chemical_mapper: Dict[str, str], uniprot_mapper: Dict[str, str], hgnc_symbol: str
) -> Optional[str]:
    """Mapping HGNC symbol to ChEMBL identifiers.

    :param chemical_mapper: A dictionary mapping the UNIPROT identifiers to ChEMBL
    :param uniprot_mapper: A dictionary mapping the HGNC identifiers to UNIPROT
    :param hgnc_symbol: A HGNC symbol
    """
    uniprot_id = uniprot_mapper.get(hgnc_symbol)
    return uniprot_to_chembl(chemical_mapper=chemical_mapper, uniprot_id=uniprot_id)


def uniprot_to_chembl(chemical_mapper: dict, uniprot_id: str) -> Optional[str]:
    """Mapping UniProt identifiers to ChEMBL identifiers.

    :param chemical_mapper: A dictionary mapping the UNIPROT identifiers to ChEMBL
    :param uniprot_id: UNIPROT identifier of a protein
    """
    target_chembl = chemical_mapper.get(uniprot_id)
    return target_chembl


"""Chemical mapper functions"""


def get_chembl_structures(chembl_ids: Iterable[str], chunk_size: int = 50) -> pd.DataFrame:
    """Get the preferred name and standard InChIKey of ChEMBL compounds.

    :param chembl_ids: ChEMBL compound identifiers (e.g. ``CHEMBL941``)
    :param chunk_size: Number of compounds requested per ChEMBL web service call
    :returns: Columns ``chembl``, ``name`` and ``inchi_key``. Compounds without a
        structure (e.g. biologics) have no InChIKey; the name falls back to the id.
    """
    from chembl_webresource_client.new_client import new_client

    ids = sorted({i for i in chembl_ids if isinstance(i, str) and i})
    rows = []
    for start in range(0, len(ids), chunk_size):
        chunk = ids[start : start + chunk_size]
        found = new_client.molecule.filter(molecule_chembl_id__in=chunk).only(
            ["molecule_chembl_id", "pref_name", "molecule_structures"]
        )
        for molecule in found:
            structures = molecule.get("molecule_structures") or {}
            rows.append(
                {
                    "chembl": molecule["molecule_chembl_id"],
                    "name": molecule.get("pref_name") or molecule["molecule_chembl_id"],
                    "inchi_key": structures.get("standard_inchi_key"),
                }
            )

    df = pd.DataFrame(rows, columns=["chembl", "name", "inchi_key"])
    missing = set(ids) - set(df["chembl"])
    if missing:
        logger.warning(f"{len(missing)} ChEMBL ids not found in ChEMBL, e.g. {sorted(missing)[:5]}")
        df = pd.concat(
            [df, pd.DataFrame({"chembl": sorted(missing), "name": sorted(missing), "inchi_key": None})],
            ignore_index=True,
        )
    return df.sort_values("chembl", ignore_index=True)
