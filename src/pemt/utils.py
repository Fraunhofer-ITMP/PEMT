# -*- coding: utf-8 -*-

import logging
from typing import Dict, Optional
from urllib.error import URLError

import pandas as pd
from pubchempy import get_compounds

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)
pubchempy_logger = logging.getLogger('pubchempy')
pubchempy_logger.setLevel(logging.WARNING)

"""Protein mapper functions"""


def get_hgnc_id(
) -> Dict[str, str]:
    """Mapping dictionary for HGNC symbol to HGNC identifiers"""
    protein_mapping = pd.read_csv(
        f'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_status&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit',
        sep='\t',
        index_col='Approved symbol'
    ).to_dict()['HGNC ID']

    return protein_mapping


def hgnc_to_chembl(
    chemical_mapper: Dict[str, str],
    uniprot_mapper: Dict[str, str],
    hgnc_symbol: str
) -> Optional[str]:
    """Mapping HGNC symbol to ChEMBL identifiers.

    :param chemical_mapper: A dictionary mapping the UNIPROT identifiers to ChEMBL
    :param uniprot_mapper: A dictionary mapping the HGNC identifiers to UNIPROT
    :param hgnc_symbol: A HGNC symbol
    """
    uniprot_id = uniprot_mapper.get(hgnc_symbol)
    return uniprot_to_chembl(chemical_mapper=chemical_mapper, uniprot_id=uniprot_id)


def uniprot_to_chembl(
    chemical_mapper: dict,
    uniprot_id: str
) -> Optional[str]:
    """Mapping UniProt identifiers to ChEMBL identifiers.

    :param chemical_mapper: A dictionary mapping the UNIPROT identifiers to ChEMBL
    :param uniprot_id: UNIPROT identifier of a protein
    """
    target_chembl = chemical_mapper.get(uniprot_id)
    return target_chembl


"""Chemical mapper functions"""


def get_chemical_names(
    chembl_id: str
) -> str:
    """Method to get chemical name from ChEMBL id.

    :param chembl_id: ChEMBL identifier of a compound
    """
    try:
        chemical_name = get_compounds(chembl_id, 'name')[0].synonyms[0]
    except (IndexError, URLError):
        chemical_name = chembl_id
    return chemical_name
