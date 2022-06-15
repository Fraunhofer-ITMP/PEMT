# -*- coding: utf-8 -*-

"""Script for extracting experimental bioassay information from ChEMBL."""

import json
import logging
import os
from collections import defaultdict
from typing import List

import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

from pemt.constants import MAPPER_DIR
from pemt.utils import hgnc_to_chembl, uniprot_to_chembl

logger = logging.getLogger(__name__)

# Change logging level for packages
chembl_logger = logging.getLogger("chembl_webresource_client")
chembl_logger.setLevel(logging.WARNING)

activity = new_client.activity

tqdm.pandas()


def get_chemical_overview(
    file_path: str
) -> None:
    """Method to report incomplete information in the chemical enrichment.

    :param file_path: Path of the JSON file storing the gene and chemical information.
    """
    gene_chemical_dict = json.load(open(file_path))

    counter_dict = defaultdict(int)

    for gene, counter in tqdm(gene_chemical_dict.items()):
        counter_dict[gene] += len(counter)

    counter_dict = dict(sorted(counter_dict.items(), key=lambda item: item[1], reverse=True))

    df = pd.DataFrame(counter_dict, index=[0]).transpose()
    value_count_dict = df[0].value_counts().to_dict()
    logger.warning(f'{value_count_dict[0]} genes found with no relevant chemical bioassay information.')


def target_to_chemical(
    chemical_mapping: dict,
    protein: str,
    protein_mapping: dict = None,
    is_uniprot: bool = False
) -> List[dict]:
    """Method to retrieve bioactive chemicals, from proteins, based on biochemical/ functional bioassays.
    A chemical is considered active if it has a pChEMBL > 6.

    :param chemical_mapping: A dictionary mapping the UNIPROT identifiers to ChEMBL identifiers
    :param protein: The protein name or identifier
    :param protein_mapping: A dictionary mapping the HGNC symbols to UNIPROT identifiers.
    By default, the value is set to None.
    :param is_uniprot: Boolean indicating whether the protein is an HGNC symbol or UNIPROT identifier.
    If using UniProt ids for protein, set the value to "True" and the protein_mapping parameter can be omitted.
    If using HGNC symbols, then the protein mapping dictionary needs to be provided.
    """
    chemicals = []

    if not is_uniprot:
        try:
            assert protein_mapping is not None
        except AssertionError:
            raise ValueError(f'HGNC symbol given without passing the HGNC to UNIPROT mapping file. \
            Either pass the mapping file to hgnc_mapping variable or set the parameter is_uniprot=True')
        target_chembl = hgnc_to_chembl(
            uniprot_mapper=protein_mapping,
            chemical_mapper=chemical_mapping,
            hgnc_symbol=protein
        )
    else:
        target_chembl = uniprot_to_chembl(
            chemical_mapper=chemical_mapping,
            uniprot_id=protein
        )

    if not target_chembl:
        return chemicals

    prot_activity_data = activity.filter(
        target_chembl_id=target_chembl,
        assay_type_iregex='(B|F)',
    ).only(['pchembl_value', 'molecule_chembl_id'])

    if len(prot_activity_data) < 1:
        return chemicals

    logger.debug(f'Analysing {len(prot_activity_data)} chemicals')
    for i in prot_activity_data:
        pchembl_val = i['pchembl_value']

        if pd.isna(pchembl_val):
            continue

        if float(pchembl_val) < 6:
            continue

        chemicals.append(i['molecule_chembl_id'], )

    return chemicals


def extract_chemicals(
    analysis_name: str,
    gene_list: list = None,
    gene_file_path: str = None,
    file_separator: str = 'comma',
    is_uniprot: bool = False,
    chembl_version: str = '30'
):
    """Enrich genes with chemical data from CheMBL bioassays.

    :param analysis_name: The name of the analysis you want to run. This name would be used to save the resultant file
    :param gene_list: The list of gene you want to extract chemicals for.
    :param gene_file_path: The path of the gene file
    :param file_separator: The separator used within the file. This can be 'comma', 'tab', or 'semicolon'.  By default,
    the file separator is set to csv.
    :param is_uniprot: A boolean value indicating whether the given gene list or file containing uniprot ids or HGNC
    symbols. By default, the value is set to False indicating that a "symbol" column is present with the respective
    HGNC symbols. If set to True, the file with "uniprot" column is expected.
    """

    # Load chembl target mapper files
    chembl_mapper = pd.read_csv(
        f'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{chembl_version}/chembl_uniprot_mapping.txt',
        dtype=str,
        skiprows=1,
        sep='\t',
        names=['uniprot', 'chembl_id', 'name', 'type']
    )
    chembl_mapper = chembl_mapper[['uniprot', 'chembl_id']]
    chembl_mapper.set_index('uniprot', inplace=True)
    chembl_mapper = chembl_mapper.to_dict()['chembl_id']

    hgnc_mapper = pd.read_csv(
        f'{MAPPER_DIR}/hgnc_mapper.tsv',
        sep='\t',
        index_col='Approved symbol'
    ).to_dict()['UniProt ID(supplied by UniProt)']

    # Loop to get and store the genes-chemical information from ChEMBL
    if os.path.exists(f'{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json'):
        gene_chemical_dict = json.load(open(f'{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json'))
    else:
        gene_chemical_dict = defaultdict()

    new_count = 0

    # Extract the gene
    if file_separator == 'comma':
        _separator=','
    elif file_separator == 'semicolon':
        _separator = ';'
    else:
        assert file_separator == 'tab'
        _separator = '\t'


    if gene_file_path:
        df = pd.read_csv(
            gene_file_path,
            sep=_separator
        )
        if is_uniprot:
            if 'uniprot' not in list(df.columns):
                raise ValueError(
                    f'Please rename columns to : "uniprot" in case of uniprot id or "symbol" in case of HGNC symbols'
                )
        elif not is_uniprot:
            if 'symbol' not in list(df.columns):
                raise ValueError(
                    f'Please rename columns to : "uniprot" in case of uniprot id or "symbol" in case of HGNC symbols'
                )

    if gene_file_path and _separator and is_uniprot:
        proteins = set(pd.read_csv(
            gene_file_path,
            sep=_separator
        )['uniprot'].tolist())
    elif gene_file_path and _separator and not is_uniprot:
        proteins = set(pd.read_csv(
            gene_file_path,
            sep=_separator
        )['symbol'].tolist())
    else:
        proteins = gene_list

    # Loop to get chemicals related to target
    for identifier in tqdm(proteins, desc='Extracting chemicals for targets'):
        if identifier in gene_chemical_dict:
            continue

        new_count += 1
        # Extract chemical-target data from ChEMBL
        chemical_list = target_to_chemical(
            protein=identifier,
            protein_mapping=hgnc_mapper,
            chemical_mapping=chembl_mapper,
            is_uniprot=is_uniprot
        )
        gene_chemical_dict[identifier] = chemical_list

    # Save dict for re-use
    if new_count > 0:
        with open(f'{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json', 'w') as f:
            json.dump(gene_chemical_dict, f, ensure_ascii=False, indent=2)

    # Get genes with no chemical hits
    get_chemical_overview(f'{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json')

    return gene_chemical_dict
