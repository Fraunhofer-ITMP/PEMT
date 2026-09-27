# -*- coding: utf-8 -*-

"""Script for extracting experimental bioassay information from ChEMBL."""

import datetime
import json
import logging
import os
from collections import defaultdict
from typing import List

import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

from pemt.constants import MAPPER_DIR
from pemt.utils import get_chembl_release, get_single_protein_targets, load_hgnc, symbols_to_uniprot

logger = logging.getLogger(__name__)

# Change logging level for packages
chembl_logger = logging.getLogger("chembl_webresource_client")
chembl_logger.setLevel(logging.WARNING)

activity = new_client.activity

tqdm.pandas()

os.makedirs(MAPPER_DIR, exist_ok=True)


def get_chemical_overview(file_path: str) -> None:
    """Method to report incomplete information in the chemical enrichment.

    :param file_path: Path of the JSON file storing the gene and chemical information.
    """
    gene_chemical_dict = json.load(open(file_path))

    counter_dict = defaultdict(int)

    for gene, counter in tqdm(gene_chemical_dict.items()):
        counter_dict[gene] += len(counter)

    counter_dict = dict(
        sorted(counter_dict.items(), key=lambda item: item[1], reverse=True)
    )

    df = pd.DataFrame(counter_dict, index=[0]).transpose()
    without_chemicals = int((df[0] == 0).sum()) if not df.empty else 0
    if without_chemicals:
        logger.warning(
            f"{without_chemicals} genes found with no relevant chemical bioassay information."
        )


def target_to_chemical(
    protein: str,
    protein_mapping: dict = None,
    is_uniprot: bool = False,
) -> List[str]:
    """Method to retrieve bioactive chemicals, from proteins, based on biochemical/ functional bioassays.
    A chemical is considered active if it has a pChEMBL >= 6.

    :param protein: The protein name or identifier
    :param protein_mapping: A dictionary mapping the HGNC symbols to UNIPROT identifiers (a list, or a
    string with several ids separated by "," or "|"), e.g. from :func:`pemt.utils.symbols_to_uniprot`.
    By default, the value is set to None.
    :param is_uniprot: Boolean indicating whether the protein is an HGNC symbol or UNIPROT identifier.
    If using UniProt ids for protein, set the value to "True" and the protein_mapping parameter can be omitted.
    If using HGNC symbols, then the protein mapping dictionary needs to be provided.
    :returns: ChEMBL ids of the active chemicals, without duplicates.
    """
    if is_uniprot:
        uniprot_ids = [protein]
    else:
        if protein_mapping is None:
            raise ValueError(
                "HGNC symbol given without passing the HGNC to UNIPROT mapping file. "
                "Either pass the mapping file to protein_mapping or set the parameter is_uniprot=True"
            )
        uniprot_ids = protein_mapping.get(protein) or []
        if isinstance(uniprot_ids, str):  # e.g. "O43687, Q9P0M2" or "O43687|Q9P0M2"
            uniprot_ids = [u.strip() for u in uniprot_ids.replace("|", ",").split(",") if u.strip()]

    targets = sorted({t for uniprot_id in uniprot_ids for t in get_single_protein_targets(uniprot_id)})

    chemicals = {}
    for target_chembl in targets:
        prot_activity_data = activity.filter(
            target_chembl_id=target_chembl,
            assay_type_iregex="(B|F)",
            pchembl_value__gte=6,
        ).only(["pchembl_value", "molecule_chembl_id"])

        for i in prot_activity_data:
            pchembl_val = i["pchembl_value"]
            if pchembl_val is None or pd.isna(pchembl_val) or float(pchembl_val) < 6:
                continue
            chemicals[i["molecule_chembl_id"]] = None  # dict keeps first-seen order

    logger.debug(f"{protein}: {len(chemicals)} active chemicals")
    return list(chemicals)


def extract_chemicals(
    analysis_name: str,
    gene_list: list = None,
    gene_file_path: str = None,
    file_separator: str = "comma",
    is_uniprot: bool = False,
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

    # Loop to get and store the genes-chemical information from ChEMBL
    if os.path.exists(f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json"):
        gene_chemical_dict = json.load(
            open(f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json")
        )
    else:
        gene_chemical_dict = defaultdict()

    new_count = 0

    # Extract the gene
    if file_separator == "comma":
        _separator = ","
    elif file_separator == "semicolon":
        _separator = ";"
    else:
        assert file_separator == "tab"
        _separator = "\t"

    if gene_file_path:
        df = pd.read_csv(gene_file_path, sep=_separator)
        if is_uniprot:
            if "uniprot" not in list(df.columns):
                raise ValueError(
                    f'Please rename columns to : "uniprot" in case of uniprot id or "symbol" in case of HGNC symbols'
                )
        elif not is_uniprot:
            if "symbol" not in list(df.columns):
                raise ValueError(
                    f'Please rename columns to : "uniprot" in case of uniprot id or "symbol" in case of HGNC symbols'
                )

    if gene_file_path and _separator and is_uniprot:
        proteins = set(pd.read_csv(gene_file_path, sep=_separator)["uniprot"].tolist())
    elif gene_file_path and _separator and not is_uniprot:
        proteins = set(pd.read_csv(gene_file_path, sep=_separator)["symbol"].tolist())
    else:
        proteins = gene_list

    proteins = sorted({str(p).strip() for p in proteins if isinstance(p, str) and str(p).strip()})

    # Record which releases this run used (ChEMBL changes a few times a year).
    info_file = f"{MAPPER_DIR}/{analysis_name}_run_info.json"
    run_info = json.load(open(info_file)) if os.path.exists(info_file) else {}
    chembl_release = get_chembl_release()
    if run_info.get("chembl_release") and chembl_release and run_info["chembl_release"] != chembl_release:
        logger.warning(
            f"Cached results of '{analysis_name}' come from {run_info['chembl_release']}, ChEMBL now serves "
            f"{chembl_release}; delete {analysis_name}_gene_to_chemicals.json to redo them"
        )
    run_info.setdefault("chembl_release", chembl_release)
    run_info["input_type"] = "uniprot" if is_uniprot else "symbol"

    # Gene symbols -> UniProt accessions via the current HGNC complete set
    hgnc_mapper = None
    if not is_uniprot:
        todo = [p for p in proteins if p not in gene_chemical_dict]
        if todo:
            hgnc = load_hgnc(MAPPER_DIR)
            hgnc_mapper = symbols_to_uniprot(todo, hgnc)
            hgnc_file = os.path.join(MAPPER_DIR, "hgnc_complete_set.txt")
            run_info["hgnc_file_date"] = datetime.date.fromtimestamp(os.path.getmtime(hgnc_file)).isoformat()
            run_info["symbol_to_uniprot"] = {**run_info.get("symbol_to_uniprot", {}), **hgnc_mapper}

    # Loop to get chemicals related to target
    for identifier in tqdm(proteins, desc="Extracting chemicals for targets"):
        if identifier in gene_chemical_dict:
            continue

        new_count += 1
        # Extract chemical-target data from ChEMBL
        chemical_list = target_to_chemical(
            protein=identifier,
            protein_mapping=hgnc_mapper,
            is_uniprot=is_uniprot,
        )
        gene_chemical_dict[identifier] = chemical_list

        if new_count == 50:
            with open(f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json", "w") as f:
                json.dump(gene_chemical_dict, f, ensure_ascii=False, indent=2)
            new_count = 0

    # Save dict for re-use
    if new_count > 0 or not os.path.exists(f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json"):
        with open(f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json", "w") as f:
            json.dump(gene_chemical_dict, f, ensure_ascii=False, indent=2)

    with open(info_file, "w") as f:
        json.dump(run_info, f, ensure_ascii=False, indent=2)

    # Get genes with no chemical hits
    get_chemical_overview(f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json")

    return gene_chemical_dict
