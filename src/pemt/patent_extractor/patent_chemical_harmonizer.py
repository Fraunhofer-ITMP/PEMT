# -*- coding: utf-8 -*-

"""Script for harmonizing the ChEMBL chemicals with patent chemicals."""

import json
import logging
import os
from collections import defaultdict

import pandas as pd
from pubchempy import get_synonyms
from tqdm import tqdm

from pemt.constants import MAPPER_DIR, PATENT_DIR
from pemt.utils import get_chemical_names

logger = logging.getLogger(__name__)


def get_surechembl_id(
    chemical_id: str, chemical_name: str, chemical_mapper: dict
) -> str:
    """Method to get SureChEMBL identifier of chemicals.

    :param chemical_id: ChEMBL identifier of the chemical
    :param chemical_name: Name of the chemical
    :param chemical_mapper: Dictionary of identifier mapping between ChEMBL and SureChEMBL
    """
    surechembl_id = chemical_mapper.get(chemical_id)

    if surechembl_id:
        return surechembl_id

    try:
        synm_dict = get_synonyms(chemical_name, namespace="name")[0]
    except IndexError:
        return None

    try:
        surechembl_id = [
            synonym for synonym in synm_dict["Synonym"] if synonym.startswith("SCHEMBL")
        ][0]
    except IndexError:
        return None

    return surechembl_id


def harmonize_chemicals(analysis_name: str, from_genes: bool = True) -> None:
    """Method that allows mapping from ChEMBL to SureChEMBL identifiers.

    :param analysis_name: The name of the analysis you want to run. This name would be used to save the resultant file.
    :param from_genes: Boolean indicating where the process needs to get chemicals based on genes or not.
    """

    # Load cached data if it exists
    if os.path.exists(f"{PATENT_DIR}/{analysis_name}_chemicals.tsv"):
        chemical_df = pd.read_csv(
            f"{PATENT_DIR}/{analysis_name}_chemicals.tsv",
            sep="\t",
            dtype=str,
        )

        cache_dict = pd.read_csv(
            f"{PATENT_DIR}/{analysis_name}_chemicals.tsv",
            sep="\t",
            dtype=str,
            index_col="chembl",
        ).to_dict()["schembl_id"]

    else:
        chemical_df = pd.DataFrame(columns=["chembl", "schembl_id", "name"])
        cache_dict = {}

    # Load chembl - schembl mapper
    with open(f"{MAPPER_DIR}/chemical_mapper.json") as f:
        chemical_mapper = json.load(f)

    # Add caching system
    cache = 0
    genes_skipped = 0

    if os.path.exists(f"{MAPPER_DIR}/{analysis_name}_chemical_names.json"):
        chemical_names = json.load(
            open(f"{MAPPER_DIR}/{analysis_name}_chemical_names.json")
        )
    else:
        chemical_names = defaultdict(str)

    if from_genes:
        if os.path.exists(f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json"):
            gene_chemical_dict = json.load(
                open(f"{MAPPER_DIR}/{analysis_name}_gene_to_chemicals.json")
            )
        else:
            raise FileNotFoundError(
                f"Please ensure that you run the experimental data extractor file first."
            )

        for genes in tqdm(
            gene_chemical_dict, desc="Harmonizing chemicals for patent retrieval"
        ):
            # Extract chemical-target data from ChEMBL
            chemical_list = gene_chemical_dict[genes]

            if len(chemical_list) < 1:
                genes_skipped += 1
                continue

            for chembl_id in chemical_list:
                # Get name for chemical and store in dict
                if chembl_id not in chemical_names:
                    chemical_names[chembl_id] = get_chemical_names(chembl_id)

                # Store chembl to surechembl mapping in separate file
                if chembl_id not in cache_dict:
                    surechembl_id = get_surechembl_id(
                        chemical_id=chembl_id,
                        chemical_name=chemical_names[chembl_id],
                        chemical_mapper=chemical_mapper,
                    )

                    if not surechembl_id:
                        continue

                    cache += 1  # add new data only
                    cache_dict[chembl_id] = surechembl_id

                    tmp_df = pd.DataFrame(
                        {
                            "chembl": chembl_id,
                            "schembl_id": surechembl_id,
                            "name": chemical_names[chembl_id],
                        },
                        index=[0],
                    )

                    chemical_df = pd.concat([chemical_df, tmp_df], ignore_index=True)

                    if cache == 10:
                        chemical_df.to_csv(
                            f"{PATENT_DIR}/{analysis_name}_chemicals.tsv",
                            sep="\t",
                            index=False,
                        )
                        cache = 0

                        # Save chemical mapping dict for re-use
                        with open(
                            f"{MAPPER_DIR}/{analysis_name}_chemical_names.json", "w"
                        ) as f:
                            json.dump(chemical_names, f, ensure_ascii=False, indent=2)
    else:

        for _, row in tqdm(
            chemical_df.iterrows(),
            total=chemical_df.shape[0],
            desc="Harmonzing chemicals for patent retrival",
        ):
            chembl_id = row["chembl"]
            # Get name for chemical and store in dict
            if chembl_id not in chemical_names:
                chemical_names[chembl_id] = get_chemical_names(chembl_id)

            # Store chembl to surechembl mapping in separate file
            if chembl_id not in cache_dict:
                surechembl_id = get_surechembl_id(
                    chemical_id=chembl_id,
                    chemical_name=chemical_names[chembl_id],
                    chemical_mapper=chemical_mapper,
                )

                if not surechembl_id:
                    continue

                cache += 1  # add new data only
                cache_dict[chembl_id] = surechembl_id

                tmp_df = pd.DataFrame(
                    {
                        "chembl": chembl_id,
                        "schembl_id": surechembl_id,
                        "name": chemical_names[chembl_id],
                    },
                    index=[0],
                )

                chemical_df = pd.concat([chemical_df, tmp_df], ignore_index=True)

                if cache == 10:
                    chemical_df.to_csv(
                        f"{PATENT_DIR}/{analysis_name}_chemicals.tsv",
                        sep="\t",
                        index=False,
                    )
                    cache = 0

                    # Save chemical mapping dict for re-use
                    with open(
                        f"{MAPPER_DIR}/{analysis_name}_chemical_names.json", "w"
                    ) as f:
                        json.dump(chemical_names, f, ensure_ascii=False, indent=2)

    chemical_df.dropna(subset=["schembl_id"], inplace=True)

    chemical_df.to_csv(
        f"{PATENT_DIR}/{analysis_name}_chemicals.tsv", sep="\t", index=False
    )

    with open(f"{MAPPER_DIR}/{analysis_name}_chemical_names.json", "w") as f:
        json.dump(chemical_names, f, ensure_ascii=False, indent=2)
