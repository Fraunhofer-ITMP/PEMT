# -*- coding: utf-8 -*-

"""Script for showing a working example of the PEMT workflow."""

import logging
import os
import xml.etree.ElementTree as ET

import pandas as pd
from tqdm import tqdm

from pemt.chemical_extractor.experimental_data_extraction import extract_chemicals
from pemt.constants import PATENT_DIR
from pemt.patent_extractor.patent_chemical_harmonizer import harmonize_chemicals
from pemt.patent_extractor.patent_enrichment import extract_patent

logger = logging.getLogger(__name__)

chromedriver_path = "C:/Users/Yojana.Gadiya/Downloads/chromedriver"


def get_top_disease():
    """Get top 5 rare diseases based on the epidemiology found in Orphanet"""

    # Data available at http://www.orphadata.org/cgi-bin/epidemio.html
    tree = ET.parse("orphanet_epidem_data.xml")
    root = tree.getroot()

    orphanet_data = []

    for disease_data in tqdm(root.find("DisorderList"), desc="Reading XML file"):
        df_data = {}
        for data in disease_data:
            if data.tag == "OrphaCode":
                df_data["orphanet id"] = data.text
            elif data.tag == "Name":
                df_data["name"] = data.text
            elif data.tag == "PrevalenceList":
                df_data["prevalence count"] = int(data.attrib.get("count"))

        orphanet_data.append(df_data)

    orphanet_df = pd.DataFrame(orphanet_data)
    orphanet_df.sort_values(by="prevalence count", ascending=False, inplace=True)
    orphanet_df.to_csv("epidemologic_distribution_orphanet.tsv", sep="\t", index=False)


def get_top_gene():
    """Get genes linked to top 5 diseases (based on epidemiological data)"""

    top_diseases = pd.read_csv(
        "epidemologic_distribution_orphanet.tsv", sep="\t", dtype=str
    )

    # Data available at http://www.orphadata.org/cgi-bin/index.php
    tree = ET.parse("orphanet_disgene_data.xml")
    root = tree.getroot()

    orphanet_genes = []

    for disease_data in tqdm(root.find("DisorderList"), desc="Reading XML file"):
        df_data = {}
        for data in disease_data:
            if data.tag == "OrphaCode":
                df_data["orphanet id"] = data.text
            elif data.tag == "Name":
                df_data["name"] = data.text
            elif data.tag == "DisorderGeneAssociationList":
                for rel_data in data:
                    for rel in rel_data:
                        if rel.tag == "Gene":
                            df_data["Symbol"] = rel.find("Symbol").text

        orphanet_genes.append(df_data)

    orphanet_gene_df = pd.DataFrame(orphanet_genes)

    # Merge the data to add prelevance count to each row
    gene_prelevance_df = orphanet_gene_df.merge(top_diseases, on="orphanet id")
    gene_prelevance_df["prevalence count"] = gene_prelevance_df[
        "prevalence count"
    ].astype(int)

    # Sort column by prelevance
    gene_prelevance_df.sort_values(by="prevalence count", ascending=False, inplace=True)

    gene_prelevance_df.to_csv("orphanet_gene_by_prevelance.tsv", sep="\t", index=False)


def run_from_gene_pipeline(name: str, genes: list, os: str):
    extract_chemicals(
        analysis_name=name,
        gene_list=genes,
        is_uniprot=False,
    )

    harmonize_chemicals(analysis_name=name, from_genes=True)

    patent_df = extract_patent(
        analysis_name=name,
        chrome_driver_path=chromedriver_path,
        os_system=os,
        patent_year=2000,
    )

    if patent_df.empty:
        logger.info(f"No patents found!")
        return None

    # Since the original patent data has chemical with no patents, we remove those entries from the data
    patent_df = patent_df[~patent_df["patent_id"].isna()]
    patent_df.to_csv(
        f"{PATENT_DIR}/cleaned_{name}_patent_data.tsv", sep="\t", index=False
    )

    logger.info(f"Done with retrival of patents")
    logger.info(f"Data file can be found under {PATENT_DIR}")
    return


def run_from_patent_pipeline(name: str, chemical_data: str, os: str):
    df = pd.read_csv(chemical_data, sep="\t", dtype=str)

    df.to_csv(f"{PATENT_DIR}/{name}_chemicals.tsv", sep="\t", index=False)

    harmonize_chemicals(analysis_name=name, from_genes=False)

    patent_df = extract_patent(
        analysis_name=name,
        chrome_driver_path=chromedriver_path,
        os_system=os,
        patent_year=2000,
    )

    if patent_df.empty:
        logger.info(f"No patents found!")
        return None

    # Since the original patent data has chemical with no patents, we remove those entries from the data
    patent_df = patent_df[~patent_df["patent_id"].isna()]
    patent_df.to_csv(
        f"{PATENT_DIR}/cleaned_{name}_patent_data.tsv", sep="\t", index=False
    )

    logger.info(f"Done with retrival of patents")
    logger.info(f"Data file can be found under {PATENT_DIR}")
    return


def main():
    """Main function to demonstrate usecase of PEMT."""

    if not os.path.exists("epidemologic_distribution_orphanet.tsv"):
        get_top_disease()

    if not os.path.exists("orphanet_gene_by_prevelance.tsv"):
        get_top_gene()

    df = pd.read_csv(
        "orphanet_gene_by_prevelance.tsv",
        sep="\t",
        usecols=["Symbol", "prevalence count"],
        dtype=str,
    )
    df = df[df["prevalence count"].astype(int) > 9]
    gene_list = set(df["Symbol"].to_list())

    run_from_gene_pipeline(name="rare disease", genes=list(gene_list), os="windows")


if __name__ == "__main__":
    main()
