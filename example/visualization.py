# -*- coding: utf-8 -*-

"""Script for analysing the data from PEMT."""

import json
import os
import xml.etree.ElementTree as ET
from collections import defaultdict

import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from pemt.constants import DATA_DIR, MAPPER_DIR, PATENT_DIR

PLOT_DIR = f"{DATA_DIR}/plots"
os.makedirs(PLOT_DIR, exist_ok=True)


def gene_chemical_analysis():
    """Method analysing the gene-chemical data."""
    gene_dict = json.load(open(f"{MAPPER_DIR}/rare disease_gene_to_chemicals.json"))
    print(f"Number of genes - {len(gene_dict)}")

    chemicals = [
        chemical for _, chemcials in gene_dict.items() for chemical in chemcials
    ]
    chemicals = set(chemicals)
    print(f"Number of chemicals - {len(chemicals)}")

    new_data = []
    for gene, chemicals in tqdm(gene_dict.items()):
        new_data.append({"gene": gene, "chemicals": len(chemicals)})

    df = pd.DataFrame(new_data)
    df.sort_values(by="chemicals", ascending=False, inplace=True)
    df = df[df["chemicals"] > 0]

    # grid lines
    plt.figure(figsize=(12, 5))
    sns.set_style("darkgrid")

    bx = sns.barplot(x="gene", y="chemicals", data=df, color="c", dodge=False)
    for item in bx.get_xticklabels():
        item.set_rotation(45)

    bx.set_xlabel("Gene symbol", fontsize=14, fontweight="bold")
    bx.set_ylabel("Number of chemicals", fontsize=14, fontweight="bold")

    # Set bar labels
    for container in bx.containers:
        bx.bar_label(container)

    plt.tight_layout()
    plt.savefig(f"{PLOT_DIR}/gene_chemical_analysis.jpg")
    plt.show()


def chemical_patent_counter():
    """Method counting the patents linked to chemcials."""
    chemical_patent_df = pd.read_csv(
        f"{PATENT_DIR}/rare disease_patent_data.tsv", sep="\t"
    )
    chemical_patent_df = chemical_patent_df[chemical_patent_df["patent_id"].notna()]

    # Create dataframe of chemical with patent counts
    chemical_counter = chemical_patent_df["chembl"].value_counts().reset_index()
    chemical_counter.columns = ["chembl", "patent counts"]

    patents = set(chemical_patent_df["patent_id"].unique())

    print(f"Number of chemicals - {chemical_counter.shape[0]}")
    print(f"Number of patents - {len(patents)}")
    print(chemical_counter["patent counts"].value_counts())


def patent_analysis():
    """Method counting the patents linked to chemicals."""
    chemical_patent_df = pd.read_csv(
        f"{PATENT_DIR}/rare disease_patent_data.tsv", sep="\t"
    )
    chemical_patent_df = chemical_patent_df[chemical_patent_df["patent_id"].notna()]

    patents = set(chemical_patent_df["patent_id"].unique())

    countries = defaultdict(int)
    patent_position = defaultdict(int)

    for patent in tqdm(patents):
        country = patent.split("-")[0]
        position = patent.split("-")[-1]
        countries[country] += 1
        patent_position[position] += 1

    print(countries)
    print(patent_position)

    patent_year = defaultdict(set)
    ipc_codes = defaultdict(set)
    assignees = set()

    for row in tqdm(chemical_patent_df.values):
        (chemical_id, surechembl_id, patent_id, date, ipc, assignee) = row

        year = date.split("-")[0]
        patent_year[year].add(patent_id)
        code = ipc.split()[0]
        ipc_codes[code].add(patent_id)
        assignees.add(assignee)

    patent_year = {year: len(patents) for year, patents in patent_year.items()}
    ipc_codes = {ipc: len(patents) for ipc, patents in ipc_codes.items()}

    print(ipc_codes)

    year_df = pd.DataFrame(patent_year, index=[0]).transpose().reset_index()
    year_df.columns = ["year", "count"]
    year_df.sort_values(by="year", inplace=True)

    # grid lines
    plt.figure(figsize=(12, 5))
    sns.set_style("darkgrid")

    bx = sns.barplot(x="year", y="count", data=year_df, color="c", dodge=False)
    for item in bx.get_xticklabels():
        item.set_rotation(45)

    bx.set_xlabel("Year", fontsize=14, fontweight="bold")
    bx.set_ylabel("Number of patents", fontsize=14, fontweight="bold")

    # Set bar labels
    for container in bx.containers:
        bx.bar_label(container)

    plt.tight_layout()
    plt.savefig(f"{PLOT_DIR}/patent_year_analysis.jpg")
    plt.show()


def orphanet_data_analysis():
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
    diseases = set(orphanet_gene_df["orphanet id"].tolist())
    genes = set(orphanet_gene_df["Symbol"].tolist())

    print(f"Orphanet has {len(diseases)} diseases and {len(genes)} genes.")


def gene_data_generator():
    """Method to annotate patents with gene."""
    chemical_patent_df = pd.read_csv(
        f"{PATENT_DIR}/rare disease_patent_data.tsv", sep="\t"
    )
    chemical_patent_df = chemical_patent_df[chemical_patent_df["patent_id"].notna()]

    # Genes chemical data
    gene_dict = json.load(open(f"{MAPPER_DIR}/rare disease_gene_to_chemicals.json"))

    chemical_dict = defaultdict(set)

    for gene, chemical_ids in tqdm(
        gene_dict.items(), desc="Creating the chemical-gene dictionary"
    ):
        for chemical in chemical_ids:
            chemical_dict[chemical].add(gene)

    chemical_patent_df["genes"] = chemical_patent_df["chembl"].map(
        lambda x: ", ".join(chemical_dict.get(x, []))
    )

    chemical_patent_df.to_csv(
        f"{PATENT_DIR}/gene_enumerated_patent_data.tsv", sep="\t", index=False
    )


def assignee_analysis():
    """Method to normalize the assignee to their assignee type."""
    df = pd.read_csv(
        f"{PATENT_DIR}/gene_enumerated_patent_data.tsv",
        sep="\t",
        usecols=["patent_id", "assignee"],
    )

    assignee_dict = {
        "COHEN FREDERICK": "Private",
        "GENENTECH INC": "Industry",
        "ARVINAS OPERATIONS INC": "Industry",
        "CULLGEN SHANGHAI INC": "Industry",
        "PORTER DALE": "Private",
        "ARVINAS INC": "Industry",
        "UNIV JOHNS HOPKINS": "Academia",
        "WARNER LAMBERT CO": "Industry",
        "THERAVANCE INC": "Industry",
        "ADVANCED MEDICINE INC": "Industry",
        "EISAI R&D MAN CO LTD": "Industry",
        "THERAVANCE BIOPHARMA R&D IP, LLC": "Industry",
        "ESAI R & D MAN CO LTD": "Industry",
        "GRIFFIN JOHN H": "Private",
        "US HEALTH": "Industry",
        "IMCLONE SYSTEMS INC": "Industry",
        "FUJIFILM CELLULAR DYNAMICS INC": "Industry",
        "ALMAC DISCOVERY LTD": "Industry",
        "SCHERING CORP": "Industry",
        "WONG ALBERT J": "Private",
        "INST QUIMIC DE SARRIA CETS FUNDACIO PRIVADA": "Academia",
        "HILLAN KENNETH J": "Private",
        "SESHAGIRI SOMASEKAR": "Private",
        "AUCKLAND UNISERVICES LTD": "Academia",
        "AXXIMA PHARMACEUTICALS AG": "Industry",
        "BHIDE RAJEEV S": "Private",
        "SQUIBB BRISTOL MYERS CO": "Industry",
        "BRISTOL-MYERS SQUIBB COMPANY": "Industry",
        "HOFFMANN LA ROCHE": "Industry",
        "JANSSEN PHARMACEUTICA NV": "Industry",
        "ALECTOS THERAPEUTICS INC": "Industry",
        "ACADEMIA SINICA": "Academia",
        "BAYER SCHERING PHARMA AG": "Industry",
        "HARVARD COLLEGE": "Academia",
        "VASCULAR BIOGENICS LTD": "Industry",
        "YUAN JUNYING": "Private",
        "PRENDERGAST GEORGE C": "Private",
        "INST MEDIZINTECHNOLOGIE MAGDEB": "Academia",
        "ELAN PHARM INC": "Industry",
        "SEVERSON BILL": "Private",
        "UNIV UTAH RES FOUND": "Academia",
        "UNIV ROCHESTER": "Academia",
        "FABIUS BIOTECHNOLOGY": "Industry",
        "KANONCURE INC": "Industry",
        "KADMON CORP LLC": "Industry",
        "ABBOTT LABORATORIES": "Industry",
        "ABBVIE INC.": "Industry",
        "NOVARTIS AG": "Industry",
        "CHAN HO MAN": "Private",
        "BROGDON JENNIFER": "Private",
        "CERO THERAPEUTICS INC": "Industry",
        "ASHWELL MARK A": "Private",
        "ARQULE INC": "Industry",
        "UNIV PENNSYLVANIA": "Academia",
        "ALI SYED M": "Private",
        "CHIRON CORP": "Industry",
        "ASTEX THERAPEUTICS LTD": "Industry",
        "LEXICON PHARMACEUTICALS INC": "Industry",
        "MANSFIELD ROBERT K": "Private",
        "BUSCHMANN NICOLE": "Private",
        "BLANCHARD STÉPHANIE": "Private",
        "UNIV NEW YORK": "Academia",
        "S BIO PTE LTD": "Industry",
        "DESIGN THERAPEUTICS INC": "Industry",
        "UNIV QUEENSLAND": "Academia",
        "IND RES LTD": "Industry",
        "EINSTEIN COLL MED": "Academia",
        "ATHYRIUM OPPORTUNITIES III CO-INVEST 1 LP": "Industry",
        "FURNEAUX RICHARD H": "Private",
        "SINCLAIR DAVID A": "Private",
        "VIFOR INT AG": "Industry",
        "UNIV COLUMBIA": "Academia",
        "UNIV BOSTON": "Academia",
        "UNIV SOUTH FLORIDA": "Academia",
    }

    df["assignee_type"] = df["assignee"].map(assignee_dict)

    df.drop_duplicates(inplace=True)

    print(df["assignee_type"].value_counts())


def gene_patent_analysis():
    """Method for generating the gene-patent data plot."""
    df = pd.read_csv(
        f"{PATENT_DIR}/gene_enumerated_patent_data.tsv",
        sep="\t",
        usecols=["patent_id", "genes"],
    )

    df.drop_duplicates(inplace=True)

    # Grouping patents by genes
    df = df.groupby(by=["genes"]).count().reset_index()
    df.sort_values(by="patent_id", inplace=True, ascending=False)

    # grid lines
    plt.figure(figsize=(12, 5))
    sns.set_style("darkgrid")

    bx = sns.barplot(x="genes", y="patent_id", data=df, color="c", dodge=False)
    for item in bx.get_xticklabels():
        item.set_rotation(45)

    bx.set_xlabel("Gene symbol", fontsize=14, fontweight="bold")
    bx.set_ylabel("Number of patents", fontsize=14, fontweight="bold")

    # Set bar labels
    for container in bx.containers:
        bx.bar_label(container)

    plt.tight_layout()
    plt.savefig(f"{PLOT_DIR}/gene_patent_analysis.jpg")
    plt.show()


if __name__ == "__main__":
    gene_patent_analysis()
