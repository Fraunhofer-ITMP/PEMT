# -*- coding: utf-8 -*-

"""Command line interface."""

import json
import logging
from collections import defaultdict

import click
import pandas as pd

from pemt.chemical_extractor.experimental_data_extraction import extract_chemicals
from pemt.constants import MAPPER_DIR, PATENT_DIR
from pemt.patent_extractor.patent_chemical_harmonizer import harmonize_chemicals
from pemt.patent_extractor.patent_enrichment import extract_patent
from pemt.surechembl import SECTIONS, SureChEMBLBulk

logger = logging.getLogger(__name__)


@click.group()
def main():
    """Run PEMT."""
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


input_data = click.option(
    "--data",
    help="Path to tab-separated gene data file",
    type=click.Path(file_okay=True, dir_okay=False, exists=True),
    required=True,
)
input_data_type = click.option(
    "--input-type",
    help="Separator of the gene data file: 'comma' (csv), 'tab' (tsv) or 'semicolon'",
    type=click.Choice(["comma", "tab", "semicolon"], case_sensitive=False),
    default="comma",
    show_default=True,
)
analysis_name = click.option(
    "--name",
    help="Name of the analysis that is to be run",
    type=str,
    required=True,
)
has_uniprot = click.option(
    "--uniprot/--no-uniprot",
    default=True,
    help="Boolean value indicating whether the gene data file has uniprot ids or not.",
)
surechembl_source = click.option(
    "--surechembl-source",
    help=(
        "SureChEMBL bulk data to use: 'latest' (default, read from the EBI server), a release "
        "date such as 2026-09-22, a URL, or a local folder with the downloaded Parquet files"
    ),
    type=str,
    default="latest",
    show_default=True,
)
patent_sections = click.option(
    "--sections",
    help=(
        "Only count a chemical if it appears in these patent sections. Repeat the option for "
        "several sections, e.g. --sections claims --sections abstract. Default: any section"
    ),
    type=click.Choice(list(SECTIONS), case_sensitive=False),
    multiple=True,
)
patent_year = click.option(
    "--year",
    help="The year from which you want to retrive patents from",
    type=int,
    default=2000,
)
from_chemical = click.option(
    "--chemical/--no-chemical",
    default=True,
    help="Boolean value indicating whether the chemical data is provided by the user or not",
)
chemcial_data = click.option(
    "--chemical-data",
    help="Path to tab-separated chemical data file with single column of chembl_ids",
    type=click.Path(),
    default="",
)


@main.command(help="Extract chemicals for genes of interest")
@analysis_name
@input_data
@input_data_type
@has_uniprot
def run_chemical_extractor(
    name: str, data: str, input_type: str, uniprot: bool
) -> None:
    """Extracting chemicals for genes with experiemtal data."""
    click.echo(f"Starting the chemical extractor pipeline for {name}")

    if uniprot:
        with_uniprot = True
    else:
        with_uniprot = False

    gene_chemical_dict = extract_chemicals(
        analysis_name=name,
        gene_file_path=data,
        file_separator=input_type,
        is_uniprot=with_uniprot,
    )

    click.echo(
        f"Completed the chemical extractor pipeline for {len(gene_chemical_dict)} genes."
    )
    click.echo(f"Data file can be found under {MAPPER_DIR}")


@main.command(help="Extract patent for filtered chemicals")
@analysis_name
@patent_year
@surechembl_source
@patent_sections
@from_chemical
@chemcial_data
def run_patent_extractor(
    name: str,
    year: int,
    surechembl_source: str,
    sections: tuple,
    chemical: bool,
    chemical_data: str,
) -> None:
    """Extracting patent from chemical data."""
    click.echo("Starting to pre-process the chemical data for patent retrieval")

    with SureChEMBLBulk(surechembl_source) as bulk:
        if chemical:
            df = pd.read_csv(chemical_data, sep="\t", dtype=str)
            df.to_csv(f"{PATENT_DIR}/{name}_chemicals.tsv", sep="\t", index=False)
            harmonize_chemicals(analysis_name=name, from_genes=False, bulk=bulk)
        else:
            harmonize_chemicals(analysis_name=name, bulk=bulk)

        click.echo(f"Starting the patent extractor pipeline for {name}")
        patent_df = extract_patent(
            analysis_name=name,
            patent_year=year,
            sections=[s.lower() for s in sections] or None,
            bulk=bulk,
        )

    if patent_df.empty:
        click.echo("No patents found!")
        return None

    _save_patent_outputs(name, patent_df, with_genes=not chemical)


@main.command(help="Run the PEMT tool with gene data")
@analysis_name
@input_data
@input_data_type
@has_uniprot
@patent_year
@surechembl_source
@patent_sections
def run_pemt(
    name: str,
    data: str,
    input_type: str,
    uniprot: bool,
    year: int,
    surechembl_source: str,
    sections: tuple,
) -> None:
    """Runs the PEMT tool with all the components together."""
    click.echo(f"Starting to run PEMT workflow for {name}")
    click.echo("Running the chemical extractor pipeline")

    gene_chemical_dict = extract_chemicals(
        analysis_name=name,
        gene_file_path=data,
        file_separator=input_type,
        is_uniprot=bool(uniprot),
    )

    click.echo(
        f"Completed running the chemical extractor pipeline for {len(gene_chemical_dict)} genes."
    )

    with SureChEMBLBulk(surechembl_source) as bulk:
        click.echo("Pre-processing the chemical data for patent retrieval")
        harmonize_chemicals(analysis_name=name, bulk=bulk)

        click.echo("Running the patent extractor pipeline")
        patent_df = extract_patent(
            analysis_name=name,
            patent_year=year,
            sections=[s.lower() for s in sections] or None,
            bulk=bulk,
        )

    if patent_df.empty:
        click.echo("No patents found!")
        return None

    _save_patent_outputs(name, patent_df, with_genes=True)


def _save_patent_outputs(name: str, patent_df: pd.DataFrame, with_genes: bool) -> None:
    """For gene based runs, add the genes to the patent table and save it.

    The patent table itself (``<name>_patent_data.tsv``) is written by ``extract_patent``.
    """
    if with_genes:
        with open(f"{MAPPER_DIR}/{name}_gene_to_chemicals.json") as f:
            gene_chemical_data = json.load(f)

        chemical_to_gene_mapper = defaultdict(set)
        for gene, chemicals in gene_chemical_data.items():
            for chemical in chemicals:
                chemical_to_gene_mapper[chemical].add(gene)

        patent_df = patent_df.assign(
            genes=patent_df["chembl"].map(lambda x: ", ".join(sorted(chemical_to_gene_mapper[x])))
        )
        patent_df.to_csv(f"{PATENT_DIR}/{name}_gene_patent_data.tsv", sep="\t", index=False)

    click.echo("Done with retrieval of patents")
    click.echo(f"Data file can be found under {PATENT_DIR}")


if __name__ == "__main__":
    main()
