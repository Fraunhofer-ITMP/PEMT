# -*- coding: utf-8 -*-

"""Command line interface."""

import json
import logging
import pandas as pd
from collections import defaultdict
from tqdm import tqdm

from pemt.constants import PATENT_DIR, MAPPER_DIR
from pemt.patent_extractor.patent_enrichment import extract_patent
from pemt.patent_extractor.patent_chemical_harmonizer import harmonize_chemicals
from pemt.chemical_extractor.experimental_data_extraction import extract_chemicals

import click
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@click.group()
def main():
    """Run PET."""
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")

input_data = click.option(
    '--data',
    help="Path to tab-separated gene data file",
    type=click.Path(file_okay=True, dir_okay=False, exists=True),
    required=True,
)
input_data_type = click.option(
    '--input-type',
    help="Type of data file i.e. 'tab' for tsv or 'comma' for csv files",
    type=str,
    default=','
)
analysis_name = click.option(
    '--name',
    help="Name of the analysis that is to be run",
    type=str,
    required=True,
)
has_uniprot = click.option(
    '--uniprot/--no-uniprot',
    default=True,
    help="Boolean value indicating whether the gene data file has uniprot ids or not.",
)
system_name = click.option(
    '--os',
    type=click.Choice(['linux', 'mac', 'windows'], case_sensitive=False),
    help="The OS system on which is the script is running",
    multiple=False
)
chromedriver_path = click.option(
    '--chromedriver-path',
    help="The path where the chromedriver can be found on the users computer",
    type=str,
    required=True,
)
patent_year = click.option(
    '--year',
    help="The year from which you want to retrive patents from",
    type=int,
    default=2000
)
from_chemical = click.option(
    '--chemical/--no-chemical',
    default=True,
    help="Boolean value indicating whether the chemical data is provided by the user or not",
)
chemcial_data = click.option(
    '--chemical-data',
    help="Path to tab-separated chemical data file with single column of chembl_ids",
    type=click.Path(),
    default=''
)

@main.command(help='Extract chemicals for genes of interest')
@analysis_name
@input_data
@input_data_type
@has_uniprot
def run_chemical_extractor(name: str, data: str, input_type: str, uniprot: bool) -> None:
    """Extracting chemicals for genes with experiemtal data."""
    click.echo(
        f"Starting the chemical extractor pipeline for {name}"
    )

    if uniprot:
        with_uniprot = True
    else:
        with_uniprot=False

    gene_chemical_dict = extract_chemicals(
        analysis_name=name,
        gene_file_path=data,
        file_separator=input_type,
        is_uniprot=with_uniprot
    )

    click.echo(f"Completed the chemical extractor pipeline for {len(gene_chemical_dict)} genes.")
    click.echo(f"Data file can be found under {MAPPER_DIR}")


@main.command(help='Extract patent for filtered chemicals')
@analysis_name
@system_name
@chromedriver_path
@patent_year
@from_chemical
@chemcial_data
def run_patent_extractor(
    name: str, os: str, chromedriver_path: str, year: str, chemical: bool, chemical_data: str
) -> None:
    """Extracting patent from chemical data."""
    click.echo(
        f"Starting to pre-process the chemical data for patent retrieval"
    )

    if chemical:
        df = pd.read_csv(chemical_data, sep='\t', dtype=str)

        df.to_csv(f'{PATENT_DIR}/{name}_chemicals.tsv', sep='\t', index=False)

        harmonize_chemicals(analysis_name=name, from_genes=False)
    else:
        harmonize_chemicals(analysis_name=name)

    click.echo(
        f"Starting the patent extractor pipeline for {name}"
    )

    patent_df = extract_patent(
        analysis_name=name,
        chrome_driver_path=chromedriver_path,
        os_system=os,
        patent_year=year
    )

    if patent_df.empty:
        click.echo(f'No patents found!')
        return None

    # Since the original patent data has chemical with no patents, we remove those entries from the data
    patent_df = patent_df[~patent_df['patent_id'].isna()]
    patent_df.to_csv(f'{PATENT_DIR}/cleaned_{name}_patent_data.tsv', sep='\t', index=False)

    if chemical:
        click.echo(f"Done with retrival of patents")
        click.echo(f"Data file can be found under {PATENT_DIR}")
        return

    gene_chemical_data = json.load(open(f'{MAPPER_DIR}/{name}_gene_to_chemicals.json'))
    chemical_to_gene_mapper = defaultdict(set)

    for gene, chemicals in tqdm(gene_chemical_data.items()):
        for chemical in chemicals:
            chemical_to_gene_mapper[chemical].add(gene)

    patent_df['genes'] = patent_df['chembl'].map(lambda x: ', '.join(chemical_to_gene_mapper[x]))

    patent_df.to_csv(f'{PATENT_DIR}/{name}_gene_patent_data.tsv', sep='\t', index=False)

    click.echo(f"Done with retrival of patents")
    click.echo(f"Data file can be found under {PATENT_DIR}")


@main.command(help='Run the PET tool with gene data')
@analysis_name
@input_data
@input_data_type
@has_uniprot
@chromedriver_path
@system_name
@patent_year
def run_pet(
    name: str, data: str, input_type: str, uniprot: bool, chromedriver_path: str, os: str, year: str
) -> None:
    """Runs the PET tool with all the components together."""
    click.echo(f"Starting to run PET workflow for {name}")

    click.echo(
        f"Running the chemical extractor pipeline"
    )

    if uniprot:
        with_uniprot = True
    else:
        with_uniprot=False

    gene_chemical_dict = extract_chemicals(
        analysis_name=name,
        gene_file_path=data,
        file_separator=input_type,
        is_uniprot=with_uniprot
    )

    click.echo(f"Completed running the chemical extractor pipeline for {len(gene_chemical_dict)} genes.")

    click.echo(f"Ppre-processing the chemical data for patent retrieval")

    harmonize_chemicals(analysis_name=name)

    click.echo(f"Running the patent extractor pipeline")

    patent_df = extract_patent(
        analysis_name=name,
        chrome_driver_path=chromedriver_path,
        os_system=os,
        patent_year=year
    )

    if patent_df.empty:
        click.echo(f'No patents found!')
        return None

    # Since the original patent data has chemical with no patents, we remove those entries from the data
    patent_df = patent_df[~patent_df['patent_id'].isna()]
    patent_df.to_csv(f'{PATENT_DIR}/cleaned_{name}_patent_data.tsv', sep='\t', index=False)

    chemical_to_gene_mapper = defaultdict(set)

    for gene, chemicals in tqdm(gene_chemical_data.items()):
        for chemical in chemicals:
            chemical_to_gene_mapper[chemical].add(gene)

    patent_df['genes'] = patent_df['chembl'].map(lambda x: ', '.join(chemical_to_gene_mapper[x]))

    patent_df.to_csv(f'{PATENT_DIR}/{name}_gene_patent_data.tsv', sep='\t', index=False)

    click.echo(f"Done with retrival of patents")


if __name__ == '__main__':
    main()
