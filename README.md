<p align="center">
  <img style="width: 200px; height: 200px;" src="docs/PEMT%20Logo.jpg">
</p>


<h1 align="center">
  PEMT: A tool for extracting patent literature in drug discovery
  <br/>
    <a href='https://pemt.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/pemt/badge/?version=latest' alt='Documentation Status' />
    </a>
    <a href="https://pypi.org/project/PEMT/">
        <img src="https://img.shields.io/pypi/v/PEMT" alt="PEMT on PyPI">
    </a>
    <a href="https://github.com/Fraunhofer-ITMP/PET/blob/master/LICENSE">
        <img src="https://img.shields.io/pypi/l/PEMT" alt="MIT">
    </a>

 [![DOI:10.1093/bioinformatics/btac716](http://img.shields.io/badge/DOI-110.1093/bioinformatics/btac716-B31B1B.svg)](https://doi.org/10.1093/bioinformatics/btac716)

</h1>

> [!WARNING]
> **Currently SureCHEMBL has undergone a major restructuring. The tool might not function in that case. The tool will be updated in soon!!**

## Table of Contents

* [General Info](#general-info)
* [Installation](#installation)
* [Documentation](#documentation)
* [Input Data](#input-data-formats)
* [Usage](#usage)
* [Issues](#issues)
* [Disclaimer](#disclaimer)

## General Info

PEMT is a patent extractor tool that enables users to retrieve patents relevant to drug discovery. The overall workflow of the tool can be seen in the figure below:

<p align="center">
  <img src="docs/source/framework.jpg">
</p>

## Installation

[comment]: <> (The code can be installed from [PyPI]&#40;https://pypi.org/project/clep/&#41; with:)

```shell
$ pip install pemt
```

The most recent code can be installed from the source on [GitHub](https://github.com/Fraunhofer-ITMP/PEMT) with:

```shell
$ pip install git+https://github.com/Fraunhofer-ITMP/PEMT.git
```

For developers, the project uses [uv](https://docs.astral.sh/uv/) to manage the environment. Clone the repository and create the environment (Python version from `.python-version`, dependencies from `uv.lock`) with:

```shell
$ git clone https://github.com/Fraunhofer-ITMP/PEMT.git
$ cd PEMT
$ uv sync
```

Run the command line tool and the tests inside that environment with:

```shell
$ uv run pemt --help
$ uv run pytest
```

## Documentation

Read the [official docs](https://pemt.readthedocs.io/en/latest/) for more information.

## Input Data Formats

### Data

For running PEMT from the gene level, you need the input file with the following structure:

| symbol | uniprot |
| ------ | -------- |
| HGNC_Symbol_1 | Uniprot_ID_1
| HGNC_Symbol_2 | Uniprot_ID_2
| HGNC_Symbol_3 | Uniprot_ID_3

For running PEMT from the chemical level, you need the input file with the following structure:

| chembl |
| ------ |
| ChEMBL_ID_1
| ChEMBL_ID_2
| ChEMBL_ID_3

**Note:** The data must be in a comma or tab separated file format. If not so, the file should have at least one of the columns shown above.

**How genes are mapped to ChEMBL:**

- **Gene symbols** (`--no-uniprot`) are resolved with the [HGNC complete set](https://www.genenames.org/download/statistics-and-files/), which PEMT downloads and refreshes every 30 days (`data/mapper/hgnc_complete_set.txt`). A symbol is matched as an approved symbol first, then as a previous symbol (e.g. `MLL` → `KMT2A`), then as an alias; ambiguous aliases are skipped. Genes with several protein products use all of their UniProt accessions.
- **UniProt accessions** are used as given.
- Each accession is looked up in the current ChEMBL release and only its *single protein* target is used, not complexes, fusion proteins or degrader targets that contain the protein.

The ChEMBL release, the HGNC file date and the symbol → UniProt mapping of each run are saved in `data/mapper/<ANALYSIS NAME>_run_info.json`.


## Usage

Patents are retrieved from the [SureChEMBL bulk data](https://chembl.gitbook.io/surechembl/downloads/bulk-data), which is queried with [DuckDB](https://duckdb.org/). By default PEMT reads the latest release directly from the EBI server, so no download is needed; looking up the patents of a run takes roughly 15 minutes regardless of how many chemicals it has. Use `--surechembl-source` to pick a dated release (e.g. `2026-09-22`, for reproducible results) or a local folder with the downloaded Parquet files.

As mentioned above, the tool has a two-step approach. Each of these steps can be run individually as well as together as show belwo:

1. **Chemical enrichment**
The following command links chemicals to genes of interest based on causality. In this command it is necessary to indicate whether the file contains uniprot ids or not with the `--uniprot` or `--no-uniprot` parameter.

```shell
$ pemt run-chemical-extractor --name=<ANALYSIS NAME> --data=<DATA FILE PATH> --input-type=<DATA FILE SEPARATOR> --uniprot

```

2. **Patent enrichment**
The following command interlinks chemicals to patent literature publicly available.

```shell
$ pemt run-patent-extractor --name=<ANALYSIS NAME> --no-chemical
```

By default a chemical counts for a patent wherever it is mentioned. To only count patents that name the chemical in specific sections, add `--sections` (repeatable), e.g. `--sections claims`.

> **Tip:** most mentions of a chemical are in a patent's description only (about 88% in a test run on ABL1, KIT and FLT3). For well-known drugs and common compounds (e.g. imatinib, or oleic acid, which is active on FLT3 in ChEMBL and appears in over 400,000 patents) this is dominated by passing mentions such as combination therapies or formulations. Use `--sections claims` for a patent landscape of chemicals that are actually claimed.

We also allow the flexibility to start the pipeline from this step, if the user has list of chemicals in the right format as indicated above. The user then has to use the tag `--chemical` and provide a respective `--chemical-data` path.

3. **PEMT workflow**
The following command generates the patent enrichment on the gene data where the gene data file is a TSV file containing uniprot identifiers.

```shell
$ pemt run-pemt --name=<ANALYSIS NAME> --data=<DATA FILE PATH> --input-type=<DATA FILE SEPARATOR>
```

## Issues

If you have difficulties using PEMT, please open an issue at our [GitHub](https://github.com/Fraunhofer-ITMP/PEMT) repository.

## Citation

If you have found PEMT useful in your work, please consider citing: [**PEMT: A patent enrichment tool for drug discovery**](https://doi.org/10.1093/bioinformatics/btac716).

> Yojana Gadiya, Andrea Zaliani, Philip Gribbon, Martin Hofmann-Apitius, PEMT: a patent enrichment tool for drug discovery, *Bioinformatics*, 2022;, btac716, [https://doi.org/10.1093/bioinformatics/btac716](https://doi.org/10.1093/bioinformatics/btac716)

## Disclaimer

PEMT is a scientific tool that has been developed in an academic capacity, and thus comes with no warranty or guarantee of maintenance, support, or back-up of data.


## Funding
This project has been funded by EOSC-Life which has received funding from the European Union's Horizon 2020 programme under grant agreement number 824087.
