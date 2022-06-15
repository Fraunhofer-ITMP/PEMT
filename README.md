<p align="center">
  <img src="../PEMT/docs/source/logoIME2.png">
</p>

<h1 align="center">
  PEMT: A tool for extracting patent literature in drug discovery
  <br/>

[comment]: <> (  <a href='https://travis-ci.com/github/hybrid-kg'>)

[comment]: <> (     <img src="https://travis-ci.com/hybrid-kg/clep.svg?branch=master" />)

[comment]: <> (  </a>)

[comment]: <> (  <a href='https://clep.readthedocs.io/en/latest/?badge=latest'>)

[comment]: <> (    <img src='https://readthedocs.org/projects/clep/badge/?version=latest' alt='Documentation Status' />)

[comment]: <> (  </a>)

[comment]: <> (  <a href="https://zenodo.org/badge/latestdoi/209278408">)

[comment]: <> (    <img src="https://zenodo.org/badge/209278408.svg" alt="DOI">)

[comment]: <> (  </a>)

[comment]: <> (  <a href="https://pypi.org/project/clep/">)

[comment]: <> (    <img src="https://img.shields.io/pypi/v/clep" alt="CLEP on PyPI">)

[comment]: <> (  </a>)

[comment]: <> (  <img src="https://img.shields.io/pypi/pyversions/clep" alt="CLEP Python versions">)

[comment]: <> (  <a href="https://github.com/Fraunhofer-ITMP/PEMT/blob/master/LICENSE">)

[comment]: <> (    <img src="https://img.shields.io/pypi/l/clep" alt="MIT">)

[comment]: <> (  </a>)
</h1>

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
  <img src="../PEMT/docs/source/framework.jpg">
</p>

## Installation

[comment]: <> (The code can be installed from [PyPI]&#40;https://pypi.org/project/clep/&#41; with:)

```bash
$ pip install pemt
```

The most recent code can be installed from the source on [GitHub](https://github.com/Fraunhofer-ITMP/PEMT) with:

```bash
$ pip install git+https://github.com/Fraunhofer-ITMP/PEMT.git
```

Alternatively, for developer the tool can be installed in an editable mode as shown below:

```bash
$ git clone https://github.com/Fraunhofer-ITMP/PEMT.git
$ conda create --name pemt python=3.8
$ conda activate pemt
$ cd PEMT
$ pip install pemt
```

For developers, the repository can be cloned from [GitHub](https://github.com/Fraunhofer-ITMP/PEMT) and installed in editable mode with:

```bash
$ git clone https://github.com/Fraunhofer-ITMP/PEMT.git
$ cd PEMT
$ pip install -e .
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


## Usage

In-order to use PEMT, an installation of [**chromedriver**](https://chromedriver.chromium.org/) is **required**.

As mentioned above, the tool has a two-step approach. Each of these steps can be run individually as well as together as show belwo:

1. **Chemical enrichment**
The following command links chemicals to genes of interest based on causality. In this command it is necessary to indicate whether the file contains uniprot ids or not with the `--uniprot` or `--no-uniprot` parameter.

```bash
$ pemt run-chemical-extractor --name=<ANALYSIS NAME> --data=<DATA FILE PATH> --input-type=<DATA FILE SEPARATOR> --uniprot

```

2. **Patent enrichment**
The following command interlinks chemicals to patent literature publicly available. 

```bash
$ pemt run_patent_extractor --data --name=<ANALYSIS NAME> --data=<DATA FILE PATH> --input-type=<DATA FILE SEPARATOR> --chromedriver-path=<PATH TO CHROMEDRIVER> --os=<OS NAME>
```

We also allow the flexibility to start the pipeline from this step, if the user has list of chemicals in the right format as indicated above. The user then has to use the tag `--chemical` and provide a respective `--chemical-data` path.

3. **PEMT workflow**
The following command generates the patent enrichment on the gene data where the gene data file is a TSV file containing uniprot identifiers.

```bash
$ pemt run-pemt --name=<ANALYSIS NAME> --data=<DATA FILE PATH> --input-type=<DATA FILE SEPARATOR> --chromedriver-path=<PATH TO CHROMEDRIVER> --os=<OS NAME>
```

## Issues

If you have difficulties using PEMT, please open an issue at our [GitHub](https://github.com/Fraunhofer-ITMP/PEMT) repository.

[comment]: <> (### Citation)

[comment]: <> (If you have found PEMT useful in your work, please consider citing:)

[comment]: <> ([**CLEP: A Hybrid Data- and Knowledge- Driven Framework for Generating Patient Representations**]&#40;https://doi.org/10.1093/bioinformatics/btab340)

[comment]: <> (&#41;.<br />)

[comment]: <> (Bharadhwaj, V. S., Ali, M., Birkenbihl, C., Mubeen, S., Lehmann, J., Hofmann-Apitius, M., Hoyt, C. T., & Domingo-Fernandez, D. &#40;2020&#41;.<br />)

[comment]: <> (*Bioinformatics*, btab340. )

[comment]: <> (### Graphics)

[comment]: <> (The CLEP logo and framework graphic was designed by Carina Steinborn.)

## Disclaimer

PEMT is a scientific tool that has been developed in an academic capacity, and thus comes with no warranty or guarantee of maintenance, support, or back-up of data.
