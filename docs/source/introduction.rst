.. _intro:

.. _GitHub: https://github.com/Fraunhofer-ITMP/PEMT

Welcome to PEMT's documentation!
===================================
**Release notes** : https://github.com/Fraunhofer-ITMP/PEMT/releases

.. image:: ./logo.jpg
    :align: center

.. raw:: html

   <h1 align="center">
     <img src="https://github.com/Fraunhofer-ITMP/PEMT/actions/workflows/tests.yml/badge.svg" alt="Tests" />
     <img src='https://readthedocs.org/projects/pemt/badge/?version=latest' alt='Documentation Status' />
     <img src='https://img.shields.io/github/license/Fraunhofer-ITMP/PEMT?color=blue' alt='GitHub License' />
   </h1>

PEMT: A patent enrichment tool for drug discovery.
---------------------------------------------------------------------------------------------
PEMT takes a two-step approach to collect patent documents relevant for drug discovery.

1. The ``chemical_extractor`` module extraction of chemicals that directly regulate (i.e. activation or inhibition) genes of interest based on functional or biochemical assays found within ChEMBL.

2. The ``patent_extractor`` module interlinking these chemicals to patent documents in SureChEMBL, a patent database, using its bulk data release (queried with DuckDB).

General info
-------------
PEMT is a patent extractor tool that enables users to retrieve patents relevant to drug discovery. The framework is depicted in the graphic below

.. image:: ./framework.jpg
    :align: center


Installation
------------

You can install PEMT package from pypi.

.. code:: shell

   # Use pip to install the latest release
   $ python3 -m pip install pemt

You may instead want to use the development version from Github, by running

.. code:: shell

   $ python3 -m pip install git+https://github.com/Fraunhofer-ITMP/PEMT.git

For contributors, the repository can be cloned from `GitHub`_ and set up with `uv <https://docs.astral.sh/uv/>`_:

.. code:: shell

   $ git clone https://github.com/Fraunhofer-ITMP/PEMT.git
   $ cd PEMT
   $ uv sync
   $ uv run pytest

Dependency
--------------
- Python 3.10+

Mandatory
~~~~~~~~~

- Pandas
- CheMBL Webresource
- DuckDB (to query the SureChEMBL bulk data)


For API information to use this library, see the :ref:`dev-guide`.

Issues
-------

If you have difficulties using PEMT, please open an issue at our `GitHub`_ repository.


Disclaimer
-----------

PEMT is a scientific tool that has been developed in an academic capacity, and thus comes with no warranty or guarantee of maintenance, support, or back-up of data.
