Welcome to PEMT's documentation!
===================================
**Release notes** : https://github.com/Fraunhofer-ITMP/PEMT/releases

.. raw:: html

   <h1 align="center">
     <img src='https://readthedocs.org/projects/pemt/badge/?version=latest' alt='Documentation Status' />
     <img src="https://img.shields.io/pypi/l/PEMT" alt="MIT">
     <img src="https://img.shields.io/pypi/v/PEMT" alt="PEMT on PyPI">
   </h1>

PEMT: A patent enrichment tool for drug discovery.
---------------------------------------------------------------------------------------------
PEMT takes a two-step approach to collect patent documents relevant for drug discovery.

1. The ``chemical_extractor`` module extraction of chemicals that directly regulate (i.e. activation or inhibition) genes of interest based on functional or biochemical assays found within ChEMBL.

2. The ``patent_extractor`` module interlinking these chemicals to patent documents by systematically querying SureChEMBL, a patent database.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction


.. toctree::
   :maxdepth: 3

   cli
   development


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
