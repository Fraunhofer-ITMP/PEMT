.. _dev-guide:

Developmental Guide
=====================================

.. module:: pemt

PEMT runs in two steps. The chemical extractor finds active chemicals for genes in
ChEMBL; the patent extractor maps those chemicals to SureChEMBL and retrieves the
patents that mention them from the SureChEMBL bulk data.


Core Module APIs
-----------------

Gene and Chemical Mapping
~~~~~~~~~~~~~~~~~~~~~~~~~

Gene symbols are mapped to UniProt accessions with the HGNC complete set, UniProt
accessions to ChEMBL single protein targets, and ChEMBL compounds to their names and
InChIKeys.

.. autofunction:: pemt.utils.load_hgnc

.. autofunction:: pemt.utils.symbols_to_uniprot

.. autofunction:: pemt.utils.get_single_protein_targets

.. autofunction:: pemt.utils.get_chembl_structures

.. autofunction:: pemt.utils.get_chembl_release

Chemical Extractor
~~~~~~~~~~~~~~~~~~

.. autofunction:: pemt.chemical_extractor.experimental_data_extraction.extract_chemicals

.. autofunction:: pemt.chemical_extractor.experimental_data_extraction.target_to_chemical

Chemical Harmonizer
~~~~~~~~~~~~~~~~~~~

.. autofunction:: pemt.patent_extractor.patent_chemical_harmonizer.harmonize_chemicals

Patent Extractor
~~~~~~~~~~~~~~~~

.. autofunction:: pemt.patent_extractor.patent_enrichment.extract_patent

SureChEMBL Bulk Data
~~~~~~~~~~~~~~~~~~~~

Patents are retrieved from the SureChEMBL bulk data (Apache Parquet), queried with
DuckDB either directly from the EBI server or from a local folder.

.. automodule:: pemt.surechembl
   :no-members:

.. autofunction:: pemt.surechembl.resolve_source

.. autoclass:: pemt.surechembl.SureChEMBLBulk
   :members: compounds_for_inchikeys, patents_for_compounds, close
