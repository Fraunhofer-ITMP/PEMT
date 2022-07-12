.. _dev-guide:

Developmental Guide
=====================================

.. module:: pemt


Core Module APIs
-----------------

Gene Harmonizer
~~~~~~~~~~~~~~~

.. automodule:: pemt.utils
   :members:

.. autofunction:: pemt.utils.uniprot_to_chembl()

.. autofunction:: pemt.utils.hgnc_to_chembl()

Chemical Extractor
~~~~~~~~~~~~~~~

.. autofunction:: pemt.chemical_extractor.experimental_data_extraction.extract_chemicals()

Chemical Harmonizer
~~~~~~~~~~~~~~~

.. autofunction:: pemt.patent_extractor.patent_chemical_harmonizer.harmonize_chemicals()


Patent Extractor
~~~~~~~~~~~~~~

.. autofunction:: pemt.patent_extractor.patent_enrichment.extract_patent()


