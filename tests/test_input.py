# -*- coding: utf-8 -*-

"""Tests for CDF patient incorporation."""

import os
import unittest

from pemt.chemical_extractor.experimental_data_extraction import extract_chemicals
from pemt.constants import MAPPER_DIR

TEST_FOLDER = os.path.dirname(os.path.realpath(__file__))
DUMMY_DATA = os.path.join(TEST_FOLDER, 'resources', 'dummy_gene.tsv')


class TestInput(unittest.TestCase):
    """Tests for input method."""

    def test_chemical_extractor(self):
        """Test chemical extractor."""

        output = extract_chemicals(
            analysis_name='test',
            gene_file_path=DUMMY_DATA,
            file_separator='comma',
            is_uniprot=True
        )

        """Test dummy protein are returned"""
        self.assertEqual(
            len(output),
            11
        )

        """Test dummy protein are outputted in a file"""
        self.assertTrue(os.path.exists(f'{MAPPER_DIR}/test_gene_to_chemicals.json'))