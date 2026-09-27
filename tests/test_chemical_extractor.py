# -*- coding: utf-8 -*-

"""Tests for the ChEMBL chemical extractor (ChEMBL is faked except in ``live`` tests)."""

import json

import pytest

import pemt.chemical_extractor.experimental_data_extraction as extractor


class FakeQuery(list):
    def only(self, fields):
        return self


class FakeActivity:
    """Mimics new_client.activity for two targets."""

    DATA = {
        "CHEMBL1862": [  # ABL1
            {"molecule_chembl_id": "CHEMBL941", "pchembl_value": "8.1"},
            {"molecule_chembl_id": "CHEMBL941", "pchembl_value": "7.5"},  # duplicate compound
            {"molecule_chembl_id": "CHEMBL_WEAK", "pchembl_value": "5.2"},  # below 6
            {"molecule_chembl_id": "CHEMBL_NOVALUE", "pchembl_value": None},
        ],
        "CHEMBL1936": [{"molecule_chembl_id": "CHEMBL941", "pchembl_value": "6.0"}],  # KIT
    }

    def __init__(self):
        self.calls = []

    def filter(self, **kwargs):
        self.calls.append(kwargs)
        return FakeQuery(self.DATA.get(kwargs["target_chembl_id"], []))


@pytest.fixture()
def fake_chembl(monkeypatch):
    targets = {"P00519": ["CHEMBL1862"], "P10721": ["CHEMBL1936"]}
    monkeypatch.setattr(extractor, "get_single_protein_targets", lambda u: targets.get(u, []))
    fake = FakeActivity()
    monkeypatch.setattr(extractor, "activity", fake)
    return fake


def test_target_to_chemical_uniprot(fake_chembl):
    assert extractor.target_to_chemical("P00519", is_uniprot=True) == ["CHEMBL941"]
    call = fake_chembl.calls[0]
    assert call["target_chembl_id"] == "CHEMBL1862"
    assert call["pchembl_value__gte"] == 6


def test_target_to_chemical_symbol(fake_chembl):
    assert extractor.target_to_chemical("KIT", protein_mapping={"KIT": "P10721"}) == ["CHEMBL941"]
    assert extractor.target_to_chemical("NOTAGENE", protein_mapping={}) == []


def test_target_to_chemical_needs_mapping_for_symbols(fake_chembl):
    with pytest.raises(ValueError, match="HGNC symbol"):
        extractor.target_to_chemical("KIT")


def test_chemical_overview_when_every_gene_has_chemicals(tmp_path):
    path = tmp_path / "x.json"
    path.write_text(json.dumps({"ABL1": ["CHEMBL941"], "KIT": ["CHEMBL941"]}))
    extractor.get_chemical_overview(str(path))  # used to raise KeyError


@pytest.mark.live
def test_live_single_protein_targets():
    """ABL1, KIT and FLT3 must map to the protein itself, not fusion/degrader targets."""
    from pemt.utils import get_single_protein_targets

    assert get_single_protein_targets("P00519") == ["CHEMBL1862"]
    assert get_single_protein_targets("P10721") == ["CHEMBL1936"]
    assert get_single_protein_targets("P36888") == ["CHEMBL1974"]
    assert get_single_protein_targets("Q8V4Y0") == []
