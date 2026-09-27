# -*- coding: utf-8 -*-

"""Tests for mapping gene symbols to UniProt (HGNC) and recording the ChEMBL release."""

import json
import os
import time

import pandas as pd
import pytest

import pemt.chemical_extractor.experimental_data_extraction as extractor
from pemt.utils import load_hgnc, symbols_to_uniprot

UNREACHABLE = "http://127.0.0.1:9/hgnc_complete_set.txt"  # fails immediately

HGNC_ROWS = [
    # hgnc_id, symbol, status, prev_symbol, alias_symbol, uniprot_ids
    ("HGNC:76", "ABL1", "Approved", "", "JTK7|p150", "P00519"),
    ("HGNC:7132", "KMT2A", "Approved", "MLL|TRX1", "HRX|ALL-1", "Q03164"),
    ("HGNC:1437", "CALCA", "Approved", "CALC1", "CGRP", "P01258|P06881"),
    ("HGNC:1", "GENEX", "Approved", "", "SHARED", "P11111"),
    ("HGNC:2", "GENEY", "Approved", "", "SHARED", "P22222"),
    ("HGNC:3", "MIR21", "Approved", "", "", ""),  # non-coding: no protein
    ("HGNC:4", "OLDGONE", "Entry Withdrawn", "", "", "P99999"),
]


@pytest.fixture()
def hgnc_dir(tmp_path):
    df = pd.DataFrame(
        HGNC_ROWS,
        columns=[
            "hgnc_id",
            "symbol",
            "status",
            "prev_symbol",
            "alias_symbol",
            "uniprot_ids",
        ],
    )
    df["name"] = "x"  # extra columns in the real file are ignored
    df.to_csv(tmp_path / "hgnc_complete_set.txt", sep="\t", index=False)
    return tmp_path


def test_load_hgnc_uses_fresh_cache(hgnc_dir):
    df = load_hgnc(str(hgnc_dir), url=UNREACHABLE)  # no download needed
    assert list(df.columns) == [
        "hgnc_id",
        "symbol",
        "status",
        "prev_symbol",
        "alias_symbol",
        "uniprot_ids",
    ]
    assert len(df) == len(HGNC_ROWS)


def test_load_hgnc_falls_back_to_old_copy(hgnc_dir, caplog):
    old = time.time() - 90 * 86400
    os.utime(hgnc_dir / "hgnc_complete_set.txt", (old, old))
    df = load_hgnc(str(hgnc_dir), url=UNREACHABLE)
    assert len(df) == len(HGNC_ROWS)
    assert "using the 90 day old copy" in caplog.text


def test_load_hgnc_without_copy_or_network(tmp_path):
    with pytest.raises(RuntimeError, match="Could not download"):
        load_hgnc(str(tmp_path), url=UNREACHABLE)


def test_symbols_to_uniprot(hgnc_dir):
    hgnc = load_hgnc(str(hgnc_dir), url=UNREACHABLE)
    mapping = symbols_to_uniprot(
        [
            "ABL1",
            "abl1 ",
            "MLL",
            "CGRP",
            "CALCA",
            "SHARED",
            "MIR21",
            "OLDGONE",
            "NOTAGENE",
        ],
        hgnc,
    )
    assert mapping == {
        "ABL1": ["P00519"],
        "abl1 ": ["P00519"],  # case and whitespace are ignored
        "MLL": ["Q03164"],  # previous symbol of KMT2A
        "CGRP": ["P01258", "P06881"],  # alias, gene with two proteins
        "CALCA": ["P01258", "P06881"],
        "SHARED": [],  # alias of two genes: ambiguous, skipped
        "MIR21": [],  # no protein
        "OLDGONE": [],  # withdrawn entries are not used
        "NOTAGENE": [],
    }


class FakeQuery(list):
    def only(self, fields):
        return self


class FakeActivity:
    DATA = {
        "CHEMBL1862": [{"molecule_chembl_id": "CHEMBL941", "pchembl_value": "8.1"}],
        "CHEMBL_CALCA_1": [{"molecule_chembl_id": "CHEMBL_A", "pchembl_value": "7"}],
        "CHEMBL_CALCA_2": [{"molecule_chembl_id": "CHEMBL_B", "pchembl_value": "7"}],
    }

    def filter(self, **kwargs):
        return FakeQuery(self.DATA.get(kwargs["target_chembl_id"], []))


@pytest.fixture()
def fake_chembl(monkeypatch, hgnc_dir):
    targets = {
        "P00519": ["CHEMBL1862"],
        "P01258": ["CHEMBL_CALCA_1"],
        "P06881": ["CHEMBL_CALCA_2"],
    }
    monkeypatch.setattr(
        extractor, "get_single_protein_targets", lambda u: targets.get(u, [])
    )
    monkeypatch.setattr(extractor, "activity", FakeActivity())
    monkeypatch.setattr(extractor, "get_chembl_release", lambda: "ChEMBL_TEST")
    monkeypatch.setattr(
        extractor, "load_hgnc", lambda d: load_hgnc(str(hgnc_dir), url=UNREACHABLE)
    )
    monkeypatch.setattr(extractor, "MAPPER_DIR", str(hgnc_dir))
    return hgnc_dir


def test_gene_with_several_proteins_uses_all_of_them(fake_chembl):
    assert extractor.target_to_chemical(
        "CALCA", protein_mapping={"CALCA": ["P01258", "P06881"]}
    ) == [
        "CHEMBL_A",
        "CHEMBL_B",
    ]
    # the old table stored them as one string; that still works
    assert extractor.target_to_chemical(
        "CALCA", protein_mapping={"CALCA": "P01258, P06881"}
    ) == [
        "CHEMBL_A",
        "CHEMBL_B",
    ]


def test_extract_chemicals_from_symbols(fake_chembl):
    result = extractor.extract_chemicals(
        "sym", gene_list=["ABL1", "MLL", "CALCA", "NOTAGENE"]
    )
    assert dict(result) == {
        "ABL1": ["CHEMBL941"],
        "CALCA": ["CHEMBL_A", "CHEMBL_B"],
        "MLL": [],
        "NOTAGENE": [],
    }
    info = json.loads((fake_chembl / "sym_run_info.json").read_text())
    assert info["chembl_release"] == "ChEMBL_TEST"
    assert info["input_type"] == "symbol"
    assert info["symbol_to_uniprot"]["MLL"] == ["Q03164"]
    assert "hgnc_file_date" in info


@pytest.mark.live
def test_live_hgnc_download(tmp_path):
    hgnc = load_hgnc(str(tmp_path))
    mapping = symbols_to_uniprot(["ABL1", "MLL", "KIT"], hgnc)
    assert mapping["ABL1"] == ["P00519"]
    assert mapping["MLL"] == ["Q03164"]
    assert mapping["KIT"] == ["P10721"]


@pytest.mark.live
def test_live_chembl_release():
    from pemt.utils import get_chembl_release

    assert get_chembl_release().startswith("ChEMBL_")
