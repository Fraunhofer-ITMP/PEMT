# -*- coding: utf-8 -*-

"""Tests for querying the SureChEMBL bulk data with DuckDB.

The ``bulk_dir`` fixture (``conftest.py``) writes a tiny copy of the bulk data with the
same schema as the real Parquet files, so these tests run offline in a second. The one test marked ``live``
reads the real files from the EBI server; run it with ``uv run pytest -m live``.
"""

import datetime

import pandas as pd
import pytest

from pemt.constants import VALID_CODES
from pemt.surechembl import BULK_DATA_URL, PATENT_COLUMNS, SureChEMBLBulk, resolve_source

IMATINIB_KEY = "KTUFNOKKBVMGRW-UHFFFAOYSA-N"


@pytest.fixture()
def bulk(bulk_dir):
    with SureChEMBLBulk(str(bulk_dir), progress_bar=False) as b:
        yield b


def test_resolve_source(bulk_dir):
    assert resolve_source(None) == BULK_DATA_URL + "latest/"
    assert resolve_source("latest") == BULK_DATA_URL + "latest/"
    assert resolve_source("2026-09-22") == BULK_DATA_URL + "2026-09-22/"
    assert resolve_source("https://example.org/x") == "https://example.org/x/"
    assert resolve_source(str(bulk_dir)).startswith(str(bulk_dir))


def test_resolve_source_missing_files(tmp_path):
    with pytest.raises(FileNotFoundError, match="patents.parquet"):
        resolve_source(str(tmp_path))


def test_compounds_for_inchikeys(bulk):
    df = bulk.compounds_for_inchikeys([IMATINIB_KEY, "NOTINDATA-UHFFFAOYSA-N", None])
    assert df.to_dict("records") == [{"inchi_key": IMATINIB_KEY, "surechembl_id": 3827}]
    assert bulk.compounds_for_inchikeys([]).empty


def test_all_patents_without_filters(bulk):
    df = bulk.patents_for_compounds([3827])
    assert list(df.columns) == PATENT_COLUMNS
    assert set(df["patent_number"]) == {
        "WO-2011041462-A2", "EP-2405272-A1", "US-5153197-A", "EP-2389136-A1", "CN-104817498-A",
    }


def test_pemt_filters(bulk):
    df = bulk.patents_for_compounds([3827], ipc_prefixes=VALID_CODES, min_year=2000)
    assert list(df["patent_number"]) == ["WO-2011041462-A2", "EP-2389136-A1"]

    first = df.iloc[0]
    assert first["ipc"] == "A61P 35/00"  # only the matching code is kept
    assert first["ipc"].split()[0] in VALID_CODES  # what visualization.py relies on
    assert first["assignee"] == "PHARMA A"
    assert first["sections"] == "description; claims"
    assert first["family_id"] == 100
    assert pd.Timestamp(first["publication_date"]).date() == datetime.date(2011, 4, 7)


def test_ipc_prefix_ignores_spaces_and_length(bulk):
    df = bulk.patents_for_compounds([3827], ipc_prefixes=["C07D401"])
    assert list(df["patent_number"]) == ["EP-2389136-A1"]


def test_sections_filter(bulk):
    df = bulk.patents_for_compounds([3827], sections=["claims"])
    assert set(df["patent_number"]) == {"WO-2011041462-A2", "US-5153197-A"}
    assert set(df["sections"]) == {"claims"}


def test_several_compounds_and_id_formats(bulk):
    df = bulk.patents_for_compounds(["SCHEMBL3827", 9, "42"], ipc_prefixes=VALID_CODES, min_year=2000)
    assert df.groupby("compound_id")["patent_number"].apply(list).to_dict() == {
        9: ["EP-2389136-A1"],
        3827: ["WO-2011041462-A2", "EP-2389136-A1"],
    }


def test_empty_and_invalid_input(bulk):
    assert list(bulk.patents_for_compounds([]).columns) == PATENT_COLUMNS
    with pytest.raises(ValueError, match="Unknown sections"):
        bulk.patents_for_compounds([3827], sections=["claim"])
    with pytest.raises(ValueError, match="Invalid IPC prefix"):
        bulk.patents_for_compounds([3827], ipc_prefixes=["A61P'; DROP"])


@pytest.mark.live
def test_live_latest_release_schema():
    """Reads only Parquet footers from the EBI server (a few seconds)."""
    with SureChEMBLBulk(progress_bar=False) as b:
        for table, required in {
            "compounds": {"id", "inchi_key"},
            "patents": {"id", "patent_number", "publication_date", "family_id", "ipcr", "assignee"},
            "patent_compound_map": {"patent_id", "compound_id", "field_id"},
        }.items():
            columns = {row[0] for row in b.con.execute(f"DESCRIBE {table}").fetchall()}
            assert required <= columns, f"{table} is missing {required - columns}"
