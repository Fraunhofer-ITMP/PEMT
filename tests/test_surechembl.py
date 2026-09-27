# -*- coding: utf-8 -*-

"""Tests for querying the SureChEMBL bulk data with DuckDB.

The fixture writes a tiny copy of the bulk data with the same schema as the real
Parquet files, so these tests run offline in a second. The one test marked ``live``
reads the real files from the EBI server; run it with ``uv run pytest -m live``.
"""

import datetime

import duckdb
import pandas as pd
import pytest

from pemt.constants import VALID_CODES
from pemt.surechembl import BULK_DATA_URL, PATENT_COLUMNS, SureChEMBLBulk, resolve_source

IMATINIB_KEY = "KTUFNOKKBVMGRW-UHFFFAOYSA-N"


@pytest.fixture(scope="module")
def bulk_dir(tmp_path_factory):
    """Write a miniature SureChEMBL release to a temporary folder."""
    path = tmp_path_factory.mktemp("surechembl")
    con = duckdb.connect()
    con.execute("""
        CREATE TABLE compounds (id BIGINT, smiles VARCHAR, inchi VARCHAR, inchi_key VARCHAR, mol_weight DOUBLE);
        INSERT INTO compounds VALUES
            (3827, 'C', 'InChI=1S/imatinib', 'KTUFNOKKBVMGRW-UHFFFAOYSA-N', 493.6),
            (9,    'C', 'InChI=1S/other',    'AAAAAAAAAAAAAA-UHFFFAOYSA-N', 100.0),
            (42,   'C', 'InChI=1S/nopat',    'BBBBBBBBBBBBBB-UHFFFAOYSA-N', 50.0);

        CREATE TABLE patents (
            id BIGINT, patent_number VARCHAR, country VARCHAR, publication_date DATE,
            family_id BIGINT, cpc VARCHAR[], ipcr VARCHAR[], ipc VARCHAR[], ecla VARCHAR[],
            assignee VARCHAR[], title VARCHAR);
        INSERT INTO patents VALUES
            -- valid IPC, recent: kept
            (1, 'WO-2011041462-A2', 'WO', DATE '2011-04-07', 100, [], ['A61P 35/00', 'H01L 21/00'], [], [], ['PHARMA A', 'BANK AS AGENT'], 't1'),
            -- only a non-PEMT IPC code: dropped by the IPC filter
            (2, 'EP-2405272-A1', 'EP', DATE '2012-01-11', 200, [], ['H01L 21/00'], [], [], ['ELECTRONICS B'], 't2'),
            -- valid IPC but published before 2000: dropped by the year filter
            (3, 'US-5153197-A', 'US', DATE '1992-10-06', 300, [], ['C07D 233/54'], [], [], ['DU PONT'], 't3'),
            -- valid IPC, recent, same family as patent 1
            (4, 'EP-2389136-A1', 'EP', DATE '2011-11-30', 100, [], ['C07D401/14'], [], [], ['PHARMA A'], 't4'),
            -- no IPCR codes at all
            (5, 'CN-104817498-A', 'CN', DATE '2015-07-29', 500, [], NULL, [], [], NULL, 't5');

        CREATE TABLE patent_compound_map (patent_id BIGINT, compound_id BIGINT, field_id BIGINT);
        INSERT INTO patent_compound_map VALUES
            (1, 3827, 1), (1, 3827, 2),   -- imatinib in description and claims
            (2, 3827, 1),
            (3, 3827, 2),
            (4, 3827, 1),                 -- description only
            (5, 3827, 1),
            (4, 9, 2);                    -- second compound, claims

        CREATE TABLE fields (id BIGINT, field_name VARCHAR);
        INSERT INTO fields VALUES (1,'desc'),(2,'clms'),(3,'abst'),(4,'ttl'),(5,'image'),(6,'molattachment');
    """)
    for table in ("compounds", "patents", "patent_compound_map", "fields"):
        con.execute(f"COPY {table} TO '{path / (table + '.parquet')}' (FORMAT parquet)")
    con.close()
    return path


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
