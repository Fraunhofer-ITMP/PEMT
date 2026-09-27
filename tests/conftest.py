# -*- coding: utf-8 -*-

"""Shared test fixtures."""

import duckdb
import pytest


@pytest.fixture(scope="session")
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
