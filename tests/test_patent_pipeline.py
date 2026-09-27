# -*- coding: utf-8 -*-

"""Tests for the patent part of PEMT: harmonizer, patent extractor and CLI.

ChEMBL is replaced by a fake lookup and SureChEMBL by the miniature release from
``conftest.py``, so these tests run offline.
"""

import json

import pandas as pd
import pytest
from click.testing import CliRunner

import pemt.cli
import pemt.patent_extractor.patent_chemical_harmonizer as harmonizer
import pemt.patent_extractor.patent_enrichment as enrichment
from pemt.patent_extractor.patent_enrichment import PATENT_COLUMNS
from pemt.surechembl import SureChEMBLBulk

FAKE_CHEMBL = pd.DataFrame(
    [
        ("CHEMBL941", "IMATINIB", "KTUFNOKKBVMGRW-UHFFFAOYSA-N"),  # SureChEMBL 3827
        ("CHEMBL_OTHER", "OTHER", "AAAAAAAAAAAAAA-UHFFFAOYSA-N"),  # SureChEMBL 9
        ("CHEMBL_NOTINSCHEMBL", "NOT IN SURECHEMBL", "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N"),
        ("CHEMBL_BIOLOGIC", "BIOLOGIC", None),  # no structure, no InChIKey
    ],
    columns=["chembl", "name", "inchi_key"],
)


@pytest.fixture()
def chembl_calls():
    """Records every (fake) ChEMBL lookup made during a test."""
    return []


@pytest.fixture()
def workdir(tmp_path, monkeypatch, chembl_calls):
    """Point PEMT's data folders at a temporary directory and fake the ChEMBL lookup."""
    for module in (harmonizer, enrichment, pemt.cli):
        monkeypatch.setattr(module, "PATENT_DIR", str(tmp_path), raising=False)
        monkeypatch.setattr(module, "MAPPER_DIR", str(tmp_path), raising=False)

    def fake_structures(ids):
        chembl_calls.append(sorted(ids))
        return FAKE_CHEMBL[FAKE_CHEMBL["chembl"].isin(ids)].reset_index(drop=True)

    monkeypatch.setattr(harmonizer, "get_chembl_structures", fake_structures)
    return tmp_path


def write_gene_file(workdir, name="run"):
    genes = {"ABL1": ["CHEMBL941", "CHEMBL_BIOLOGIC"], "KIT": ["CHEMBL941", "CHEMBL_OTHER"], "EMPTY": []}
    (workdir / f"{name}_gene_to_chemicals.json").write_text(json.dumps(genes))


def test_harmonize_from_genes(workdir, bulk_dir, chembl_calls):
    write_gene_file(workdir)
    df = harmonizer.harmonize_chemicals("run", source=str(bulk_dir))

    mapping = dict(zip(df["chembl"], df["schembl_id"]))
    assert mapping["CHEMBL941"] == "SCHEMBL3827"
    assert mapping["CHEMBL_OTHER"] == "SCHEMBL9"
    assert pd.isna(mapping["CHEMBL_BIOLOGIC"])
    assert (workdir / "run_chemicals.tsv").exists()

    # Second run: everything is cached, so ChEMBL is not asked again.
    harmonizer.harmonize_chemicals("run", source=str(bulk_dir))
    assert len(chembl_calls) == 1


def test_run_without_any_chemicals(workdir, bulk_dir, chembl_calls):
    """Genes without ChEMBL chemicals: no lookups, an empty chemicals file, no patents."""
    (workdir / "empty_gene_to_chemicals.json").write_text(json.dumps({"GENE1": [], "GENE2": []}))
    assert harmonizer.harmonize_chemicals("empty", source=str(bulk_dir)).empty
    assert chembl_calls == []
    assert (workdir / "empty_chemicals.tsv").exists()
    assert enrichment.extract_patent("empty", source=str(bulk_dir)).empty


def test_extract_patent_without_chemical_file(workdir):
    assert list(enrichment.extract_patent("never_harmonized").columns) == PATENT_COLUMNS


def test_harmonize_requires_gene_file(workdir, bulk_dir):
    with pytest.raises(FileNotFoundError, match="experimental data extractor"):
        harmonizer.harmonize_chemicals("missing", source=str(bulk_dir))


def test_harmonize_from_user_chemicals(workdir, bulk_dir):
    pd.DataFrame({"chembl": ["CHEMBL941", "CHEMBL_NOTINSCHEMBL"]}).to_csv(
        workdir / "user_chemicals.tsv", sep="\t", index=False
    )
    df = harmonizer.harmonize_chemicals("user", from_genes=False, source=str(bulk_dir))
    assert df.set_index("chembl")["schembl_id"].to_dict().get("CHEMBL941") == "SCHEMBL3827"
    assert pd.isna(df.set_index("chembl").loc["CHEMBL_NOTINSCHEMBL", "schembl_id"])


def test_extract_patent(workdir, bulk_dir):
    write_gene_file(workdir)
    harmonizer.harmonize_chemicals("run", source=str(bulk_dir))
    df = enrichment.extract_patent("run", patent_year=2000, source=str(bulk_dir))

    assert list(df.columns) == PATENT_COLUMNS
    rows = df[["chembl", "surechembl", "patent_id", "date"]].values.tolist()
    assert rows == [
        ["CHEMBL941", "SCHEMBL3827", "WO-2011041462-A2", "2011-04-07"],
        ["CHEMBL941", "SCHEMBL3827", "EP-2389136-A1", "2011-11-30"],
        ["CHEMBL_OTHER", "SCHEMBL9", "EP-2389136-A1", "2011-11-30"],
    ]
    first = df.iloc[0]
    assert first["ipc"] == "A61P 35/00"
    assert first["assignee"] == "PHARMA A"
    assert first["sections"] == "description; claims"

    saved = pd.read_csv(workdir / "run_patent_data.tsv", sep="\t", dtype=str)
    assert saved["patent_id"].tolist() == df["patent_id"].tolist()
    state = json.loads((workdir / "run_patent_data.json").read_text())
    assert state["compounds_done"] == ["SCHEMBL3827", "SCHEMBL9"]
    assert state["settings"]["year"] == 2000


class ExplodingBulk:
    """Stands in for SureChEMBLBulk when a test expects no query to happen."""

    def __init__(self, source):
        self.source = source

    def patents_for_compounds(self, *args, **kwargs):
        raise AssertionError("SureChEMBL should not have been queried")


def test_extract_patent_uses_cache(workdir, bulk_dir):
    write_gene_file(workdir)
    harmonizer.harmonize_chemicals("run", source=str(bulk_dir))
    first = enrichment.extract_patent("run", source=str(bulk_dir))

    # Same settings: served from the cache without touching SureChEMBL.
    with SureChEMBLBulk(str(bulk_dir)) as real:
        cached = enrichment.extract_patent("run", bulk=ExplodingBulk(real.source))
    pd.testing.assert_frame_equal(first, cached, check_dtype=False)

    # Different settings: queried again, with the claims-only filter applied.
    claims = enrichment.extract_patent("run", source=str(bulk_dir), sections=["claims"])
    assert claims[["chembl", "patent_id"]].values.tolist() == [
        ["CHEMBL941", "WO-2011041462-A2"],
        ["CHEMBL_OTHER", "EP-2389136-A1"],
    ]


def test_one_chemical_with_two_surechembl_records(workdir, bulk_dir):
    """A ChEMBL chemical matching two SureChEMBL records gets one row per patent."""
    pd.DataFrame(
        {"chembl": ["CHEMBL941", "CHEMBL941"], "schembl_id": ["SCHEMBL3827", "SCHEMBL9"]}
    ).to_csv(workdir / "dup_chemicals.tsv", sep="\t", index=False)

    df = enrichment.extract_patent("dup", source=str(bulk_dir))
    assert df["patent_id"].tolist() == ["WO-2011041462-A2", "EP-2389136-A1"]
    ep = df.set_index("patent_id").loc["EP-2389136-A1"]
    assert ep["surechembl"] == "SCHEMBL3827; SCHEMBL9"
    assert ep["sections"] == "description; claims"  # 3827 in description, 9 in claims


def test_extract_patent_no_mapped_chemicals(workdir):
    pd.DataFrame({"chembl": ["CHEMBL_X"], "schembl_id": [None]}).to_csv(
        workdir / "none_chemicals.tsv", sep="\t", index=False
    )
    assert list(enrichment.extract_patent("none").columns) == PATENT_COLUMNS


def test_cli_run_patent_extractor(workdir, bulk_dir):
    chemical_file = workdir / "input.tsv"
    pd.DataFrame({"chembl": ["CHEMBL941"]}).to_csv(chemical_file, sep="\t", index=False)

    result = CliRunner().invoke(
        pemt.cli.main,
        [
            "run-patent-extractor",
            "--name", "cli",
            "--chemical",
            "--chemical-data", str(chemical_file),
            "--surechembl-source", str(bulk_dir),
            "--sections", "claims",
        ],
    )
    assert result.exit_code == 0, result.output
    out = pd.read_csv(workdir / "cli_patent_data.tsv", sep="\t", dtype=str)
    assert out["patent_id"].tolist() == ["WO-2011041462-A2"]
    assert not (workdir / "cleaned_cli_patent_data.tsv").exists()


def test_cli_help_has_no_selenium_options():
    result = CliRunner().invoke(pemt.cli.main, ["run-pemt", "--help"])
    assert result.exit_code == 0
    assert "--surechembl-source" in result.output
    assert "chromedriver" not in result.output


@pytest.mark.live
def test_live_chembl_structures():
    """Asks the real ChEMBL web service for imatinib's name and InChIKey."""
    from pemt.utils import get_chembl_structures

    df = get_chembl_structures(["CHEMBL941", "CHEMBL_DOES_NOT_EXIST"]).set_index("chembl")
    assert df.loc["CHEMBL941", "inchi_key"] == "KTUFNOKKBVMGRW-UHFFFAOYSA-N"
    assert df.loc["CHEMBL941", "name"] == "IMATINIB"
    assert pd.isna(df.loc["CHEMBL_DOES_NOT_EXIST", "inchi_key"])
