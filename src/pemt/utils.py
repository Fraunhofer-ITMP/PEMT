# -*- coding: utf-8 -*-

import logging
import os
from typing import Dict, Iterable, List, Optional

import pandas as pd

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)

"""Protein mapper functions"""

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
CHEMBL_STATUS_URL = "https://www.ebi.ac.uk/chembl/api/data/status.json"


def load_hgnc(
    cache_dir: str, max_age_days: int = 30, url: str = HGNC_URL
) -> pd.DataFrame:
    """Load the HGNC complete set, downloading it when the local copy is missing or old.

    HGNC publishes the file monthly. If a download fails but an older copy exists, the
    older copy is used with a warning, so runs keep working offline.

    :param cache_dir: Folder for ``hgnc_complete_set.txt``
    :param max_age_days: Re-download when the cached file is older than this
    :param url: Download location of the HGNC complete set (TSV)
    :returns: Columns ``hgnc_id``, ``symbol``, ``status``, ``prev_symbol``, ``alias_symbol``
        and ``uniprot_ids`` (multi-valued fields are ``|``-separated)
    """
    import time

    import requests

    path = os.path.join(cache_dir, "hgnc_complete_set.txt")
    age_days = (
        (time.time() - os.path.getmtime(path)) / 86400 if os.path.exists(path) else None
    )

    if age_days is None or age_days > max_age_days:
        try:
            logger.info(f"Downloading the HGNC complete set from {url}")
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            os.makedirs(cache_dir, exist_ok=True)
            with open(path + ".tmp", "wb") as f:
                f.write(response.content)
            os.replace(path + ".tmp", path)
        except Exception as exc:  # network problems: fall back to an older copy
            if age_days is None:
                raise RuntimeError(
                    f"Could not download the HGNC complete set from {url}: {exc}"
                ) from exc
            logger.warning(
                f"Could not update the HGNC file ({exc}); using the {age_days:.0f} day old copy"
            )

    columns = [
        "hgnc_id",
        "symbol",
        "status",
        "prev_symbol",
        "alias_symbol",
        "uniprot_ids",
    ]
    return pd.read_csv(
        path, sep="\t", dtype=str, usecols=columns, keep_default_na=False
    )


def symbols_to_uniprot(
    symbols: Iterable[str], hgnc: pd.DataFrame
) -> Dict[str, List[str]]:
    """Map gene symbols to UniProt accessions with the HGNC complete set.

    Each symbol is resolved in this order: approved symbol, previous symbol, alias. A
    previous symbol or alias is only used if it points to exactly one approved gene.
    Matching ignores case and surrounding whitespace.

    :param symbols: Gene symbols, e.g. ``["ABL1", "MLL"]``
    :param hgnc: Output of :func:`load_hgnc`
    :returns: ``{symbol: [UniProt accessions]}``; the list is empty for symbols that could
        not be resolved or whose gene has no protein (e.g. non-coding RNA genes)
    """
    approved = hgnc[hgnc["status"] == "Approved"]
    uniprot_of = {
        row.symbol.upper(): [u for u in row.uniprot_ids.split("|") if u]
        for row in approved.itertuples()
    }

    def index(column: str) -> Dict[str, set]:
        lookup: Dict[str, set] = {}
        for row in approved[approved[column] != ""].itertuples():
            for name in getattr(row, column).split("|"):
                lookup.setdefault(name.strip().upper(), set()).add(row.symbol.upper())
        return lookup

    previous, aliases = index("prev_symbol"), index("alias_symbol")

    result: Dict[str, List[str]] = {}
    renamed, unresolved, ambiguous = [], [], []
    for symbol in symbols:
        key = str(symbol).strip().upper()
        if key in uniprot_of:
            result[symbol] = uniprot_of[key]
            continue
        for lookup in (previous, aliases):
            candidates = lookup.get(key, set())
            if len(candidates) == 1:
                current = next(iter(candidates))
                result[symbol] = uniprot_of[current]
                renamed.append(f"{symbol}->{current}")
                break
            if len(candidates) > 1:
                ambiguous.append(f"{symbol} ({', '.join(sorted(candidates))})")
                result[symbol] = []
                break
        else:
            unresolved.append(symbol)
            result[symbol] = []

    if renamed:
        logger.info(
            f"{len(renamed)} symbols mapped via previous/alias symbols: {', '.join(renamed[:10])}"
        )
    if ambiguous:
        logger.warning(
            f"{len(ambiguous)} symbols are ambiguous and were skipped: {'; '.join(ambiguous[:10])}"
        )
    if unresolved:
        logger.warning(
            f"{len(unresolved)} symbols not found in HGNC: {', '.join(map(str, unresolved[:10]))}"
        )
    return result


def get_chembl_release() -> Optional[str]:
    """Return the ChEMBL release served by the web services (e.g. ``"ChEMBL_36"``)."""
    import requests

    try:
        response = requests.get(CHEMBL_STATUS_URL, timeout=30)
        response.raise_for_status()
        return response.json().get("chembl_db_version")
    except Exception as exc:
        logger.warning(f"Could not read the ChEMBL release: {exc}")
        return None


def get_single_protein_targets(uniprot_id: Optional[str]) -> List[str]:
    """Find the ChEMBL SINGLE PROTEIN target(s) for a UniProt accession.

    Many proteins also appear in ChEMBL as part of complexes, fusion proteins, protein
    families or degrader (protein-protein interaction) targets; only the target that is
    the protein itself is returned, normally exactly one.

    :param uniprot_id: UniProt accession, e.g. ``P00519`` for ABL1
    :returns: ChEMBL target ids (empty if the protein is not a ChEMBL target)
    """
    if not isinstance(uniprot_id, str) or not uniprot_id:
        return []
    from chembl_webresource_client.new_client import new_client

    targets = new_client.target.filter(
        target_components__accession=uniprot_id, target_type="SINGLE PROTEIN"
    ).only(["target_chembl_id"])
    return sorted({t["target_chembl_id"] for t in targets})


"""Chemical mapper functions"""


def get_chembl_structures(
    chembl_ids: Iterable[str], chunk_size: int = 50
) -> pd.DataFrame:
    """Get the preferred name and standard InChIKey of ChEMBL compounds.

    :param chembl_ids: ChEMBL compound identifiers (e.g. ``CHEMBL941``)
    :param chunk_size: Number of compounds requested per ChEMBL web service call
    :returns: Columns ``chembl``, ``name`` and ``inchi_key``. Compounds without a
        structure (e.g. biologics) have no InChIKey; the name falls back to the id.
    """
    from chembl_webresource_client.new_client import new_client

    ids = sorted({i for i in chembl_ids if isinstance(i, str) and i})
    rows = []
    for start in range(0, len(ids), chunk_size):
        chunk = ids[start : start + chunk_size]
        found = new_client.molecule.filter(molecule_chembl_id__in=chunk).only(
            ["molecule_chembl_id", "pref_name", "molecule_structures"]
        )
        for molecule in found:
            structures = molecule.get("molecule_structures") or {}
            rows.append(
                {
                    "chembl": molecule["molecule_chembl_id"],
                    "name": molecule.get("pref_name") or molecule["molecule_chembl_id"],
                    "inchi_key": structures.get("standard_inchi_key"),
                }
            )

    df = pd.DataFrame(rows, columns=["chembl", "name", "inchi_key"])
    missing = set(ids) - set(df["chembl"])
    if missing:
        logger.warning(
            f"{len(missing)} ChEMBL ids not found in ChEMBL, e.g. {sorted(missing)[:5]}"
        )
        df = pd.concat(
            [
                df,
                pd.DataFrame(
                    {
                        "chembl": sorted(missing),
                        "name": sorted(missing),
                        "inchi_key": None,
                    }
                ),
            ],
            ignore_index=True,
        )
    return df.sort_values("chembl", ignore_index=True)
