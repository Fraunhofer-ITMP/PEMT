# -*- coding: utf-8 -*-

"""Query the SureChEMBL bulk data (Apache Parquet) with DuckDB.

SureChEMBL publishes its full database every two weeks as Parquet files at
https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/. DuckDB can query
those files straight from the EBI server, reading only the columns a query needs, or
from a local copy of them. Both are handled the same way here: a *source* is either a
release name (``"latest"`` or a dated folder such as ``"2026-09-22"``), a URL, or a
local folder containing ``compounds.parquet``, ``patents.parquet`` and
``patent_compound_map.parquet``.

The compound-to-patent map (about 1.5 billion rows) is sorted by patent, so any lookup
by compound scans its whole ``compound_id`` column once. That costs the same for one
compound as for a thousand, so callers should pass all compounds of a run together.
"""

import logging
import os
import re
from typing import Dict, Iterable, List, Optional

import duckdb
import pandas as pd

logger = logging.getLogger(__name__)

BULK_DATA_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/"

#: Patent sections a compound can be found in (``fields.parquet`` of the bulk data).
SECTIONS: Dict[str, int] = {
    "description": 1,
    "claims": 2,
    "abstract": 3,
    "title": 4,
    "image": 5,
    "mol_attachment": 6,
}
SECTION_NAMES = {v: k for k, v in SECTIONS.items()}

TABLES = ("compounds", "patents", "patent_compound_map")

PATENT_COLUMNS = [
    "compound_id",
    "patent_number",
    "publication_date",
    "family_id",
    "ipc",
    "assignee",
    "sections",
]

_RELEASE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
_IPC_PREFIX_RE = re.compile(r"^[A-Z0-9/]+$")


def resolve_source(source: Optional[str] = None) -> str:
    """Turn a release name, URL or local folder into a base location ending in ``/``.

    :param source: ``None`` or ``"latest"`` for the newest release on the EBI server, a
        release date such as ``"2026-09-22"``, a URL, or a path to a local folder.
    """
    if source is None or source == "latest":
        return BULK_DATA_URL + "latest/"
    if _RELEASE_RE.match(source):
        return f"{BULK_DATA_URL}{source}/"
    if source.startswith(("http://", "https://", "s3://")):
        return source if source.endswith("/") else source + "/"

    path = os.path.abspath(os.path.expanduser(source))
    missing = [t for t in TABLES if not os.path.exists(os.path.join(path, f"{t}.parquet"))]
    if missing:
        raise FileNotFoundError(
            f"SureChEMBL source folder {path} is missing: "
            + ", ".join(f"{t}.parquet" for t in missing)
        )
    return path + os.sep


def _ipc_condition(variable: str, prefixes: Iterable[str]) -> str:
    """SQL condition that is true when an IPC code starts with one of ``prefixes``.

    Codes in the bulk data look like ``"A61K 31/42"``; spaces are ignored so a prefix
    such as ``"A61K31"`` also works. Prefixes are validated before being inlined.
    """
    by_length: Dict[int, List[str]] = {}
    for prefix in prefixes:
        prefix = prefix.replace(" ", "").upper()
        if not _IPC_PREFIX_RE.match(prefix):
            raise ValueError(f"Invalid IPC prefix: {prefix!r}")
        by_length.setdefault(len(prefix), []).append(prefix)
    if not by_length:
        raise ValueError("ipc_prefixes must not be empty; pass None to disable the filter")
    parts = [
        f"left(replace({variable}, ' ', ''), {length}) IN ({', '.join(repr(p) for p in sorted(codes))})"
        for length, codes in sorted(by_length.items())
    ]
    return "(" + " OR ".join(parts) + ")"


class SureChEMBLBulk:
    """A DuckDB connection to one SureChEMBL bulk data release.

    Example::

        bulk = SureChEMBLBulk()                      # latest release, read remotely
        ids = bulk.compounds_for_inchikeys(["KTUFNOKKBVMGRW-UHFFFAOYSA-N"])
        patents = bulk.patents_for_compounds(ids["surechembl_id"], ipc_prefixes=["A61P"])
    """

    def __init__(self, source: Optional[str] = None, progress_bar: bool = True):
        """Open a connection.

        :param source: See :func:`resolve_source`.
        :param progress_bar: Show DuckDB's progress bar for long scans.
        """
        self.source = resolve_source(source)
        self.con = duckdb.connect()
        if self.source.startswith(("http://", "https://", "s3://")):
            self.con.execute("INSTALL httpfs; LOAD httpfs;")
        self.con.execute(f"SET enable_progress_bar = {str(progress_bar).lower()};")
        for table in TABLES:
            location = f"{self.source}{table}.parquet".replace("'", "''")
            self.con.execute(f"CREATE VIEW {table} AS SELECT * FROM read_parquet('{location}')")
        logger.info("Using SureChEMBL bulk data from %s", self.source)

    def close(self) -> None:
        """Close the DuckDB connection."""
        self.con.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def compounds_for_inchikeys(self, inchikeys: Iterable[str]) -> pd.DataFrame:
        """Map standard InChIKeys to SureChEMBL compound ids.

        :param inchikeys: InChIKeys, e.g. from ChEMBL's ``molecule_structures``.
        :returns: Columns ``inchi_key`` and ``surechembl_id``; keys without a match are
            left out, and a key can in rare cases map to more than one compound.
        """
        keys = sorted({k for k in inchikeys if isinstance(k, str) and k})
        if not keys:
            return pd.DataFrame(columns=["inchi_key", "surechembl_id"])
        return self.con.execute(
            """
            SELECT c.inchi_key, c.id AS surechembl_id
            FROM compounds c
            WHERE c.inchi_key IN (SELECT unnest(?::VARCHAR[]))
            ORDER BY c.inchi_key, c.id
            """,
            [keys],
        ).fetchdf()

    def patents_for_compounds(
        self,
        compound_ids: Iterable[int],
        ipc_prefixes: Optional[Iterable[str]] = None,
        min_year: Optional[int] = None,
        sections: Optional[Iterable[str]] = None,
    ) -> pd.DataFrame:
        """Find the patents that mention the given compounds, with PEMT's filters applied.

        All compounds are looked up in a single scan, so pass every compound of a run
        in one call.

        :param compound_ids: SureChEMBL compound ids (integers, without ``SCHEMBL``).
        :param ipc_prefixes: Keep patents with at least one IPCR code starting with one
            of these (e.g. ``pemt.constants.VALID_CODES``). ``None`` keeps all patents.
        :param min_year: Keep patents published in or after this year.
        :param sections: Only count a compound if it appears in these sections (keys of
            :data:`SECTIONS`, e.g. ``["claims"]``). ``None`` means any section.
        :returns: One row per compound and patent with the columns in
            :data:`PATENT_COLUMNS`. ``ipc`` holds the matching IPCR codes (all codes if
            no IPC filter is set) joined by ``"; "``; ``assignee`` is the first listed
            assignee; ``sections`` lists where the compound appears in that patent.
        """
        ids = sorted({int(str(i).upper().replace("SCHEMBL", "")) for i in compound_ids})
        if not ids:
            return pd.DataFrame(columns=PATENT_COLUMNS)

        map_filters = ["m.compound_id IN (SELECT unnest(?::BIGINT[]))"]
        if sections is not None:
            unknown = set(sections) - set(SECTIONS)
            if unknown:
                raise ValueError(f"Unknown sections {sorted(unknown)}; choose from {sorted(SECTIONS)}")
            field_ids = sorted(SECTIONS[s] for s in set(sections))
            if not field_ids:
                raise ValueError("sections must not be empty; pass None for any section")
            map_filters.append(f"m.field_id IN ({', '.join(map(str, field_ids))})")

        patent_filters = []
        if min_year is not None:
            patent_filters.append(f"p.publication_date >= DATE '{int(min_year):04d}-01-01'")
        if ipc_prefixes is not None:
            matched = f"list_filter(p.ipcr, x -> {_ipc_condition('x', ipc_prefixes)})"
            patent_filters.append(f"len({matched}) > 0")
        else:
            matched = "p.ipcr"

        query = f"""
            WITH hits AS (
                SELECT m.compound_id, m.patent_id,
                       list(DISTINCT m.field_id ORDER BY m.field_id) AS field_ids
                FROM patent_compound_map m
                WHERE {' AND '.join(map_filters)}
                GROUP BY m.compound_id, m.patent_id
            )
            SELECT h.compound_id,
                   p.patent_number,
                   p.publication_date,
                   p.family_id,
                   array_to_string({matched}, '; ') AS ipc,
                   p.assignee[1] AS assignee,
                   h.field_ids
            FROM hits h
            JOIN patents p ON p.id = h.patent_id
            {'WHERE ' + ' AND '.join(patent_filters) if patent_filters else ''}
            ORDER BY h.compound_id, p.publication_date, p.patent_number
        """
        df = self.con.execute(query, [ids]).fetchdf()
        df["sections"] = df.pop("field_ids").map(
            lambda f: "; ".join(SECTION_NAMES.get(int(i), str(i)) for i in f)
        )
        return df[PATENT_COLUMNS]
