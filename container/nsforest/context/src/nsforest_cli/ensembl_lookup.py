# container/nsforest/context/src/nsforest_cli/ensembl_lookup.py
from __future__ import annotations

import time
from typing import Dict, Iterable, Optional
import requests

ENSEMBL_REST = "https://rest.ensembl.org"


def ensembl_lookup (
    ensg_ids: Iterable[str],
    *,
    sleep_sec: float = 0.05,
    timeout_sec: float = 15.0,
    prefer_dbname: str = "EntrezGene",
) -> Dict[str, str]:
    """
    Map ENSG ids -> gene symbols using Ensembl xrefs endpoint.
    No caching; each id is fetched on demand and falls back to the ENSG id if not found.

    Parameters
    ----------
    ensg_ids : iterable of str
        ENSG identifiers to resolve.
    sleep_sec : float
        Small delay between calls to be polite to the API.
    timeout_sec : float
        HTTP timeout per request.
    prefer_dbname : str
        Prefer entries with this 'dbname' for display_id (default: EntrezGene).

    Returns
    -------
    dict
        { ENSG -> symbol_or_ENSG }
    """
    headers = {"Content-Type": "application/json"}
    mapping: Dict[str, str] = {}

    for eid in [e for e in ensg_ids if e]:
        url = f"{ENSEMBL_REST}/xrefs/id/{eid}?content-type=application/json"
        symbol: Optional[str] = None
        try:
            r = requests.get(url, headers=headers, timeout=timeout_sec)
            r.raise_for_status()
            data = r.json() or []
            # Prefer a specific dbname first
            for entry in data:
                if entry.get("dbname") == prefer_dbname and entry.get("display_id"):
                    symbol = entry["display_id"]
                    break
            # Fallback to any display_id
            if symbol is None:
                for entry in data:
                    if entry.get("display_id"):
                        symbol = entry["display_id"]
                        break
        except Exception:
            symbol = None

        mapping[eid] = symbol or eid
        if sleep_sec:
            time.sleep(sleep_sec)

    return mapping

