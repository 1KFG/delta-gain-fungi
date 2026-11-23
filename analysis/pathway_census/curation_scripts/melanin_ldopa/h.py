"""Curation helpers. All requests go through kbsources (rate limit, user agent, cache)."""
import sys, json
from pathlib import Path
PKG = Path(__file__).resolve().parents[2]  # analysis/pathway_census
sys.path.insert(0, str(PKG))
from pathway_census import kbsources as k
from pathway_census.kbsources import uniprot_entry, pubmed_abstract, chebi_compound, pubchem_by_inchikey, KbFetchError

def raw(namespace, key, url, as_json=True):
    """Fetch a URL not covered by a kbsources fetcher, through its rate-limited fetch and cache."""
    def parse(r):
        if not r.strip():
            raise KbFetchError(f"{namespace} {key}: empty response")
        return json.loads(r) if as_json else r
    return k._cached(namespace, key, url, parse, None, None)
