#!/usr/bin/env python3
"""Independent verification of the curated knowledge base against public records.

Usage: /usr/bin/python3.12 bin/verify_kb.py [pathway ...] [--report FILE] [--strict]
Default: all four pathways (carotenoid, melanin_dhn, melanin_ldopa, mycosporine). The report goes
to kb/VERIFICATION_REPORT.md. Requests go through pathway_census.kbsources (cache, rate limit).

``verify_pathway(yaml_path, kb_dir, fetch=None)`` returns a DataFrame with columns
``check, subject, expected, found, status``. Status is ``pass``, ``fail`` or ``unchecked``.
``unchecked`` always carries its reason in ``found``. This script reads the written files (YAML,
FASTA, TSV) and the public records. It does not use the data of the curation scripts.

What a status means:
- fail: the public record was read and it disagrees with the file, or the record does not exist
  (HTTP 404/400, or an empty PubMed result).
- unchecked: the record could not be read (network error, no abstract, no PubChem compound, no
  accession in the YAML) or a quote check that is a review flag for tiers T1b and T1c.
- pass: the record was read and agrees.
Gene names: UniProt gene names, synonyms, ORF names and ordered locus names are read from the full
UniProt record (``kbsources.uniprot_entry`` does not read ordered locus names). The kind that
matched is written in ``found``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pathway_census import kbsources as ks  # noqa: E402
from pathway_census.defs import load_pathway  # noqa: E402
from pathway_census.kb import load_kb, species_key  # noqa: E402
from pathway_census.kbsources import KbFetchError  # noqa: E402

PKG = Path(__file__).resolve().parent.parent
PATHWAYS = ("carotenoid", "melanin_dhn", "melanin_ldopa", "mycosporine")
COLS = ["check", "subject", "expected", "found", "status"]
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
GENBANK_RE = re.compile(r"^[A-Z]{3}\d{5,7}(\.\d+)?$")
NCBI_PROTEIN_URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={acc}"
                    "&rettype=gp&retmode=text")
ABBREV = {"f", "sp", "spp", "subsp", "var", "fig", "figs", "e.g", "i.e", "vs", "al", "ca", "cf",
          "approx", "no", "nos", "st", "dr", "eq", "ref", "refs", "cv", "ssp", "pv", "ex", "syn"}


def norm_text(s: str) -> str:
    return " ".join(str(s).split())


def norm_seq(s: str) -> str:
    return re.sub(r"[\s*]", "", str(s)).upper()


def sha(s: str) -> str:
    return hashlib.sha256(norm_seq(s).encode()).hexdigest()


def is_one_sentence(text: str) -> bool:
    """True when ``text`` has no sentence break inside. A break is . ! or ? followed by space and
    a capital letter, digit or opening bracket, unless the word before the full stop is a single
    letter (genus initial), a known abbreviation, or the full stop sits in a number."""
    t = norm_text(text)
    for m in re.finditer(r"([.!?])\s+(?=[A-Z0-9(\[])", t):
        before = t[:m.start()]
        word = re.search(r"([A-Za-z.]+)$", before)
        w = word.group(1).lower() if word else ""
        if m.group(1) == "." and (len(w) == 1 or w in ABBREV or w.rstrip(".") in ABBREV):
            continue
        if m.end() >= len(t):
            continue
        return False
    return True


def read_fasta(path: Path) -> dict[str, str] | None:
    """Map each header token (split on whitespace and |) to the sequence. None if no file."""
    if not Path(path).is_file():
        return None
    out: dict[str, str] = {}
    toks: list[str] = []
    seq: list[str] = []

    def flush():
        for t in toks:
            out.setdefault(t, "".join(seq))

    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            flush()
            toks = [t for t in line[1:].replace("|", " ").split() if t]
            seq = []
        else:
            seq.append(line)
    flush()
    return out


def not_found(exc: KbFetchError) -> bool:
    m = str(exc)
    return ("HTTP Error 404" in m or "HTTP Error 400" in m or "no PubMed article" in m
            or "record version is" in m)


class Verifier:
    def __init__(self, fetch=None, cache_dir=None):
        self.fetch, self.cache_dir = fetch, cache_dir
        self.rows: list[dict] = []
        self._rec: dict = {}
        self._pm: dict = {}

    def add(self, check, subject, expected, found, status):
        assert status in ("pass", "fail", "unchecked")
        self.rows.append({"check": check, "subject": str(subject), "expected": str(expected),
                          "found": str(found), "status": status})

    # records ---------------------------------------------------------------------------------
    def record(self, acc: str) -> dict:
        """dict: kind, organism, sequence, length, names [(kind, name)], protein_names, citations.
        Raises KbFetchError."""
        if acc in self._rec:
            r = self._rec[acc]
            if isinstance(r, Exception):
                raise r
            return r
        try:
            r = self._load_record(acc)
        except KbFetchError as e:
            self._rec[acc] = e
            raise
        self._rec[acc] = r
        return r

    def _load_record(self, acc: str) -> dict:
        if UNIPROT_RE.match(acc):
            d = ks._cached("uniprot_full", acc, ks.UNIPROT_URL.format(acc=acc), json.loads,
                           self.fetch, self.cache_dir)
            if d.get("primaryAccession") != acc:
                raise KbFetchError(f"UniProt {acc}: record has primaryAccession {d.get('primaryAccession')}")
            names: list[tuple[str, str]] = []
            for g in d.get("genes", []):
                if "geneName" in g:
                    names.append(("gene name", g["geneName"]["value"]))
                names += [("gene name (synonym)", s["value"]) for s in g.get("synonyms", [])]
                names += [("ORF name", o["value"]) for o in g.get("orfNames", [])]
                names += [("ordered locus name", o["value"]) for o in g.get("orderedLocusNames", [])]
            pn: list[str] = []
            pdesc = d.get("proteinDescription", {})
            for n in [pdesc.get("recommendedName")] + pdesc.get("alternativeNames", []) + pdesc.get("submissionNames", []):
                if n and n.get("fullName", {}).get("value"):
                    pn.append(n["fullName"]["value"])
                for s in (n or {}).get("shortNames", []):
                    pn.append(s["value"])
            seq = d["sequence"]["value"]
            if len(seq) != d["sequence"]["length"]:
                raise KbFetchError(f"UniProt {acc}: sequence length does not match length field")
            return {"kind": "UniProt", "organism": d["organism"]["scientificName"], "sequence": seq,
                    "length": len(seq), "names": names, "protein_names": pn, "citations": [],
                    "reviewed": d.get("entryType", "").startswith("UniProtKB reviewed")}
        if GENBANK_RE.match(acc):
            base = acc.split(".")[0]

            def parse(raw):
                if not raw.strip():
                    raise KbFetchError(f"NCBI protein {acc}: empty response")
                return raw

            t = ks._cached("ncbi_protein", base + "_gp", NCBI_PROTEIN_URL.format(acc=base), parse,
                           self.fetch, self.cache_dir)
            ver = re.search(r"^VERSION\s+(\S+)", t, flags=re.M)
            if not ver:
                raise KbFetchError(f"NCBI protein {acc}: unusable response (no VERSION line)")
            if "." in acc and ver.group(1) != acc:
                raise KbFetchError(f"NCBI protein {acc}: record version is {ver.group(1)}")

            def q(key):
                m = re.search(r'/%s="([^"]*)"' % key, t)
                return " ".join(m.group(1).split()) if m else None

            seq = re.sub(r"[^a-z]", "", t.split("ORIGIN")[1].split("//")[0]).upper()
            length = int(re.search(r"^LOCUS\s+\S+\s+(\d+) aa", t, flags=re.M).group(1))
            if len(seq) != length:
                raise KbFetchError(f"NCBI protein {acc}: sequence length {len(seq)} != LOCUS length {length}")
            cites = [tuple(" ".join((x or "").split()) for x in r)
                     for r in re.findall(r"TITLE\s+(.*?)\n\s+JOURNAL\s+(.*?)\n(?:\s+REMARK\s+(.*?)\n)?", t, flags=re.S)]
            return {"kind": "GenBank", "organism": q("organism"), "sequence": seq, "length": length,
                    "names": [("gene name", q("gene"))] if q("gene") else [],
                    "protein_names": [q("product")] if q("product") else [], "citations": cites,
                    "version": ver.group(1)}
        raise KbFetchError(f"accession {acc}: not a UniProt or GenBank protein accession pattern")

    def pubmed(self, pmid: str) -> dict:
        if pmid not in self._pm:
            try:
                self._pm[pmid] = ks.pubmed_abstract(pmid, fetch=self.fetch, cache_dir=self.cache_dir)
            except KbFetchError as e:
                self._pm[pmid] = e
        r = self._pm[pmid]
        if isinstance(r, Exception):
            raise r
        return r

    # checks ----------------------------------------------------------------------------------
    def check_pmid_exists(self, check, subject, pmid):
        try:
            a = self.pubmed(pmid)
            self.add(check, subject, f"PMID {pmid} in PubMed", f"title: {a['title'][:80]}", "pass")
            return a
        except KbFetchError as e:
            self.add(check, subject, f"PMID {pmid} in PubMed", str(e)[:200],
                     "fail" if not_found(e) else "unchecked")
            return None

    def match_name(self, rec, gene):
        """(kind, name) of the first record name equal to ``gene``, else None."""
        for kind, n in rec["names"]:
            if n == gene:
                return kind, n
        return None

    def check_record(self, prefix, subject, acc, species, gene, fasta, region=None, step_fa="", prov_region=None):
        """exists, organism, gene, sequence checks for one accession. Returns the record or None."""
        try:
            rec = self.record(acc)
        except KbFetchError as e:
            self.add(f"{prefix}_exists", subject, f"record {acc}", str(e)[:200],
                     "fail" if not_found(e) else "unchecked")
            return None
        self.add(f"{prefix}_exists", subject, f"record {acc}",
                 f"{rec['kind']} record, organism {rec['organism']}, {rec['length']} aa", "pass")
        got, want = species_key(rec["organism"]), species_key(species)
        self.add(f"{prefix}_organism", subject, species, rec["organism"], "pass" if got == want else "fail")
        if gene:
            hit = self.match_name(rec, gene)
            if hit:
                self.add(f"{prefix}_gene", subject, gene, f"{hit[0]} '{hit[1]}'", "pass")
            else:
                self.add(f"{prefix}_gene", subject, gene,
                         "no match; names in record: " + ", ".join(n for _, n in rec["names"]), "fail")
        if fasta is None:
            self.add(f"{prefix}_sequence", subject, f"FASTA record {acc}", f"no FASTA file {step_fa}", "unchecked")
        elif acc not in fasta:
            self.add(f"{prefix}_sequence", subject, f"FASTA record {acc}", f"accession not in {step_fa}", "fail")
        else:
            fs, rs = norm_seq(fasta[acc]), norm_seq(rec["sequence"])
            if fs == rs:
                self.add(f"{prefix}_sequence", subject, f"sha256 {sha(rs)[:12]}",
                         f"equal, {len(rs)} aa, sha256 {sha(fs)[:12]}", "pass")
            elif region and prov_region is None:
                self.add(f"{prefix}_sequence", subject, "source coordinates in provenance.tsv",
                         f"region in YAML but not in provenance (YAML region {list(region)})", "fail")
            elif region:
                a, b = prov_region
                sl = rs[a - 1:b]
                if fs == sl:
                    self.add(f"{prefix}_sequence", subject, f"source record positions {a}-{b}",
                             f"FASTA equals the fetched sequence at the provenance region {a}-{b} ({len(fs)} of {len(rs)} aa)", "pass")
                else:
                    where = f"; FASTA is a substring at {rs.index(fs) + 1}" if fs in rs else ""
                    self.add(f"{prefix}_sequence", subject, f"source record positions {a}-{b} (sha256 {sha(sl)[:12]})",
                             f"FASTA sha256 {sha(fs)[:12]} ({len(fs)} aa) differs from the provenance-region slice{where}", "fail")
            else:
                self.add(f"{prefix}_sequence", subject, f"sha256 {sha(rs)[:12]} ({len(rs)} aa)",
                         f"FASTA sha256 {sha(fs)[:12]} ({len(fs)} aa)", "fail")
        if region and fasta is not None and acc in fasta:
            n = len(norm_seq(fasta[acc]))
            self.add(f"{prefix}_region", subject, f"YAML region [1, {n}] (FASTA coordinates)", f"YAML region {list(region)}",
                     "pass" if list(region) == [1, n] else "fail")
        return rec

    def check_genbank_cites(self, subject, rec, pmid):
        if not pmid:
            return
        try:
            a = self.pubmed(pmid)
        except KbFetchError as e:
            self.add("ref_record_cites_pmid", subject, f"GenBank record cites PMID {pmid}",
                     f"PubMed record not read: {str(e)[:150]}", "unchecked")
            return
        tt = a["title"].rstrip(".").casefold()
        for title, journal, remark in rec["citations"]:
            if title.rstrip(".").casefold() == tt:
                self.add("ref_record_cites_pmid", subject, f"GenBank record cites PMID {pmid}",
                         f"GenBank REFERENCE title equals PubMed title ({journal}; {remark})", "pass")
                return
        self.add("ref_record_cites_pmid", subject, f"GenBank record cites PMID {pmid}",
                 "no GenBank REFERENCE title equals the PubMed title '%s'" % a["title"], "fail")

    def names_quote(self, rec, gene, quote):
        """First candidate name found in the quote as a whole token. Names of 3 or more characters
        match case-insensitively (a hyphen may be a space or absent); 1-2 characters match
        case-sensitively."""
        cands = {gene} if gene else set()
        if rec:
            cands |= {n for _, n in rec["names"]} | set(rec["protein_names"])
        for c in sorted(cands, key=len, reverse=True):
            if not c:
                continue
            if len(c) < 3:
                pat, fl = r"(?<![A-Za-z0-9])" + re.escape(c) + r"(?![A-Za-z0-9])", 0
            else:
                pat = r"(?<![A-Za-z0-9])" + r"[-\s]?".join(re.escape(x) for x in c.split("-")) + r"(?![A-Za-z0-9])"
                fl = re.I
            if re.search(pat, quote, flags=fl):
                return c
        return None

    def check_link(self, r):
        subj = f"{r['step_id']}|{r['reaction_id']}|{r['protein_accession'] or '-'}|{r['organism'].split(' (')[0]}"
        acc, tier = r["protein_accession"], r["tier"]
        rec = None
        if acc:
            try:
                rec = self.record(acc)
            except KbFetchError as e:
                self.add("link_record", subj, f"record {acc}", str(e)[:200], "fail" if not_found(e) else "unchecked")
            if rec:
                got, want = species_key(rec["organism"]), species_key(r["organism"])
                self.add("link_organism", subj, r["organism"], rec["organism"], "pass" if got == want else "fail")
                if r["gene_name"]:
                    hit = self.match_name(rec, r["gene_name"])
                    self.add("link_gene", subj, r["gene_name"],
                             f"{hit[0]} '{hit[1]}'" if hit else "no match; names in record: " + ", ".join(n for _, n in rec["names"]),
                             "pass" if hit else "fail")
        elif tier != "T4":
            self.add("link_record", subj, "protein accession", f"no protein accession in step_reaction.tsv (tier {tier})", "unchecked")
        pmid, quote = r["pmid"], r["quote"]
        if tier in ("T1a", "T1b", "T1c", "T2") and not (pmid and quote):
            miss = " and ".join(n for n, x in (("PMID", pmid), ("quote", quote)) if not x)
            self.add("link_evidence_present", subj, f"PMID and quote for tier {tier}", f"missing {miss}", "fail")
        if not pmid:
            return
        a = self.check_pmid_exists("link_pmid_exists", subj, pmid)
        if not quote:
            return
        # quote in abstract
        if a is None:
            self.add("link_quote_in_abstract", subj, "quote in abstract", f"abstract of PMID {pmid} not available (PubMed record not read)", "unchecked")
        elif not a["abstract"].strip():
            self.add("link_quote_in_abstract", subj, "quote in abstract", f"PubMed record {pmid} has no abstract", "unchecked")
        elif norm_text(quote) in norm_text(a["abstract"]):
            self.add("link_quote_in_abstract", subj, "quote in abstract", "found (whitespace-normalised substring)", "pass")
        else:
            self.add("link_quote_in_abstract", subj, quote[:160], "not a substring of the abstract of PMID " + pmid, "fail")
        one = is_one_sentence(quote)
        self.add("link_quote_one_sentence", subj, "one sentence", "one sentence" if one else "more than one sentence: " + quote[:120],
                 "pass" if one else "fail")
        if tier in ("T1a", "T1b", "T1c"):
            hit = self.names_quote(rec, r["gene_name"], quote)
            if hit:
                self.add("link_quote_names_protein", subj, "quote names the gene or protein", f"names '{hit}'", "pass")
            else:
                self.add("link_quote_names_protein", subj, f"gene '{r['gene_name']}' or a protein name of {acc or 'the record'}",
                         "quote does not name protein" + (" (flag for review, tier %s)" % tier if tier != "T1a" else ""),
                         "fail" if tier == "T1a" else "unchecked")

    def check_compound(self, r):
        cid, ik = r["chebi_id"], r["inchikey"]
        try:
            c = ks.chebi_compound(cid, fetch=self.fetch, cache_dir=self.cache_dir)
        except KbFetchError as e:
            for ch in ("chebi_formula", "chebi_inchikey"):
                self.add(ch, cid, "ChEBI record", str(e)[:200], "fail" if not_found(e) else "unchecked")
            self.add("pubchem_block1", cid, "PubChem record", "ChEBI record not read", "unchecked")
            return
        self.add("chebi_formula", cid, r["formula"], c["formula"], "pass" if c["formula"] == r["formula"] else "fail")
        if not ik:
            self.add("chebi_inchikey", cid, "InChIKey", f"compound.tsv has no InChIKey; ChEBI has {c['inchikey']}", "unchecked")
        else:
            self.add("chebi_inchikey", cid, ik, c["inchikey"], "pass" if c["inchikey"] == ik else "fail")
        qk = ik or c["inchikey"]
        try:
            pc = ks.pubchem_by_inchikey(qk, fetch=self.fetch, cache_dir=self.cache_dir)
        except KbFetchError as e:
            why = ("PubChem has no compound with this InChIKey" if not_found(e) else "PubChem not read")
            self.add("pubchem_block1", cid, qk[:14], f"{why}: {str(e)[:120]}", "unchecked")
            self.add("pubchem_formula", cid, r["formula"], f"{why}: {str(e)[:120]}", "unchecked")
            return
        ok = pc["inchikey"][:14] == qk[:14]
        self.add("pubchem_block1", cid, qk[:14], f"PubChem CID {pc['cid']}, {pc['inchikey'][:14]}", "pass" if ok else "fail")
        self.add("pubchem_formula", cid, r["formula"], f"PubChem CID {pc['cid']}, {pc['formula']}",
                 "pass" if pc["formula"] == r["formula"] else "fail")


def verify_pathway(yaml_path: Path, kb_dir: Path, *, fetch=None, refs_dir: Path | None = None,
                   cache_dir: Path | None = None) -> pd.DataFrame:
    """Run all checks. ``fetch`` (url -> bytes/str) replaces the network; None uses the real
    cached rate-limited fetchers of ``kbsources``."""
    p = load_pathway(Path(yaml_path), refs_dir)
    refs = Path(p.refs_dir)
    kb = load_kb(Path(kb_dir))
    v = Verifier(fetch, cache_dir)
    prov: dict[tuple[str, str], tuple[int, int]] = {}
    pf = refs / "provenance.tsv"
    if pf.is_file():
        for _, pr in pd.read_csv(pf, sep="\t", dtype=str, keep_default_na=False).iterrows():
            m = re.fullmatch(r"(\d+)-(\d+)", pr.get("region", "").strip())
            if m:
                prov[(pr["step_id"], pr["accession"])] = (int(m.group(1)), int(m.group(2)))
    for s in p.steps:
        fa = read_fasta(refs / f"{s.id}.faa")
        for r in s.refs:
            subj = f"{s.id}|{r.accession or '-'}|{r.species}"
            if r.accession:
                rec = v.check_record("ref", subj, r.accession, r.species, r.gene_name, fa, r.region, f"{s.id}.faa",
                                 prov.get((s.id, r.accession)))
                if rec and rec["kind"] == "GenBank":
                    v.check_genbank_cites(subj, rec, r.pmid)
            else:
                v.add("ref_exists", subj, "accession", "no accession in the YAML (see UNVERIFIED.md)", "unchecked")
            if r.pmid:
                v.check_pmid_exists("ref_pmid_exists", subj, r.pmid)
        dfa = read_fasta(refs / "decoys" / f"{s.id}.faa")
        for d in s.decoys:
            subj = f"{s.id}|{d.accession or '-'}|{d.species}"
            if not d.accession:
                v.add("decoy_exists", subj, "accession", "no accession in the YAML", "unchecked")
                continue
            try:
                rec = v.record(d.accession)
            except KbFetchError as e:
                v.add("decoy_exists", subj, f"record {d.accession}", str(e)[:200], "fail" if not_found(e) else "unchecked")
                continue
            v.add("decoy_exists", subj, f"record {d.accession}", f"{rec['kind']} record, organism {rec['organism']}", "pass")
            v.add("decoy_organism", subj, d.species, rec["organism"],
                  "pass" if species_key(rec["organism"]) == species_key(d.species) else "fail")
            if d.gene_name:
                hit = v.match_name(rec, d.gene_name)
                v.add("decoy_gene", subj, d.gene_name, f"{hit[0]} '{hit[1]}'" if hit else "no match; names in record: " + ", ".join(n for _, n in rec["names"]),
                      "pass" if hit else "fail")
            if dfa is None:
                v.add("decoy_sequence", subj, "FASTA record", f"no FASTA file decoys/{s.id}.faa", "unchecked")
            elif d.accession not in dfa:
                v.add("decoy_sequence", subj, "FASTA record", f"accession not in decoys/{s.id}.faa", "fail")
            else:
                same = norm_seq(dfa[d.accession]) == norm_seq(rec["sequence"])
                v.add("decoy_sequence", subj, f"sha256 {sha(rec['sequence'])[:12]}", f"FASTA sha256 {sha(dfa[d.accession])[:12]}", "pass" if same else "fail")
            if d.pmid:
                v.check_pmid_exists("decoy_pmid_exists", subj, d.pmid)
    for _, r in kb["step_reaction"].iterrows():
        v.check_link(r)
    for _, r in kb["compound"].iterrows():
        v.check_compound(r)
    return pd.DataFrame(v.rows, columns=COLS)


def render_report(results: dict[str, pd.DataFrame]) -> str:
    out = ["# Knowledge-base verification report", "",
           "Generated by `bin/verify_kb.py` from the written YAML, FASTA and TSV files and the public "
           "records (UniProt, NCBI protein, PubMed, ChEBI, PubChem), read through the cached "
           "fetchers of `pathway_census.kbsources`. Re-run: `/usr/bin/python3.12 bin/verify_kb.py`.", "",
           "::: {.callout-important}", "## Limits of this check", "",
           "- A pass means the public record agrees with the file. It does not mean the biology is right.",
           "- A quote check proves the sentence is in the PubMed abstract. It does not prove the sentence supports the tier.",
           "- `unchecked` means no source was read. Each reason is listed below. None of them is a pass.",
           "- `quote does not name protein` is a fail for tier T1a and an unchecked review flag for T1b and T1c.",
           "- Gene names also match ORF names and ordered locus names. The kind that matched is in the `found` column.",
           "- The protein-name match uses UniProt protein names. A generic name (for example `Polyketide synthase`) matches any quote that contains it, so a T1 quote about a different enzyme of the same kind can pass. Known limit.",
           "- `pubchem_block1` and `pubchem_formula` compare compound.tsv with the PubChem compound found by the same InChIKey. They show that PubChem has that compound and agrees on the first block and the formula.",
           "- Tiers, quotes and step sets were not changed by this check.", ":::", "",
           "## Counts per pathway", "", "| pathway | checks | pass | fail | unchecked |", "|---|---|---|---|---|"]
    tot = [0, 0, 0, 0]
    for name, df in results.items():
        c = df.status.value_counts()
        row = [len(df), int(c.get("pass", 0)), int(c.get("fail", 0)), int(c.get("unchecked", 0))]
        tot = [a + b for a, b in zip(tot, row)]
        out.append(f"| {name} | " + " | ".join(map(str, row)) + " |")
    out.append("| all | " + " | ".join(map(str, tot)) + " |")
    out += ["", "## Counts per check", "", "| check | pass | fail | unchecked |", "|---|---|---|---|"]
    allr = pd.concat([d.assign(pathway=n) for n, d in results.items()], ignore_index=True) if results else pd.DataFrame(columns=COLS + ["pathway"])
    if len(allr):
        piv = allr.pivot_table(index="check", columns="status", values="subject", aggfunc="count", fill_value=0)
        for chk, r in piv.iterrows():
            out.append(f"| {chk} | {int(r.get('pass', 0))} | {int(r.get('fail', 0))} | {int(r.get('unchecked', 0))} |")

    def md(s):
        return str(s).replace("|", "\\|").replace("\n", " ")

    out += ["", "## Gene-name passes by kind", "",
            "| pathway | gene name | gene name (synonym) | ORF name | ordered locus name |", "|---|---|---|---|---|"]
    kinds = ["gene name", "gene name (synonym)", "ORF name", "ordered locus name"]
    for name, df in results.items():
        g = df[df.check.isin(["ref_gene", "link_gene", "decoy_gene"]) & (df.status == "pass")]
        k = g.found.str.split(" '").str[0].value_counts()
        out.append(f"| {name} | " + " | ".join(str(int(k.get(x, 0))) for x in kinds) + " |")
    out.append("")
    out.append("Rows are ref_gene, link_gene and decoy_gene passes together.")

    fails = allr[allr.status == "fail"] if len(allr) else allr
    out += ["", f"## Failures ({len(fails)})", ""]
    if len(fails):
        out += ["| pathway | check | subject | expected | found |", "|---|---|---|---|---|"]
        out += [f"| {r.pathway} | {r.check} | {md(r.subject)} | {md(r.expected)} | {md(r.found)} |" for r in fails.itertuples()]
    else:
        out.append("None.")
    unc = allr[allr.status == "unchecked"] if len(allr) else allr
    out += ["", f"## Unchecked ({len(unc)})", ""]
    if len(unc):
        out += ["| pathway | check | subject | reason |", "|---|---|---|---|"]
        out += [f"| {r.pathway} | {r.check} | {md(r.subject)} | {md(r.found)} |" for r in unc.itertuples()]
    else:
        out.append("None.")
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pathways", nargs="*", default=list(PATHWAYS))
    ap.add_argument("--report", type=Path, default=PKG / "kb" / "VERIFICATION_REPORT.md")
    ap.add_argument("--strict", action="store_true", help="exit 1 when any check fails")
    a = ap.parse_args(argv)
    results = {}
    for name in a.pathways:
        results[name] = verify_pathway(PKG / "pathways" / f"{name}.yaml", PKG / "kb" / name)
        c = results[name].status.value_counts()
        print(f"{name}: {len(results[name])} checks, pass {c.get('pass', 0)}, fail {c.get('fail', 0)}, unchecked {c.get('unchecked', 0)}")
    a.report.write_text(render_report(results), encoding="utf-8")
    print("report:", a.report)
    nfail = sum(int((d.status == "fail").sum()) for d in results.values())
    return 1 if (a.strict and nfail) else 0


if __name__ == "__main__":
    sys.exit(main())
