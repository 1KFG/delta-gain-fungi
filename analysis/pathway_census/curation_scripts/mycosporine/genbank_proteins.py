"""GenBank protein records of the Phaffia rhodozyma UCD 67-385 mycosporine cluster (nucleotide
record OQ547787; proteins WEV88679.1 DDGS, WEV88680.1 OMT, WEV88681.1 ATPG). UniProtKB has no
entry for these three sequences (UniParc lists them without a UniProtKB accession), so the
records are read from NCBI through the cached fetch.

Import: gp(acc) -> dict(accession, gene, product, organism, strain, taxid, length, sequence,
coded_by, title, journal). Run: prints the records and, with --fasta FILE, writes a FASTA file."""
from helpers import *
import re

ACCS = ("WEV88679.1", "WEV88680.1", "WEV88681.1")

def gp_text(acc):
    return raw("ncbi_protein", acc.split(".")[0] + "_gp",
               f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={acc.split('.')[0]}&rettype=gp&retmode=text", as_json=False)

def gp(acc):
    t = gp_text(acc)
    def q(key):
        m = re.search(r'/%s="([^"]*)"' % key, t)
        return " ".join(m.group(1).split()) if m else None
    ver = re.search(r"^VERSION\s+(\S+)", t, flags=re.M).group(1)
    assert ver == acc, (acc, ver)
    seq = re.sub(r"[^a-z]", "", t.split("ORIGIN")[1].split("//")[0]).upper()
    length = int(re.search(r"^LOCUS\s+\S+\s+(\d+) aa", t, flags=re.M).group(1))
    assert len(seq) == length, (acc, len(seq), length)
    refs = re.findall(r"TITLE\s+(.*?)\n\s+JOURNAL\s+(.*?)\n(?:\s+REMARK\s+(.*?)\n)?", t, flags=re.S)
    return {"accession": ver, "gene": q("gene"), "product": q("product"), "organism": q("organism"),
            "strain": q("strain"), "taxid": re.search(r'/db_xref="taxon:(\d+)"', t).group(1),
            "length": length, "sequence": seq, "coded_by": q("coded_by"),
            "refs": [tuple(" ".join((x or "").split()) for x in r) for r in refs]}

if __name__ == "__main__":
    out = []
    for a in ACCS:
        d = gp(a)
        print({k: v for k, v in d.items() if k != "sequence"})
        out.append(f">{a} {d['gene']} {d['organism']} {d['strain']}\n{d['sequence']}\n")
    if "--fasta" in sys.argv:
        open(sys.argv[sys.argv.index("--fasta") + 1], "w").write("".join(out))
