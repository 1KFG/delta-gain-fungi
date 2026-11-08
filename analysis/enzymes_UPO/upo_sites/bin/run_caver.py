#!/usr/bin/env python3.12
"""
Run CAVER 3.0.2 on a peroxygenase structure, starting from the heme iron, and report the
bottleneck radius of the heme access channel.

This exists to reproduce the measurement behind Karich et al. 2017's ring-size rule: AaeUPO has
a ~7 A heme access channel and converts PAHs up to about 6 A minimum width, while MroUPO has a
~11 A channel and handles bulkier substrates. Our own hard-sphere bottleneck measure
(heme_channel.py) reproduced AaeUPO within 1 A but gave MroUPO 6.4 A instead of 11 A, so it could
not reproduce the contrast the rule depends on. CAVER is the tool the published numbers came from.

Input PDBs are cleaned first: one chain, protein atoms plus the heme, waters and buffer
heteroatoms removed. The heme is kept because it forms one wall of the channel.
"""
from __future__ import annotations
import argparse, csv, os, shutil, subprocess, sys, tempfile

SOLVENT_HET = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD",
               "CL", "NA", "K", "ACT", "DMS", "TRS", "FMT", "NO3", "CA", "MN"}
HEME = {"HEM", "HEC", "HEB"}


def clean_pdb(src, dst, chain="A"):
    """Write chain + heme only. Returns the heme iron coordinates."""
    fe = None
    with open(src) as fh, open(dst, "w") as out:
        for line in fh:
            rec = line[:6]
            if rec not in ("ATOM  ", "HETATM"):
                continue
            if line[16] not in (" ", "A"):
                continue
            resn = line[17:20].strip()
            if line[21] != chain:
                continue
            if rec == "HETATM" and (resn in SOLVENT_HET or resn not in HEME):
                continue
            el = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
            if el in ("H", "D"):
                continue
            if el == "FE" and resn in HEME:
                fe = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            out.write(line)
        out.write("END\n")
    return fe


def parse_tunnels(outdir):
    """Return (bottleneck_radius, throughput) of the best tunnel, or (None, None)."""
    cands = [
        os.path.join(outdir, "analysis", "tunnel_characteristics.csv"),
        os.path.join(outdir, "analysis", "tunnel_profiles.csv"),
    ]
    for p in cands:
        if not os.path.exists(p):
            continue
        with open(p, newline="") as fh:
            rows = list(csv.reader(fh))
        if not rows:
            continue
        hdr = [h.strip().lower() for h in rows[0]]
        bi = next((i for i, h in enumerate(hdr) if "bottleneck" in h and "radius" in h), None)
        ti = next((i for i, h in enumerate(hdr) if "throughput" in h), None)
        if bi is None:
            continue
        best = None
        for r in rows[1:]:
            if len(r) <= bi:
                continue
            try:
                b = float(r[bi])
                t = float(r[ti]) if ti is not None and len(r) > ti else 0.0
            except ValueError:
                continue
            if best is None or t > best[1]:
                best = (b, t)
        if best:
            return best
    return (None, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", nargs="+", required=True)
    ap.add_argument("--label", nargs="*", default=None)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--caver-home", required=True, help=".../caver_3.0/caver")
    ap.add_argument("--probe", type=float, default=0.9)
    ap.add_argument("--shell-radius", type=float, default=3.0)
    ap.add_argument("--shell-depth", type=float, default=4.0)
    ap.add_argument("--xmx", default="4000m")
    ap.add_argument("--out", default=None)
    ap.add_argument("--keep", default=None, help="keep CAVER output under this directory")
    a = ap.parse_args()

    jar = os.path.join(a.caver_home, "caver.jar")
    lib = os.path.join(a.caver_home, "lib")
    for p in (jar, lib):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    labels = a.label if a.label else [os.path.basename(p) for p in a.pdb]
    results = []
    print("structure\tbottleneck_radius_A\tchannel_diameter_A\tthroughput", flush=True)
    for pdb, lab in zip(a.pdb, labels):
        work = tempfile.mkdtemp(prefix="caver.")
        pdbdir = os.path.join(work, "in"); os.makedirs(pdbdir)
        fe = clean_pdb(pdb, os.path.join(pdbdir, "model.pdb"), a.chain)
        if fe is None:
            print(f"{lab}\tNO_HEME_FE\t\t", flush=True)
            shutil.rmtree(work, ignore_errors=True); continue
        conf = os.path.join(work, "config.txt")
        with open(conf, "w") as fh:
            fh.write(f"starting_point_coordinates {fe[0]:.3f} {fe[1]:.3f} {fe[2]:.3f}\n")
            fh.write(f"probe_radius {a.probe}\n")
            fh.write(f"shell_radius {a.shell_radius}\n")
            fh.write(f"shell_depth {a.shell_depth}\n")
            fh.write("seed 1\n")
            fh.write("generate_summary yes\n")
            fh.write("save_dynamics_visualization no\n")
            fh.write("compute_tunnel_residues no\n")
        outd = os.path.join(work, "out")
        r = subprocess.run(["java", f"-Xmx{a.xmx}", "-cp", lib, "-jar", jar,
                            "-home", a.caver_home, "-pdb", pdbdir,
                            "-conf", conf, "-out", outd],
                           capture_output=True, text=True, timeout=3600)
        b, t = parse_tunnels(outd)
        if b is None:
            tail = (r.stderr or r.stdout).strip().splitlines()[-3:]
            print(f"{lab}\tFAILED\t\t  # {' | '.join(tail)}", flush=True)
        else:
            print(f"{lab}\t{b:.2f}\t{2*b:.2f}\t{t:.3g}", flush=True)
            results.append((lab, b, 2 * b, t))
        if a.keep:
            os.makedirs(a.keep, exist_ok=True)
            dst = os.path.join(a.keep, lab.replace("/", "_"))
            shutil.rmtree(dst, ignore_errors=True)
            if os.path.exists(outd):
                shutil.copytree(outd, dst)
        shutil.rmtree(work, ignore_errors=True)

    if a.out:
        with open(a.out, "w") as fh:
            fh.write("structure\tbottleneck_radius_A\tchannel_diameter_A\tthroughput\n")
            for lab, b, d, t in results:
                fh.write(f"{lab}\t{b:.3f}\t{d:.3f}\t{t:.4g}\n")


if __name__ == "__main__":
    main()
