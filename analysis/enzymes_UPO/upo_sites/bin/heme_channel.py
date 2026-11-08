#!/usr/bin/env python3.12
"""
Measure the heme access channel bottleneck of a peroxygenase structure.

WHY THIS EXISTS, AND WHY IT IS NOT CAVER
----------------------------------------
Karich et al. 2017 explain PAH substrate range mechanically: AaeUPO has a ~7 A heme access
channel and converts PAHs up to about 6 A minimum width, while MroUPO has a ~11 A channel and
handles bulkier ones. That is a ring-size rule, and applying it to BFD enzymes needs a channel
width per enzyme.

CAVER 3.0 now requires registration (its direct download returns 404), so rather than block on
an install this computes the same quantity directly. The measurement is validated against the two
published structures whose widths Karich reports -- see --validate.

METHOD
------
1. Protein heavy atoms and the heme are treated as hard spheres with van der Waals radii.
   Waters and other solvent heteroatoms are discarded; the heme is KEPT, since it forms one wall
   of the channel.
2. A regular grid is built over the structure. A grid point is "solid" if it lies inside any
   atom's van der Waals sphere.
3. An exact Euclidean distance transform of the free space gives, at every point, the radius of
   the largest probe sphere that fits there (the clearance).
4. The bottleneck radius is the largest r for which a connected path of points with
   clearance >= r joins the heme iron to bulk solvent (the grid boundary). This is found by
   binary search on r using connected-component labelling, which is exact to grid resolution and
   much faster than a shortest-path search.

Reported channel width is the DIAMETER, 2 x bottleneck radius, to match how Karich reports it.

LIMITATIONS, read before using a number
---------------------------------------
- A single static structure. Real channels breathe; side chains move. This is the closed-state
  bottleneck, which will read narrower than a substrate actually experiences.
- Grid resolution sets the precision; the default 0.4 A grid gives about +/- 0.4 A on the
  diameter.
- The published 7 A / 11 A values are themselves from CAVER with its own probe and parameters.
  Agreement with them is the calibration (--validate), not a guarantee of transferability.
- Only two enzymes anchor the width-to-activity relation. Any ring-size prediction from this
  carries that n=2 calibration and must say so.
"""
from __future__ import annotations
import argparse, sys
import numpy as np
from scipy import ndimage

# Bondi van der Waals radii (A)
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "FE": 2.00,
       "MG": 1.73, "ZN": 1.39, "SE": 1.90, "H": 1.20}
DEFAULT_R = 1.70
SOLVENT_HET = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD",
               "CL", "NA", "K", "ACT", "DMS", "TRS", "FMT", "NO3", "CA", "MN"}


def parse_pdb(path, chain=None, keep_het=("HEM", "HEC", "HEB")):
    """Return (coords, radii, fe_coord). Waters and buffer heteroatoms are dropped."""
    xyz, rad, fe = [], [], None
    with open(path) as fh:
        for line in fh:
            rec = line[:6]
            if rec not in ("ATOM  ", "HETATM"):
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            resn = line[17:20].strip()
            ch = line[21]
            if chain and ch != chain:
                continue
            if rec == "HETATM":
                if resn in SOLVENT_HET or resn not in keep_het:
                    continue
            el = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
            if el == "H" or el == "D":
                continue
            c = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            xyz.append(c)
            rad.append(VDW.get(el, DEFAULT_R))
            if el == "FE" and resn in keep_het:
                fe = c
    if not xyz:
        raise SystemExit(f"no atoms parsed from {path} (chain={chain})")
    return np.array(xyz), np.array(rad), fe


def bottleneck(xyz, rad, fe, grid=0.4, pad=8.0, rmax=12.0, tol=0.02,
               seed_radius=9.0, verbose=False):
    """Largest probe radius connecting the heme iron to bulk solvent."""
    lo = xyz.min(0) - (rad.max() + pad)
    hi = xyz.max(0) + (rad.max() + pad)
    shape = np.ceil((hi - lo) / grid).astype(int) + 1
    if shape.prod() > 60_000_000:
        raise SystemExit(f"grid too large ({shape}); increase --grid")

    solid = np.zeros(shape, dtype=bool)
    # rasterise each atom's vdW sphere
    for c, r in zip(xyz, rad):
        i0 = np.floor((c - r - lo) / grid).astype(int)
        i1 = np.ceil((c + r - lo) / grid).astype(int) + 1
        i0 = np.maximum(i0, 0); i1 = np.minimum(i1, shape)
        if np.any(i1 <= i0):
            continue
        sl = tuple(slice(a, b) for a, b in zip(i0, i1))
        g = np.ogrid[tuple(slice(a, b) for a, b in zip(i0, i1))]
        d2 = sum((lo[k] + g[k] * grid - c[k]) ** 2 for k in range(3))
        solid[sl] |= d2 <= r * r

    # clearance: distance from each free point to the nearest solid point
    clearance = ndimage.distance_transform_edt(~solid, sampling=grid)

    # Start point: the DISTAL POCKET above the heme, not merely the closest free voxel.
    # The iron is buried under the porphyrin, so the nearest free voxel is usually a sub-angstrom
    # crevice between heme atoms and makes the whole search report one grid unit. Instead take the
    # free voxel within `seed_radius` of the iron that has the LARGEST clearance: that is the
    # substrate-binding cavity the channel actually leads to.
    fi0 = np.round((fe - lo) / grid).astype(int)
    fi0 = np.clip(fi0, 0, shape - 1)
    rr = int(np.ceil(seed_radius / grid))
    sl = tuple(slice(max(0, fi0[k] - rr), min(shape[k], fi0[k] + rr + 1)) for k in range(3))
    sub = clearance[sl]
    if sub.size == 0 or sub.max() <= 0:
        return 0.0
    loc = np.unravel_index(np.argmax(sub), sub.shape)
    fi = np.array([sl[k].start + loc[k] for k in range(3)])
    seed_clear = clearance[tuple(fi)]

    def connected(r):
        mask = clearance >= r
        if not mask[tuple(fi)]:
            return False
        lab, _ = ndimage.label(mask)
        start = lab[tuple(fi)]
        if start == 0:
            return False
        b = np.concatenate([lab[0].ravel(), lab[-1].ravel(),
                            lab[:, 0].ravel(), lab[:, -1].ravel(),
                            lab[:, :, 0].ravel(), lab[:, :, -1].ravel()])
        return bool(np.any(b == start))

    if verbose:
        print(f"    seed clearance {seed_clear:.2f} A at {fi}, "
              f"max clearance in box {clearance.max():.2f} A", file=sys.stderr)
    if not connected(0.0):
        return 0.0
    lo_r, hi_r = 0.0, rmax
    while hi_r - lo_r > tol:
        mid = 0.5 * (lo_r + hi_r)
        if connected(mid):
            lo_r = mid
        else:
            hi_r = mid
    return lo_r


def run(path, chain, grid, label=None):
    xyz, rad, fe = parse_pdb(path, chain)
    if fe is None:
        print(f"{label or path}\tNO_HEME_FE", flush=True)
        return None
    r = bottleneck(xyz, rad, fe, grid=grid, verbose=True)
    print(f"{label or path}\t{len(xyz)}\t{2*r:.2f}", flush=True)
    return 2 * r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", nargs="*", default=[])
    ap.add_argument("--chain", default="A")
    ap.add_argument("--grid", type=float, default=0.4)
    ap.add_argument("--label", nargs="*", default=None)
    ap.add_argument("--validate", action="store_true",
                    help="measure the two structures whose channel widths Karich 2017 reports "
                         "(AaeUPO 2YOR ~7 A, MroUPO 5FUJ ~11 A) and show the comparison")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    rows = []
    if a.validate:
        import os
        here = os.path.dirname(os.path.abspath(__file__))
        ref = os.path.join(here, "..", "..", "data", "structures", "ref")
        print("structure\tatoms\tchannel_diameter_A", flush=True)
        exp = {"2yor": ("AaeUPO", 7.0), "5fuj": ("MroUPO", 11.0)}
        got = {}
        for p, (nm, e) in exp.items():
            f = os.path.join(ref, f"{p}.pdb")
            if not os.path.exists(f):
                print(f"{nm}\tMISSING {f}"); continue
            got[nm] = run(f, a.chain, a.grid, label=f"{nm} ({p.upper()})")
        print("\ncalibration against Karich et al. 2017:")
        for nm, (_, e) in [(v[0], (k, v[1])) for k, v in exp.items()]:
            if nm in got and got[nm] is not None:
                print(f"  {nm:8} measured {got[nm]:5.2f} A   published ~{e:4.1f} A   "
                      f"difference {got[nm]-e:+.2f} A")
        return

    labels = a.label if a.label else [None] * len(a.pdb)
    print("structure\tatoms\tchannel_diameter_A", flush=True)
    for p, l in zip(a.pdb, labels):
        rows.append((l or p, run(p, a.chain, a.grid, label=l)))
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("structure\tchannel_diameter_A\n")
            for n, v in rows:
                fh.write(f"{n}\t{'' if v is None else f'{v:.2f}'}\n")


if __name__ == "__main__":
    main()
