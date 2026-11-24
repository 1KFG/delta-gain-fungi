# refs

`mibig_json_4.0.tar.gz`: the MIBiG 4.0 entry JSON files (3,013 files, 15 MB unpacked, 0.95 MB packed). `build_ebc_table.py` reads them for compound names, PubChem ids, quality and status.

Use:
    tar -C $SCRATCH -xzf refs/mibig_json_4.0.tar.gz
    python3 build_ebc_table.py results $SCRATCH/mibig_json_4.0 out.tsv

Source: copied on 2026-11-23 from a session scratchpad copy. Where that copy came from is not recorded here (UNVERIFIED); compare with the MIBiG download before citing.
