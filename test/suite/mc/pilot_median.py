#!/usr/bin/env python3
"""Print the pilot median of mdt_wald or mdt_nis from ptp_pilot_s*.csv (rc==0 rows)."""
import csv
import glob
import statistics
import sys

outdir, col = sys.argv[1], sys.argv[2]
key = {"mdt_wald": "mdt_wald", "mdt_nis": "mdt_nis"}[col]
vals = []
for f in sorted(glob.glob(f"{outdir}/ptp_pilot_s*.csv")):
    with open(f) as fh:
        for row in csv.DictReader(fh):
            if row["rc"] == "0" and row[key] not in (".", ""):
                vals.append(float(row[key]))
if not vals:
    sys.exit(1)
print(f"{statistics.median(vals):.10g}")
