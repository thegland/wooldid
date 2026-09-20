*! n05 — SKEP-N3 repro (read-only; no estimation).  The 2026-09-16 prereg retracts an
*! explanation "so it does not propagate"; the file that records the same measurement
*! still carries it.  Run this and read the two greps side by side.
clear all
set more off
dis "===== the RETRACTION (PREREG_FEATURE_FINECLUS_2026-09-16.md)"
!grep -n "Do not cite that passage" /Users/hegland/Documents/Research/wooldid_dev/test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md
!grep -n "A wrong analogy, retracted here" /Users/hegland/Documents/Research/wooldid_dev/test/suite/PREREG_FEATURE_FINECLUS_2026-09-16.md
dis "===== the STALE COPY (test/suite/mc/mc_fineclus.do, header)"
!grep -n "same effect the help already describes" /Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/mc_fineclus.do
!grep -n "stops seeing the treated cells" /Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/mc_fineclus.do
dis "===== DONE"
