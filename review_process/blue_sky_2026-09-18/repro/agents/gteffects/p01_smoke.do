clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"

which wooldid

* ---------------- A: balanced, clean -------------------------------------
mkpanel , n(60) t(10) onsets("2004 2006 2008") ysd(0.3)
dis "N=" _N
cap noi wooldid y unit year ttre, att att_it att_i att_itime gteffects esprelength(3) espostlength(3)
dis "A_RC=" _rc
dis "--- e(b) colnames ---"
dis "`: colnames e(b)'"
dis "--- gteffects colnames ---"
dis "`: colnames e(gteffects)'"
gtcheck , tag(A_balanced_unit) pre es
gtdefcheck , tag(A_balanced_unit)

* ---------------- B: i = cohort ------------------------------------------
cap noi wooldid y ttre year ttre, att att_it att_i att_itime gteffects esprelength(3) espostlength(3)
dis "B_RC=" _rc
gtcheck , tag(B_balanced_cohort) pre es
gtdefcheck , tag(B_balanced_cohort)

dis "DONE_P01"
