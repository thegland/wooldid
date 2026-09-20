* probe3: dump e(escomp) and e(gteffects) side by side, aw=k vs physical duplication,
* so the divergences probe2 flagged can be attributed column by column.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

use "p2_A.dta", clear
qui wooldid y i t ttre [aw=k], att att_it esprelength(2) espostlength(3)
matrix EA = e(escomp)
dis "A_refN=" e(escomp_refN) " A_refW=" e(escomp_refW) " A_trtN=" e(escomp_trtN) " A_trtW=" e(escomp_trtW)
dis "ESCOMP_A colnames: " "`: colnames EA'"
matlist EA, format(%14.6g) names(all)

use "p2_B.dta", clear
qui wooldid y i t ttre , att att_it esprelength(2) espostlength(3)
matrix EB = e(escomp)
dis "B_refN=" e(escomp_refN) " B_refW=" e(escomp_refW) " B_trtN=" e(escomp_trtN) " B_trtW=" e(escomp_trtW)
matlist EB, format(%14.6g) names(all)

* now gteffects
use "p2_A.dta", clear
qui wooldid y i t ttre [aw=k], att gteffects
matrix GA = e(gteffects)
dis "GTE_A colnames: " "`: colnames GA'"
matlist GA, format(%14.6g) names(all)

use "p2_B.dta", clear
qui wooldid y i t ttre , att gteffects
matrix GB = e(gteffects)
matlist GB, format(%14.6g) names(all)

dis "PROBE3_DONE"
