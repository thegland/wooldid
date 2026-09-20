clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260922
set sortseed 20260922
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-estimands"

use "dataE.dta", clear
dis "=== third engine arm: unconditionalse (also explicit dummies) ==="
cap noi wooldid y i t ttre, att semielasticity
dis "J1_RC=" _rc
dis "J1_default_reghdfe=" %21x e(b)[1,1]
use "dataE.dta", clear
cap noi wooldid y i t ttre, att semielasticity customvce(robust)
dis "J2_RC=" _rc
dis "J2_customvce=" %21x e(b)[1,1]
use "dataE.dta", clear
cap noi wooldid y i t ttre, att semielasticity unconditionalse
dis "J3_RC=" _rc
dis "J3_unconditionalse=" %21x e(b)[1,1]

dis "=== is there ANY correct path for a cell-weighted estimand + semielasticity? ==="
use "dataE.dta", clear
cap noi wooldid y i t ttre, att_it semielasticity
dis "K1_RC=" _rc
dis "K1_default_att_it=" %21x e(b)[1,1]
use "dataE.dta", clear
cap noi wooldid y i t ttre, att_it semielasticity unconditionalse
dis "K2_RC=" _rc
use "dataE.dta", clear
cap noi wooldid y i t ttre, att_it semielasticity customvce(robust)
dis "K3_RC=" _rc
dis "K3_customvce_att_it=" %21x e(b)[1,1]

dis "=== does anything caution on the default semielasticity run? ==="
use "dataE.dta", clear
cap noi wooldid y i t ttre, att semielasticity verbose
dis "L1_RC=" _rc
dis "ALLDONE_RUN7"
