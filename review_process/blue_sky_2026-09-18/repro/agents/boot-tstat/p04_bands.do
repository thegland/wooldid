clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid

local gG 20
local gUnits 3
local gT 8
local gonsets "2003 2004 2005"
local gnever 5
local gskew 1
include "`WD'/gen_small.do"
dis "NOBS=" _N
save skew20, replace

*==================================================================
* S1: bootstrap, NO boottstat, default display(percentile), skewed y.
*     Is the simultaneous band ever NARROWER than the pointwise interval?
*==================================================================
dis "=== S1 percentile display + band ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(399) rseed(555) att esprelength(3) espostlength(3) simultaneousbands dots(0)
dis "S1_RC=" _rc
dis "S1_ENGINE=" "`e(esband_engine)'"
dis "S1_FAMS=" "`e(esband_families)'"
dis "S1_CRIT=" e(esband_crit_att)
local BFAM att
local BTAG S1
include "`WD'/bandcmp.do"

*==================================================================
* S2: same run but bootdisplay(wald) -- b +/- zc*sd, the family the band
*     is built to match.
*==================================================================
dis "=== S2 wald display + band ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(399) rseed(555) att esprelength(3) espostlength(3) simultaneousbands bootdisplay(wald) dots(0)
dis "S2_RC=" _rc
dis "S2_ENGINE=" "`e(esband_engine)'"
dis "S2_CRIT=" e(esband_crit_att)
local BFAM att
local BTAG S2
include "`WD'/bandcmp.do"

*==================================================================
* S3: boottstat -> studentized band
*==================================================================
dis "=== S3 boottstat + band ==="
cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(399) rseed(555) boottstat att esprelength(3) espostlength(3) simultaneousbands dots(0)
dis "S3_RC=" _rc
dis "S3_ENGINE=" "`e(esband_engine)'"
dis "S3_CRIT=" e(esband_crit_att)
dis "S3_N=" e(esband_n_att)
local BFAM att
local BTAG S3
include "`WD'/bandcmp.do"

*==================================================================
* S4: ANALYTIC band -- family independence (D-9).  att alone vs att+att_it.
*==================================================================
dis "=== S4a analytic band, att only ==="
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(3) simultaneousbands
dis "S4a_RC=" _rc
dis "S4a_ENGINE=" "`e(esband_engine)'"
dis "S4a_FAMS=" "`e(esband_families)'"
mat S4A = e(esband_att)
scalar c4a = e(esband_crit_att)
dis "S4a_CRIT=" %21x c4a

dis "=== S4b analytic band, att + att_it ==="
cap noi wooldid y state t ttre , cluster(state) att att_it esprelength(3) espostlength(3) simultaneousbands
dis "S4b_RC=" _rc
dis "S4b_FAMS=" "`e(esband_families)'"
mat S4B = e(esband_att)
scalar c4b = e(esband_crit_att)
dis "S4b_CRIT=" %21x c4b
dis "S4_CRIT_EQUAL=" (c4a == c4b)
mata:
  A = st_matrix("S4A"); B = st_matrix("S4B")
  if (rows(A)==rows(B) & cols(A)==cols(B)) printf("S4_BANDMAXDIFF=%g\n", max(abs(A-B)))
  else printf("S4_SHAPE_MISMATCH\n")
end

dis "=== S4c analytic band, att + subgroups ==="
gen byte reg2 = mod(state,2)
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(3) simultaneousbands subgroups(reg2)
dis "S4c_RC=" _rc
dis "S4c_FAMS=" "`e(esband_families)'"
scalar c4c = e(esband_crit_att)
dis "S4c_CRIT=" %21x c4c
dis "S4_CRIT_EQUAL_AC=" (c4a == c4c)
dis "P04_DONE"
