clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"
local OO "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/out"
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/lib_dump.do"

*==================== NOISE FLOOR ====================
* same arm twice: must be bit-identical
qui use "`DD'/dC.dta", clear
cap noisily wooldid y cohort t ttre, att att_i att_it att_itime jointtests esprelength(3) espostlength(3) nofastcollapse
dis "NF1_RC=" _rc
wddump , file("`OO'/n0.A.txt")
qui use "`DD'/dC.dta", clear
cap noisily wooldid y cohort t ttre, att att_i att_it att_itime jointtests esprelength(3) espostlength(3) nofastcollapse
dis "NF2_RC=" _rc
wddump , file("`OO'/n0.B.txt")

* and the collapse arm twice
qui use "`DD'/dC.dta", clear
cap noisily wooldid y cohort t ttre, att att_i att_it att_itime jointtests esprelength(3) espostlength(3)
wddump , file("`OO'/n1.A.txt")
qui use "`DD'/dC.dta", clear
cap noisily wooldid y cohort t ttre, att att_i att_it att_itime jointtests esprelength(3) espostlength(3)
wddump , file("`OO'/n1.B.txt")

*==================== BASE SWEEP (i = cohort) ====================
local BASE wooldid y cohort t ttre
local EST att att_i att_it att_itime jointtests esprelength(3) espostlength(3)

wdcell , dta("`DD'/dC.dta") tag(c01_plain)   base(`BASE') out("`OO'") opts(`EST')
wdcell , dta("`DD'/dC.dta") tag(c02_robust)  base(`BASE') out("`OO'") opts(`EST' robust)
wdcell , dta("`DD'/dC.dta") tag(c03_clreg)   base(`BASE') out("`OO'") opts(`EST' cluster(region))
wdcell , dta("`DD'/dC.dta") tag(c04_aw)      base(wooldid y cohort t ttre [aw=w])   out("`OO'") opts(`EST')
wdcell , dta("`DD'/dC.dta") tag(c05_fw)      base(wooldid y cohort t ttre [fw=fwt]) out("`OO'") opts(`EST')
wdcell , dta("`DD'/dC.dta") tag(c06_pw)      base(wooldid y cohort t ttre [pw=w])   out("`OO'") opts(`EST')
wdcell , dta("`DD'/dC.dta") tag(c07_impx)    base(`BASE') out("`OO'") opts(att impvalues cluster(region) esprelength(3) espostlength(3))
wdcell , dta("`DD'/dC.dta") tag(c08_impnest) base(`BASE') out("`OO'") opts(att impvalues cluster(clnest) esprelength(3) espostlength(3))
wdcell , dta("`DD'/dC.dta") tag(c09_ctrl)    base(`BASE') out("`OO'") opts(`EST' controls(x1 x2))
wdcell , dta("`DD'/dC.dta") tag(c10_sg)      base(`BASE') out("`OO'") opts(att subgroups(sgc) jointtests esprelength(3) espostlength(3))
dis "P01_DONE"
