clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 120
which wooldid

*----- helper: half-width ratio check on a matrix with estimate/se/lb/ub -----
cap prog drop hw
prog def hw
  syntax , M1(name) M2(name) C1(string) C2(string) tag(string)
  local n1 : colfullnames `M1'
  local n2 : colfullnames `M2'
  dis "@@`tag'_names1 [`n1']"
  dis "@@`tag'_names2 [`n2']"
  local ie  = 0
  local ilb1 = 0
  local iub1 = 0
  local k = 0
  foreach nm of local n1 {
    local ++k
    if "`nm'" == "estimate" local ie = `k'
    if "`nm'" == "`c1'" local ilb1 = `k'
  }
  dis "@@`tag'_idx e=`ie' lb1=`ilb1'"
end

*=========== D: gteffects intervals under level(95) vs level(80) ===========
use "dat.dta", clear
qui wooldid y uid t ttre, att robust gteffects level(95)
mat G95 = e(gteffects)
scalar df = e(df_r)
qui wooldid y uid t ttre, att robust gteffects level(80)
mat G80 = e(gteffects)
local n95 : colfullnames G95
local n80 : colfullnames G80
dis "@@GN95 [`n95']"
dis "@@GN80 [`n80']"
scalar t95 = invttail(df, .025)
scalar t80 = invttail(df, .10)
dis "@@GDF df=" df " t95=" t95 " t80=" t80 " target_ratio=" t80/t95
* locate columns by name
local ke = 0
local k95lb = 0
local k95ub = 0
local kse = 0
local k = 0
foreach nm of local n95 {
  local ++k
  if "`nm'"=="estimate" local ke = `k'
  if "`nm'"=="se" local kse = `k'
  if "`nm'"=="lb_95ci" local k95lb = `k'
  if "`nm'"=="ub_95ci" local k95ub = `k'
}
local k = 0
local k80lb = 0
local k80ub = 0
foreach nm of local n80 {
  local ++k
  if "`nm'"=="lb_80ci" local k80lb = `k'
  if "`nm'"=="ub_80ci" local k80ub = `k'
}
dis "@@GIDX e=`ke' se=`kse' lb95=`k95lb' ub95=`k95ub' lb80=`k80lb' ub80=`k80ub'"
local worst = 0
local nr = rowsof(G95)
local nbad = 0
forvalues r = 1/`nr' {
  if `k95lb'==0 | `k80lb'==0 continue
  local h95 = (G95[`r',`k95ub'] - G95[`r',`k95lb'])/2
  local h80 = (G80[`r',`k80ub'] - G80[`r',`k80lb'])/2
  if missing(`h95') | missing(`h80') continue
  if `h95' == 0 continue
  local rat = `h80'/`h95'
  local d = abs(`rat' - t80/t95)/(t80/t95)
  if `d' > `worst' local worst = `d'
  if `d' > 1e-6 local ++nbad
  if `r' <= 3 dis "@@GROW r=`r' h95=" %12.6g `h95' " h80=" %12.6g `h80' " rat=" %12.8f `rat' " reldev=" %12.4e `d'
}
dis "@@GWORST worst_reldev=" %12.4e `worst' " nbad=`nbad' of `nr'"

*=========== E: main results / ES rows under level(95) vs level(80) ========
use "dat.dta", clear
qui wooldid y uid t ttre, att robust esprelength(2) espostlength(3) level(95) jointtests
mat M95 = e(wdidmainresults)
scalar dfm = e(df_r)
qui wooldid y uid t ttre, att robust esprelength(2) espostlength(3) level(80) jointtests
mat M80 = e(wdidmainresults)
local m95 : colfullnames M95
local m80 : colfullnames M80
local r95 : rowfullnames M95
dis "@@MN95 [`m95']"
dis "@@MN80 [`m80']"
dis "@@MR [`r95']"
scalar tm95 = invttail(dfm, .025)
scalar tm80 = invttail(dfm, .10)
dis "@@MDF df=" dfm " target_ratio=" tm80/tm95
mat list M95, format(%12.6g)
mat list M80, format(%12.6g)

*=========== F: non-integer level column naming ============================
use "dat.dta", clear
qui wooldid y uid t ttre, att robust level(92.5)
mat Q = e(wdidmainresults)
local qn : colfullnames Q
dis "@@QN [`qn'] elevel=" e(level)
qui wooldid y uid t ttre, att robust gteffects level(99.91)
mat Q2 = e(gteffects)
local q2n : colfullnames Q2
dis "@@Q2N [`q2n'] elevel=" e(level)
qui wooldid y uid t ttre, att robust level(10.01)
local q3n : colfullnames e(wdidmainresults)
dis "@@Q3N [`q3n'] elevel=" e(level)

dis "S9_DONE"
