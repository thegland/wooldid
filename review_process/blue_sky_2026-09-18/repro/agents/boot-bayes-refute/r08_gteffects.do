* r08 (own hunt): the safety gate protects att_it/att_i/att_itime.  But e(gteffects)
* is a table of PER-CELL effects, and with i = unit each cell sits wholly inside one
* cluster, so the Dirichlet multiplier is constant on it - the same geometry that makes
* the cell-weighted estimands collapse.  gteffects is available with plain `att', which
* needs NO safety(off).  Are the gteffects cell standard errors under bootstrap(bayes)
* understated relative to bootstrap(standard) on the same panel and seed, with no gate
* and no caution?  Same question for the event-study rows.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

set seed 424242
set sortseed 424242
clear
set obs 28
gen int st = _n
gen double a_i = rnormal()*1.2
gen double shk = rnormal()
gen int onset = .
replace onset = 2006 if st <=  7
replace onset = 2009 if st >  7 & st <= 14
replace onset = 2012 if st > 14 & st <= 20
expand 14
bysort st: gen int year = 2000 + _n
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.04*(year-2000) + (0.8 + 1.5*shk)*post + rnormal()*0.6
tempfile dd
qui save `dd'

foreach eng in analytic standard bayes {
  use `dd', clear
  if "`eng'" == "analytic" local opt
  else local opt bootstrap(`eng') reps(200) rseed(77777) dots(0)
  cap noi qui wooldid y st year onset, att cluster(st) gteffects espostlength(3) esprelength(2) `opt'
  dis "`eng'_RC=" _rc
  if _rc == 0 {
    matrix G = e(gteffects)
    if "`eng'" == "analytic" {
      dis "gteffects colnames: " "`: colfullnames G'"
      dis "gteffects nrows=" rowsof(G) " ncols=" colsof(G)
      matrix list G, format(%9.4f) nohalf
      matrix M = e(wdidmainresults)
      dis "main rownames: " "`: rowfullnames M'"
    }
    matrix `eng'G = G
    matrix `eng'M = e(wdidmainresults)
  }
}
* locate the se column of e(gteffects) by name
local gc : colfullnames analyticG
dis "GTCOLS: `gc'"
local sepos = 0
local j = 0
foreach c of local gc {
  local ++j
  if strpos(lower("`c'"), "se") > 0 & `sepos'==0 local sepos = `j'
}
dis "se column index = `sepos'"
if `sepos' > 0 {
  local nr = rowsof(analyticG)
  dis "row  analytic_se   standard_se   bayes_se   bayes/standard"
  local sa = 0
  local ss = 0
  local sb = 0
  forvalues r = 1/`nr' {
    local va = analyticG[`r',`sepos']
    local vs = standardG[`r',`sepos']
    local vb = bayesG[`r',`sepos']
    if `va' < . & `vs' < . & `vb' < . {
      dis "`r'  " %10.6f `va' "  " %10.6f `vs' "  " %10.6f `vb' "   ratio=" %8.4f `vb'/`vs'
      local sa = `sa' + `va'
      local ss = `ss' + `vs'
      local sb = `sb' + `vb'
    }
  }
  dis "MEAN gteffects se: analytic=" %10.6f `sa'/`nr' "  standard=" %10.6f `ss'/`nr' "  bayes=" %10.6f `sb'/`nr'
  dis "MEAN RATIO bayes/standard over gteffects cells = " %8.4f `sb'/`ss'
}
* main results table rows (att + event study)
local mr : rowfullnames analyticM
local nr = rowsof(analyticM)
local j = 0
foreach r of local mr {
  local ++j
  dis "MAIN `r'  analytic_se=" %10.6f analyticM[`j',2] "  std_se=" %10.6f standardM[`j',2] "  bayes_se=" %10.6f bayesM[`j',2] "  ratio=" %8.4f bayesM[`j',2]/standardM[`j',2]
}
dis "R08_DONE"
