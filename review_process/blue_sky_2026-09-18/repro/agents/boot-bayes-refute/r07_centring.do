* r07 (own hunt): the Dirichlet CENTRING fixes are the freshest code in this engine
* (the bootalpha(equal) split was changed 2026-09-06, the bootalpha(size) divisor
* 2026-09-13 - ado:14288-14310).  Both fixes are about where the replication
* distribution SITS relative to the point estimate on UNBALANCED, WEIGHTED data.
* Measure (mean of draws - b) / sd of draws for att, att_it and att_i under
*   bootstrap(standard), bootstrap(bayes)+bootalpha(equal), bootstrap(bayes)+bootalpha(size)
* on a deliberately unbalanced, aweighted panel.  A displacement above ~0.3 bootstrap
* sd is a centring defect; the ado comment says the old bug gave 9.36.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

set seed 191919
set sortseed 191919
clear
set obs 26
gen int st = _n
gen double a_i = rnormal()*1.2
gen int onset = .
replace onset = 2006 if st <=  7
replace onset = 2009 if st >  7 & st <= 13
replace onset = 2012 if st > 13 & st <= 18
* strongly unbalanced panel lengths: some states enter late
gen int first = 2000 + floor(runiform()*5)
expand 14
bysort st: gen int year = 2000 + _n
drop if year < first
* strongly unbalanced weights (population), correlated with onset
gen double pop = exp(rnormal()*1.1) * (1 + 2*(onset < 2010 & !missing(onset)))
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.04*(year-2000) + 0.8*post + rnormal()*0.6
qui count
dis "ROWS=" r(N)
qui gstats sum pop, meanonly
tempfile dd
qui save `dd'
qui bysort st: gen byte tg = _n==1
qui su pop, detail
dis "pop p10=" %8.4f r(p10) " p90=" %8.4f r(p90)

foreach spec in "standard" "bayes_equal" "bayes_size" {
  use `dd', clear
  if "`spec'" == "standard"    local opt bootstrap(standard)
  if "`spec'" == "bayes_equal" local opt bootstrap(bayes) bootalpha(equal) safety(off)
  if "`spec'" == "bayes_size"  local opt bootstrap(bayes) bootalpha(size) safety(off)
  cap noi qui wooldid y st year onset [aw=pop], att att_it att_i cluster(st) `opt' reps(300) rseed(4321) dots(0)
  dis "`spec'_RC=" _rc
  if _rc == 0 {
    matrix A = e(wdidmainresults)
    dis "`spec' rownames: " "`: rowfullnames A'"
    matrix BD = e(bootdist)
    dis "`spec' bootdist cols: " "`: colfullnames BD'"
    local nr = rowsof(A)
    preserve
    clear
    qui svmat double BD, names(col)
    qui ds
    local vl `r(varlist)'
    local j = 0
    foreach v of local vl {
      local ++j
      if `j' > `nr' continue
      qui su `v'
      local m = r(mean)
      local s = r(sd)
      local n = r(N)
      local b = A[`j',1]
      dis "`spec' col`j'=`v'  n=`n'  b=" %10.6f `b' "  repmean=" %10.6f `m' "  sd=" %10.6f `s' "  (mean-b)/sd=" %9.4f (`m'-`b')/`s'
    }
    restore
  }
}
dis "R07_DONE"
