* r06 shard a: coverage Monte Carlo for bb-03(b) - the only part of that finding
* the prober marked uncertain.  30 clusters x 10 periods, i = st, cluster(st),
* TRUE att = 1, homogeneous effects, reps(79), nominal 95% - the design bb-03 states.
* Records the REPORTED percentile and normal-theory intervals for att under
* bootstrap(standard) and bootstrap(bayes), PAIRED on the same panel and seed.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid
tempname fh
file open `fh' using "r06_a.csv", write replace
file write `fh' "p,eng,b,se,pct_lo,pct_hi,wald_lo,wald_hi,analse" _n

forvalues p = 1/20 {
  set seed `=800000 + `p'*13'
  set sortseed `=800000 + `p'*13'
  clear
  set obs 30
  gen int st = _n
  gen double a_i = rnormal()
  gen int onset = .
  replace onset = 2004 if st <=  8
  replace onset = 2006 if st >  8 & st <= 15
  replace onset = 2008 if st > 15 & st <= 21
  expand 10
  bysort st: gen int year = 2000 + _n
  gen byte post = (year >= onset) & !missing(onset)
  gen double y = a_i + 0.03*(year-2000) + 1.0*post + rnormal()*0.7
  tempfile dd
  qui save `dd'
  cap noi qui wooldid y st year onset, att cluster(st)
  local ase = .
  if _rc == 0 {
    matrix A0 = e(wdidmainresults)
    local ase = A0[1,2]
  }
  foreach eng in standard bayes {
    use `dd', clear
    cap noi qui wooldid y st year onset, att cluster(st) bootstrap(`eng') reps(79) rseed(`=5000+`p'') dots(0)
    local rc = _rc
    dis "P`p'_`eng'_RC=`rc'"
    if `rc' == 0 {
      matrix A  = e(wdidmainresults)
      matrix CP = e(boot_ci_percentile)
      matrix CN = e(boot_ci_normal)
      file write `fh' "`p',`eng'," (A[1,1]) "," (A[1,2]) "," (CP[1,1]) "," (CP[2,1]) "," (CN[1,1]) "," (CN[2,1]) "," (`ase') _n
      if `p' == 1 & "`eng'"=="standard" {
        dis "rownames A: " "`: rowfullnames A'"
        dis "CP rows: " "`: rowfullnames CP'" "  cols: " "`: colfullnames CP'"
      }
    }
  }
}
file close `fh'
dis "R06_a_DONE"
