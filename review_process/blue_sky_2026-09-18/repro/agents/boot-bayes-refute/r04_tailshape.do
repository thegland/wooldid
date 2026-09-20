* r04: bb-03(b).  The claim is that at cluster()==i the bayes PERCENTILE interval for
* att sits at a different LOCATION from the bootstrap(standard) one at equal width, and
* under-covers.  Diagnose the replication distributions directly at reps(999), which
* needs no Monte Carlo, and check whether the gap is a low-reps artifact (their MC ran
* reps(79), where the 2.5% endpoint is essentially the 2nd order statistic).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes-refute"
which wooldid

forvalues p = 1/3 {
  set seed `=600000 + `p'*77'
  set sortseed `=600000 + `p'*77'
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

  foreach eng in standard bayes {
    foreach R in 79 999 {
      use `dd', clear
      cap noi qui wooldid y st year onset, att cluster(st) bootstrap(`eng') reps(`R') rseed(9090) dots(0)
      local rc = _rc
      dis "P`p'_`eng'_`R'_RC=`rc'"
      if `rc' == 0 {
        matrix RR = e(wdidmainresults)
        if `p'==1 & "`eng'"=="standard" & `R'==79 dis "colnames wdidmainresults: " "`: colfullnames RR'"
        local b = RR[1,1]
        matrix BD = e(bootdist)
        local ncol = colsof(BD)
        if `p'==1 & "`eng'"=="standard" & `R'==79 dis "bootdist colnames: " "`: colfullnames BD'"
        preserve
        clear
        qui svmat double BD, names(col)
        qui ds
        local first : word 1 of `r(varlist)'
        qui drop if missing(`first')
        qui count
        local n = r(N)
        qui su `first', detail
        local mn = r(mean)
        local sd = r(sd)
        local md = r(p50)
        local sk = r(skewness)
        qui _pctile `first', p(2.5 97.5)
        local lo = r(r1)
        local hi = r(r2)
        dis "P`p' `eng' R=`R' n=`n' b=" %10.6f `b' " repmean=" %10.6f `mn' " repmed=" %10.6f `md' " sd=" %10.6f `sd' " skew=" %8.4f `sk'
        dis "P`p' `eng' R=`R' CI=[" %10.6f `lo' "," %10.6f `hi' "]  width=" %9.6f `hi'-`lo' " centre-b=" %10.6f (`lo'+`hi')/2-`b' " (rep mean-b)/sd=" %8.4f (`mn'-`b')/`sd'
        dis "P`p' `eng' R=`R' covers1=" (`lo' <= 1 & 1 <= `hi')
        restore
      }
    }
  }
}
dis "R04_DONE"
