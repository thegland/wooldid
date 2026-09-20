clear all
set more off
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
use res04a.dta, clear
append using res04b.dta
dis "NROWS=" _N
tab arm
* reshape to wide by ds
keep ds arm p sd nge nfail
reshape wide p sd nge nfail, i(ds) j(arm) string
dis "NDS=" _N
foreach A in CTRLY CTRLX FEWAVE CCC {
  gen double sdratio_`A' = sd`A'/sdBASE
  gen double pdiff_`A'   = p`A' - pBASE
  gen byte flip_`A' = ((pBASE<=.05) != (p`A'<=.05)) if !missing(pBASE, p`A')
}
sum sdratio_* , format
sum pdiff_*, format
foreach A in CTRLY CTRLX FEWAVE CCC {
  qui count if flip_`A'==1
  dis "FLIP5PCT_`A' = " r(N) " of " _N
}
qui count if pBASE<=.05
dis "REJ_BASE = " r(N) " of " _N
foreach A in CTRLY CTRLX FEWAVE CCC {
  qui count if p`A'<=.05
  dis "REJ_`A' = " r(N) " of " _N
}
* rank of observed among the 90 (nge = #draws >= b_obs); uniform 1..90 under exactness
foreach A in BASE CTRLY CTRLX FEWAVE CCC {
  qui sum nge`A'
  dis "RANK_`A': mean=" %8.3f r(mean) " sd=" %8.3f r(sd) " min=" r(min) " max=" r(max)
}
list ds pBASE pCTRLY pCTRLX pFEWAVE pCCC sdratio_CCC, noobs
save res04all.dta, replace
dis "PROBE08 DONE"
