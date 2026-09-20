clear all
set more off
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
use res04all.dta, clear
foreach A in CTRLY CTRLX FEWAVE CCC {
  qui sum sdratio_`A', det
  dis "SDRATIO_`A': p10=" %8.4f r(p10) " p25=" %8.4f r(p25) " med=" %8.4f r(p50) " p75=" %8.4f r(p75) " p90=" %8.4f r(p90) " max=" %10.4f r(max)
  qui count if sdratio_`A' < 1
  dis "SDRATIO_`A': narrower than BASE in " r(N) " of " _N
}
foreach A in BASE CTRLY CTRLX FEWAVE CCC {
  qui sum p`A', det
  dis "P_`A': min=" %8.4f r(min) " p25=" %8.4f r(p25) " med=" %8.4f r(p50) " p75=" %8.4f r(p75) " max=" %8.4f r(max) " mean=" %8.4f r(mean)
}
qui count if pCCC > pBASE
dis "CCC p LARGER than BASE in " r(N) " of " _N
qui count if sdratio_CCC > 2
dis "CCC nullSD more than 2x BASE in " r(N) " of " _N
qui count if sdratio_CCC > 5
dis "CCC nullSD more than 5x BASE in " r(N) " of " _N
dis "PROBE09 DONE"
