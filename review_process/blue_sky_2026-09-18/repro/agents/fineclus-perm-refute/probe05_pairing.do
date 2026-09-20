* probe05_pairing.do -- the PAIRED per-assignment comparison, from the saved draw files.
* No wooldid call: pure arithmetic on e(permdist_b) already captured in probe03.
clear all
set more off
set processors 2
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
foreach P in A B {
  use D_`P'BASE.dta, clear
  foreach A in CTRL FEWAVE FECOH CCC {
    merge 1:1 _rep using D_`P'`A'.dta, nogen
  }
  local bb b_`P'BASE
  qui sum `bb'
  dis "PANEL`P'_SD_BASE = " %14.10f r(sd)
  foreach A in CTRL FEWAVE FECOH CCC {
    local ba b_`P'`A'
    qui gen double d_`A' = `ba' - `bb'
    qui count if abs(d_`A') > 1e-10
    local nd = r(N)
    qui sum d_`A'
    dis "PANEL`P'_`A': NDIFF=`nd' of " _N "  meandiff=" %12.8f r(mean) "  maxabsdiff=" %12.8f max(abs(r(min)),abs(r(max)))
    qui sum `ba'
    dis "PANEL`P'_`A': SD=" %14.10f r(sd)
  }
  * which row is the observed (identity) assignment, and is the diff zero there?
  gen double absd = abs(d_CTRL)
  sort absd
  list _rep `bb' d_CTRL d_FEWAVE d_FECOH d_CCC in 1/3, noobs
  dis "PANEL`P'_MIN_ABS_D_CTRL = " %21x absd[1]
  drop absd
  save PAIRED_`P'.dta, replace
}
dis "PROBE05 DONE"
