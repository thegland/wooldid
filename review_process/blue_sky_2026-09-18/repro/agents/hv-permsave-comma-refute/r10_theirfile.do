* r10_theirfile.do - read-only inspection of the OTHER agent's bootsave artifacts, to
* check whether their reported "observed statistic in row 0 of jointW_att = .2704913"
* is an instance of N1 (a coefficient written where the Wald statistic belongs).
clear all
set more off
set processors 2
local d "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma"
foreach f in bdraws b1 b2 {
  capture noisily use "`d'/`f'.dta", clear
  dis "FILE `f' useRC=" _rc
  if _rc == 0 {
    qui ds
    dis "`f'_COLS=" r(varlist)
    capture confirm variable jointW_att
    if _rc == 0 {
      qui su jointW_att if _rep == 0, meanonly
      dis "`f'_jointW_att_row0=" %21x r(mean) " = " r(mean)
      qui su jointW_att if _rep > 0
      dis "`f'_jointW_att_reps n=" r(N) " min=" r(min) " max=" r(max)
      qui ds _rep jointW_* se_* w_*, not
      foreach v of local r(varlist) {
        qui su `v' if _rep == 0, meanonly
        dis "  `f' coef `v' row0=" %21x r(mean)
      }
    }
    capture confirm variable advW
    if _rc == 0 {
      qui su advW if _rep == 0, meanonly
      dis "`f'_advW_row0=" %21x r(mean) " = " r(mean)
      qui su advW if _rep > 0
      dis "`f'_advW_reps n=" r(N) " min=" r(min) " max=" r(max)
    }
  }
}
dis "ALLDONE"
