* probe04_multi.do -- how BIG is the distortion, across many datasets?
* Sharp null (no treatment effect), RANDOMISED onset assignment each dataset (which is
* what makes the randomisation test's exactness testable), 6 states, M=90 enumerated.
* Arms: BASE; CTRL_Y (controls(cohbaseY) where cohbaseY = the cohort's own mean
* pre-2005 outcome -- the "control for the cohort's baseline level" a real user writes);
* CTRL_X (controls(cohidx), an arbitrary cohort-level index); FEWAVE (fe(early)).
* Every covariate is constant within the OBSERVED cohort, hence inert for the estimate.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-refute"
which wooldid
local LO = $LO
local HI = $HI
local OUT "$OUT"

tempname M
postfile `M' int ds str8 arm double p double sd double bobs int nge int nfail using `OUT', replace

forvalues d = `LO'/`HI' {
  clear
  set seed `=90000 + `d''
  set sortseed `=90000 + `d''
  set obs 6
  gen int state = _n
  gen double u = runiform()
  sort u
  gen byte coh = ceil(_n/2)
  drop u
  sort state
  gen int onset = cond(coh==1, 2005, cond(coh==2, 2009, .))
  gen int t0 = 2000 + mod(state,3)
  gen int t1 = 2013 + mod(state,2)
  gen double a_s = rnormal()*1.2
  gen double g_s = rnormal()*0.03
  expand t1-t0+1
  bys state: gen int t = t0 + _n - 1
  drop t0 t1
  gen double y = a_s + g_s*(t-2006) + 0.5*rnormal()
  * SHARP NULL: y does not depend on onset at all.
  * cohort-level covariates, both constant within coh:
  bys coh: egen double cohbaseY = mean(cond(t<2005, y, .))
  gen double cohidx = cond(coh==1, 0.80, cond(coh==2, 1.45, 2.10))
  gen byte early = (coh==1)
  compress
  tempfile pp
  save `pp'
  foreach A in BASE CTRLY CTRLX FEWAVE CCC {
    use `pp', clear
    local o
    if "`A'"=="CTRLY"  local o controls(cohbaseY)
    if "`A'"=="CTRLX"  local o controls(cohidx)
    if "`A'"=="FEWAVE" local o fe(early)
    if "`A'"=="CCC"    local o coarsecohortcontrols(c.cohbaseY)
    cap noi qui wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts `o'
    local rc = _rc
    if `rc'==0 {
      mat B = e(b)
      mat D = e(permdist_b)
      local bo = B[1,1]
      local nn = rowsof(D)
      local cnt = 0
      forvalues r = 1/`nn' {
        if D[`r',1] >= `bo' local cnt = `cnt' + 1
      }
      post `M' (`d') ("`A'") (e(perm_p_att)) (e(perm_sd_att)) (`bo') (`cnt') (e(perm_failed))
    }
    else post `M' (`d') ("`A'") (.) (.) (.) (.) (-`rc')
    dis "DS`d'_`A'_RC=`rc'"
  }
}
postclose `M'
use `OUT', clear
dis "ROWS=" _N
list in 1/12
