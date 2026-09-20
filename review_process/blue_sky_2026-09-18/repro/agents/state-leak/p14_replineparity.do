* p14: REPLICATION-LINE PARITY.
* Under bootstrap() the shell re-parses the command line (`qui wooldid_validate `0''
* at ado:312) and rebuilds a replication option string BY HAND (`wdb_ro', ado:313-365).
* If an option that shapes the POINT ESTIMATES is missing from that string, every
* replication estimates a different statistic than the observed run -- silently.
* Two observables per spec:
*   (a) e(bias_boot)/sd  -- a screen, not a verdict (G=12 gives genuine small-G bias)
*   (b) the DISCRIMINATING test: is e(b_boot) closer to the observed b WITH the option,
*       or to the observed b WITHOUT it?  Both are printed for every spec.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

* a genuinely time-varying control (unit-constant controls are absorbed by the unit FE)
cap prog drop wdgenx
prog def wdgenx
  wdgen
  qui gen double xvar = 0.4*rnormal() + 0.1*year
  qui replace y = y + 0.25*xvar
end

* print b for one run
cap prog drop wdb1
prog def wdb1
  args tag
  capture confirm matrix e(b)
  if _rc != 0 {
    dis "B1|`tag'|NO_b"
    exit
  }
  tempname B
  matrix `B' = e(b)
  local cn : colfullnames `B'
  forvalues j = 1/`=colsof(`B')' {
    local nm : word `j' of `cn'
    dis "B1|`tag'|`nm'=" %18.12g `B'[1,`j']
  }
end

cap prog drop wdbias
prog def wdbias
  args tag
  capture confirm matrix e(b_boot)
  if _rc != 0 {
    dis "BIAS|`tag'|NO_b_boot"
    exit
  }
  tempname B BB BI VV
  matrix `B'  = e(b)
  matrix `BB' = e(b_boot)
  matrix `BI' = e(bias_boot)
  matrix `VV' = e(V)
  local cn : colfullnames `B'
  dis "BIAS|`tag'|completed=" e(boot_completed) "|failed=" e(boot_failed)
  forvalues j = 1/`=colsof(`B')' {
    local nm : word `j' of `cn'
    local bo = `B'[1,`j']
    local bm = `BB'[1,`j']
    local bi = `BI'[1,`j']
    local sd = sqrt(`VV'[`j',`j'])
    local rat = .
    if `sd' > 0 & `sd' < . local rat = `bi'/`sd'
    dis "BIAS|`tag'|`nm'|b=" %18.12g `bo' "|bmean=" %18.12g `bm' "|sd=" %12.6g `sd' "|bias_over_sd=" %9.4f `rat'
  }
end

local R = 60

* ============ 1 plain att (the control: no option to drop)
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att
dis "S1a_RC=" _rc
wdb1 01_plain_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S1_RC=" _rc
wdbias 01_plain

* ============ 2 fe()
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att fe(sg)
dis "S2a_RC=" _rc
wdb1 02_fe_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att fe(sg) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S2_RC=" _rc
wdbias 02_fe

* ============ 3 controls() -- time-varying
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att controls(c.xvar)
dis "S3a_RC=" _rc
wdb1 03_controls_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att controls(c.xvar) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S3_RC=" _rc
wdbias 03_controls

* ============ 4 timetrends
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att timetrends
dis "S4a_RC=" _rc
wdb1 04_timetrends_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att timetrends bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S4_RC=" _rc
wdbias 04_timetrends

* ============ 5 esfixedbaseperiod
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att esfixedbaseperiod esprelength(3) espostlength(3)
dis "S5a_RC=" _rc
wdb1 05_esfixedbase_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att esfixedbaseperiod esprelength(3) espostlength(3) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S5_RC=" _rc
wdbias 05_esfixedbase

* ============ 6 anticipation()
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att anticipation(1)
dis "S6a_RC=" _rc
wdb1 06_anticipation_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att anticipation(1) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S6_RC=" _rc
wdbias 06_anticipation

* ============ 7 oostreatedcontrols()
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att oostreatedcontrols(2)
dis "S7a_RC=" _rc
wdb1 07_oostreated_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att oostreatedcontrols(2) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S7_RC=" _rc
wdbias 07_oostreated

* ============ 8 subgroups()
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att subgroups(sg)
dis "S8a_RC=" _rc
wdb1 08_subgroups_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att subgroups(sg) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S8_RC=" _rc
wdbias 08_subgroups

* ============ 9 customweightmultiplier()
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att customweightmultiplier(wt)
dis "S9a_RC=" _rc
wdb1 09_cwm_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att customweightmultiplier(wt) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S9_RC=" _rc
wdbias 09_cwm

* ============ 10 att_i (cell-weighted; physical re-id leg)
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att_i
dis "S10a_RC=" _rc
wdb1 10_att_i_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att_i bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S10_RC=" _rc
wdbias 10_att_i

* ============ 11 contreat()
wdgenx
qui gen double dose = 0
qui replace dose = 0.5 + 0.5*mod(unit,3) if !missing(ttre) & year >= ttre
qui replace y = y + 0.3*dose
capture noisily wooldid y unit year ttre, cluster(state) att contreat(dose)
dis "S11a_RC=" _rc
wdb1 11_contreat_noboot
wdgenx
qui gen double dose = 0
qui replace dose = 0.5 + 0.5*mod(unit,3) if !missing(ttre) & year >= ttre
qui replace y = y + 0.3*dose
capture noisily wooldid y unit year ttre, cluster(state) att contreat(dose) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S11_RC=" _rc
wdbias 11_contreat

* ============ 12 interactivecontrols() -- time-varying
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att interactivecontrols(xvar)
dis "S12a_RC=" _rc
wdb1 12_interactive_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att interactivecontrols(xvar) bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S12_RC=" _rc
wdbias 12_interactive

* ============ 13 aweight
wdgenx
capture noisily wooldid y unit year ttre [aw=wt], cluster(state) att
dis "S13a_RC=" _rc
wdb1 13_aweight_noboot
wdgenx
capture noisily wooldid y unit year ttre [aw=wt], cluster(state) att bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S13_RC=" _rc
wdbias 13_aweight

* ============ 14 an if-restricted sample
wdgenx
capture noisily wooldid y unit year ttre if state <= 10, cluster(state) att
dis "S14a_RC=" _rc
wdb1 14_if_noboot
wdgenx
capture noisily wooldid y unit year ttre if state <= 10, cluster(state) att bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S14_RC=" _rc
wdbias 14_if

* ============ 15 semielasticity (the rep line re-emits it only when validate
* normalised it to the string "eydx" -- ado:350)
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att semielasticity
dis "S15a_RC=" _rc
wdb1 15_semielas_noboot
wdgenx
capture noisily wooldid y unit year ttre, cluster(state) att semielasticity bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S15_RC=" _rc
wdbias 15_semielas

* ============ 16 poisson
wdgenx
qui replace y = exp(y/3)
capture noisily wooldid y unit year ttre, cluster(state) att poisson
dis "S16a_RC=" _rc
wdb1 16_poisson_noboot
wdgenx
qui replace y = exp(y/3)
capture noisily wooldid y unit year ttre, cluster(state) att poisson bootstrap(standard) reps(`R') dots(0) rseed(7)
dis "S16_RC=" _rc
wdbias 16_poisson

dis "P14_DONE"
