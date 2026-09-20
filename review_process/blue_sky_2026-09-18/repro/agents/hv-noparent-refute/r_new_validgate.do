* ============================================================================
* MY OWN (unseeded) probe: the generalizable class behind the seeded finding is
* not "inert options" but "VALUE-DOMAIN VALIDATION GATED BEHIND A PARENT".
* wooldid's own stated convention (sthlp:917, on contreatpoly()) is that a
* value "is validated whether or not contreat() is given". This probe sweeps
* every other option that takes a value and asks whether an out-of-domain or
* nonsensical value is caught when the option's consumer is not reached.
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r_new_validgate.do
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 606
set sortseed 606

which wooldid

clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = rnormal()
gen double dose0 = 0.5 + runiform()*2
expand 10
bysort unit: gen int t = _n
gen double y = 2 + ufe + 0.1*t + rnormal()
replace y = y + 0.8*(t>=ttre) if ttre<.
gen double dose = dose0*(t>=ttre & ttre<.)
gen double wgt = 1 + runiform()
tempfile D
save `D', replace

cap prog drop RUNIT
prog def RUNIT
  syntax , nm(string) opts(string asis)
  capture noisily wooldid y unit t ttre, `opts'
  dis "VG_`nm'_RC=" _rc
end

* --- esrelativeto(): help says (-int). A POSITIVE value, no esfixedbaseperiod --
use `D', clear
RUNIT, nm(esrel_pos_noparent) opts(att robust esrelativeto(3))
use `D', clear
RUNIT, nm(esrel_pos_withES) opts(att robust esrelativeto(3) esprelength(3) espostlength(2) esfixedbaseperiod)
use `D', clear
RUNIT, nm(esrel_neg_noES) opts(att robust esrelativeto(-3))

* --- esrelativeto() far outside the requested pre-window -----------------------
use `D', clear
RUNIT, nm(esrel_outofwindow) opts(att robust esrelativeto(-99) esprelength(2) espostlength(2) esfixedbaseperiod)

* --- pretrendpower(): a nonsense power target with no event-study window -------
use `D', clear
RUNIT, nm(ptp_nowindow) opts(att robust pretrendpower(0.8))
use `D', clear
RUNIT, nm(ptp_bad_nowindow) opts(att robust pretrendpower(5))
use `D', clear
RUNIT, nm(ptp_bad_window) opts(att robust pretrendpower(5) esprelength(3) espostlength(2))

* --- anticipation() / oostreatedcontrols() with no event-study request ---------
use `D', clear
RUNIT, nm(antic_nowindow) opts(att robust anticipation(3))
use `D', clear
RUNIT, nm(oos_nowindow) opts(att robust oostreatedcontrols(3))

* --- cfxlattice() without makecfxplots: documented CAUTION (ado:5307) ----------
use `D', clear
RUNIT, nm(cfxlat_noparent) opts(att robust contreat(dose) cfxlattice(dose0))

* --- customweightmultiplier(): help says >=0. Feed it negatives. ---------------
use `D', clear
replace wgt = -1 in 1/50
RUNIT, nm(cwm_negative) opts(att robust customweightmultiplier(wgt))

* --- regtol(): an absurd tolerance --------------------------------------------
use `D', clear
RUNIT, nm(regtol_huge) opts(att robust regtol(40))

dis "DONE_VG"
