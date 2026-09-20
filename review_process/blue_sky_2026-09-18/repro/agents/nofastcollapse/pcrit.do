clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"
local OO "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/out"
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/lib_dump.do"
local ES att esprelength(3) espostlength(3) jointtests

* ===== (1) NOISE FLOOR WITH THE CONTROL PRESENT.  Must be EXACTLY 0. =====
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/nfc_A1.txt")
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/nfc_A2.txt")
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
wddump , file("`OO'/nfc_B1.txt")
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
wddump , file("`OO'/nfc_B2.txt")
* and the REVERSED order (B first, then A) to rule out any within-session order effect
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/nfc_Arev.txt")

* ===== (2) FRAMING: is a plain ROW PERMUTATION (collapse OFF in both) as big? =====
use "`DD'/dCc.dta", clear
set seed 4242
gen double _s = runiform()
sort _s
drop _s
save "`DD'/dCc_shuf.dta", replace
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
wddump , file("`OO'/y01_orig_nofc.txt")
qui use "`DD'/dCc_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar) nofastcollapse
wddump , file("`OO'/y01_shuf_nofc.txt")
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/y01_orig_fc.txt")
qui use "`DD'/dCc_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' controls(xdollar)
wddump , file("`OO'/y01_shuf_fc.txt")
* the same pair WITHOUT any control, as the reference the help's calibration line describes
qui use "`DD'/dCc.dta", clear
qui wooldid y cohort t ttre, `ES' nofastcollapse
wddump , file("`OO'/y02_orig_nofc.txt")
qui use "`DD'/dCc_shuf.dta", clear
qui wooldid y cohort t ttre, `ES' nofastcollapse
wddump , file("`OO'/y02_shuf_nofc.txt")

* ===== (3) the two control sets dropped earlier on a syntax error =====
local BASE wooldid y cohort t ttre
wdcell , dta("`DD'/dCc.dta") tag(y03_ccc)    base(`BASE') out("`OO'") opts(`ES' coarsecohortcontrols(c.xdollar))
wdcell , dta("`DD'/dCc.dta") tag(y04_cccabs) base(`BASE') out("`OO'") opts(`ES' ccc_absorb(xdollar))
wdcell , dta("`DD'/dCc.dta") tag(y05_inter)  base(`BASE') out("`OO'") opts(`ES' interactivecontrols(xdollar))

*==== CRIT: can the residue flip the adversarial 1e-12 rank truncation? ====
use "`DD'/dCc.dta", clear
gen double xs4 = 10000*x1 + 50000
save "`DD'/dCs.dta", replace
local BASE2 wooldid y cohort t ttre
wdcell , dta("`DD'/dCs.dta") tag(x05_advc) base(`BASE2') out("`OO'") opts(att esprelength(4) espostlength(4) cluster(cohort) esfixedbaseperiod adversarial jointtests controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x06_band) base(`BASE2') out("`OO'") opts(att att_i esprelength(4) espostlength(4) cluster(cohort) esfixedbaseperiod adversarial jointtests simultaneousbands pretrendpower(0.5) controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x07_sgc)  base(`BASE2') out("`OO'") opts(att esprelength(3) espostlength(3) jointtests subgroups(sgc) controls(xs4))
wdcell , dta("`DD'/dCs.dta") tag(x05_noc)  base(`BASE2') out("`OO'") opts(att esprelength(4) espostlength(4) cluster(cohort) esfixedbaseperiod adversarial jointtests)

*==== CRIT: rows-per-cell ladder (severity upgrade criterion) ====
foreach u in 4 8 32 128 512 {
  clear
  set seed 20260917
  set sortseed 20260917
  qui set obs `=4*`u''
  qui gen long unit   = _n
  qui gen byte cohort = mod(_n-1,4)+1
  qui gen double ttre = .
  qui replace ttre = 2004 if cohort==1
  qui replace ttre = 2006 if cohort==2
  qui replace ttre = 2008 if cohort==3
  qui gen double a_i = rnormal()
  qui expand 10
  qui bysort unit: gen int t = 2000 + _n - 1
  qui gen double x1 = rnormal()
  qui gen double y = a_i + 0.4*(t-2000) + 0.5*x1 + 1.2*(ttre<. & t>=ttre) + rnormal()*0.3
  qui save "`DD'/cs`u'.dta", replace
  wdcell , dta("`DD'/cs`u'.dta") tag(cs`u') base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(3) espostlength(3) jointtests controls(x1))
  wdcell , dta("`DD'/cs`u'.dta") tag(cs`u'_noc) base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(3) espostlength(3) jointtests)
}
dis "PCRIT_DONE"
