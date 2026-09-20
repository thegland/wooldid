* p16: bayes across the option surfaces it is allowed on + reproducibility.
*  - reproducibility at a common rseed, and the reps()-prefix property
*  - event study (esprelength/espostlength)
*  - subgroups() + pairwise contrasts
*  - contreat() average marginal effects
*  - jointtests (fixed-metric bootstrap rank test)
* Each arm reports the CENTRING diagnostic (mean(draws)-b)/sd(draws): the F001 defect
* showed up there at 3.5-9.4.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

clear
set obs 30
gen int st = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2005 if st <= 8
replace onset = 2008 if st > 8 & st <= 16
replace onset = 2011 if st > 16 & st <= 22
gen byte sg = mod(st,2)
gen double dose0 = runiform()*2
expand 12
bysort st: gen int year = 2000 + _n
gen byte post = (year >= onset) & !missing(onset)
gen double dose = post*dose0
gen double y = a_i + 0.03*(year-2000) + 1.0*post + 0.4*dose + rnormal()*0.5
gen double pop = exp(rnormal()*0.4)*100
save p16_panel, replace

prog def ctr
  syntax , tag(string)
  matrix BB = e(b)
  matrix BD = e(bootdist)
  local K = colsof(BB)
  local nm : colnames BB
  mata {
    B = st_matrix("BB"); D = st_matrix("BD")
    for (k=1; k<=cols(B); k++) {
      d = select(D[.,k], D[.,k] :< .)
      if (rows(d) > 2 & sqrt(variance(d)) > 0) printf("`tag' col %2.0f  b=%12.6f  sd=%10.6f  (mean-b)/sd = %8.4f\n", k, B[1,k], sqrt(variance(d)), (mean(d)-B[1,k])/sqrt(variance(d)))
      else printf("`tag' col %2.0f  b=%12.6f  NO USABLE DRAWS\n", k, B[1,k])
    }
  }
end

dis "=== R1 reproducibility: same rseed twice ==="
use p16_panel, clear
qui wooldid y st year onset, att cluster(st) bootstrap(bayes) reps(50) rseed(4242) dots(0)
matrix D1 = e(bootdist)
use p16_panel, clear
qui wooldid y st year onset, att cluster(st) bootstrap(bayes) reps(50) rseed(4242) dots(0)
matrix D2 = e(bootdist)
mata: printf("R1 same-seed reruns bitwise identical: %f\n", (st_matrix("D1")==st_matrix("D2")))
dis "=== R2 reps() prefix property: first 50 of reps(120) == reps(50) ==="
use p16_panel, clear
qui wooldid y st year onset, att cluster(st) bootstrap(bayes) reps(120) rseed(4242) dots(0)
matrix D3 = e(bootdist)
mata: printf("R2 prefix identical: %f\n", (st_matrix("D3")[1..50,.]==st_matrix("D1")))

dis "=== S1 event study ==="
use p16_panel, clear
cap noi wooldid y st year onset, att cluster(st) esprelength(3) espostlength(3) bootstrap(bayes) reps(120) rseed(77) dots(0)
dis "S1_RC=" _rc
if _rc==0 ctr , tag(S1)

dis "=== S2 subgroups + jointtests ==="
use p16_panel, clear
cap noi wooldid y st year onset, att cluster(st) subgroups(sg) jointtests esprelength(2) espostlength(2) bootstrap(bayes) reps(120) rseed(77) dots(0)
dis "S2_RC=" _rc
if _rc==0 ctr , tag(S2)

dis "=== S3 contreat (AME) ==="
use p16_panel, clear
cap noi wooldid y st year onset, att cluster(st) contreat(dose) bootstrap(bayes) reps(120) rseed(77) dots(0)
dis "S3_RC=" _rc
if _rc==0 ctr , tag(S3)

dis "=== S4 weighted + size, centring (F001 regression check) ==="
use p16_panel, clear
cap noi wooldid y st year onset [aw=pop], att cluster(st) bootstrap(bayes) bootalpha(size) reps(200) rseed(77) dots(0)
dis "S4_RC=" _rc
if _rc==0 ctr , tag(S4)
use p16_panel, clear
cap noi wooldid y st year onset [aw=pop], att cluster(st) bootstrap(bayes) bootalpha(equal) reps(200) rseed(77) dots(0)
dis "S5_RC=" _rc
if _rc==0 ctr , tag(S5)
dis "PROBE_DONE"
