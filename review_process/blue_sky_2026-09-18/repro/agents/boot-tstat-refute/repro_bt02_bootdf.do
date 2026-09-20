*=====================================================================
* BT-02 VERDICT REPRO (skeptic's own; extends the original to 3 markout sources)
* Claim: bootdf(clusters)'s G is counted on (user if) & !missing(cluster), with NO
*   markout against the estimation sample, so clusters that contribute no row to
*   e(sample) are still in the degrees of freedom AND are still drawn.
* Look for: A_GAP / C_GAP / E_GAP nonzero while B_GAP / D_GAP are 0.
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 606060
set sortseed 606060
clear
set obs 26
gen int state = _n
gen int onset = .
qui replace onset = 2002 if state <= 6
qui replace onset = 2003 if state >= 7  & state <= 12
qui replace onset = 2004 if state >= 13 & state <= 18
gen double a_i = rnormal()
expand 3
bys state: gen int unit = _n + 100*state
gen double u_i = rnormal()*0.5
expand 8
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t-onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset
gen double xctl = rnormal() + 0.3*y
tempfile base
save `base', replace

cap prog drop gap
prog def gap
  syntax , tag(string)
  dis "`tag'_NCLUST=" e(N_clust) "  BOOTDF=" e(boot_df) "  DFR=" e(df_r) "  N=" e(N) ///
      "  b_att=" %14.10f e(b)[1,1]
  dis "`tag'_GAP=" e(boot_df) - (e(N_clust)-1)
end

dis "### A: outcome suppressed (missing) in 4 states, rows still in the data"
use `base', clear
qui replace y = . if inlist(state, 3, 9, 15, 21)
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) bootstrap(standard) reps(299) dots(0)
dis "A_RC=" _rc
gap , tag(A)
dis "### B: the SAME estimation sample, those 4 states physically dropped"
use `base', clear
qui drop if inlist(state, 3, 9, 15, 21)
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) bootstrap(standard) reps(299) dots(0)
dis "B_RC=" _rc
gap , tag(B)
dis "### C: a control variable missing in 3 states (markout via controls())"
use `base', clear
qui replace xctl = . if inlist(state, 5, 11, 17)
cap noi wooldid y state t ttre , cluster(state) att controls(xctl) esprelength(2) espostlength(2) bootstrap(standard) reps(299) dots(0)
dis "C_RC=" _rc
gap , tag(C)
dis "### D: same estimation sample as C, 3 states dropped"
use `base', clear
qui drop if inlist(state, 5, 11, 17)
cap noi wooldid y state t ttre , cluster(state) att controls(xctl) esprelength(2) espostlength(2) bootstrap(standard) reps(299) dots(0)
dis "D_RC=" _rc
gap , tag(D)
dis "### E: analysis weight 0 in 3 states (rows present, contribute nothing)"
use `base', clear
gen double wt = 1
qui replace wt = 0 if inlist(state, 4, 10, 16)
cap noi wooldid y state t ttre [aw=wt] , cluster(state) att esprelength(2) espostlength(2) bootstrap(standard) reps(299) dots(0)
dis "E_RC=" _rc
gap , tag(E)
dis "BT02_DONE"
