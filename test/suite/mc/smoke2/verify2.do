* verify2.do (v2) - the one fact mc_bootse2's design A rests on: exact shift
* equivariance of the bootstrap surfaces under y -> y - tau*D, which is what makes
* the tau=0 run's level-independent p-values the tau=2 run's coverage indicators.
* Plus two cheap side checks (V1 windowed truth on B; V4 the r(198) seen on
* "att att_it" together with an event-study window).
clear all
set more off
set linesize 250
set graphics off
qui do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"
run "/Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/mc_bootse2_dgp.doh"

di "VER ===== V1b design B windowed truth (att alone) ====="
dgpB 1 1
qui replace y = tau*(t>=ttre & ttre<.)
qui wooldid y unit t ttre, att cluster(state) esprelength(2) espostlength(1) safety(off)
di "VER B windowed att   = " %20.16f el(e(wdidmainresults),1,1)
qui wooldid y unit t ttre, att_it cluster(state) esprelength(2) espostlength(1) safety(off)
di "VER B windowed att_it= " %20.16f el(e(wdidmainresults),1,1)
qui wooldid y unit t ttre, att_it cluster(state) safety(off)
di "VER B plain att_it   = " %20.16f el(e(wdidmainresults),1,1)

di "VER ===== V4 att + att_it together WITH an event-study window ====="
cap noi wooldid y unit t ttre, att att_it cluster(state) esprelength(2) espostlength(1) safety(off)
di "VER V4 rc=" _rc
cap noi qui wooldid y unit t ttre, att att_it cluster(state) safety(off)
di "VER V4 no-window rc=" _rc

di "VER ===== V2 shift equivariance, design A ====="
dgpA 900001 2
qui gen byte dd = (t>=ttre & ttre<.)
qui wooldid y unit t ttre, att cluster(cl) esprelength(2) espostlength(1) jointtests ///
    bootstrap(standard) reps(199) rseed(4242) dots(0) safety(off)
matrix A_cp = e(boot_ci_percentile)
matrix A_cb = e(boot_ci_bc)
matrix A_wm = e(boot_wald_mainresults)
matrix A_pp = e(boot_p_percentile)
matrix A_pb = e(boot_p_bc)
scalar A_jt = e(joint_horiz)
qui gen double y0 = y - 2*dd
qui wooldid y0 unit t ttre, att cluster(cl) esprelength(2) espostlength(1) jointtests ///
    bootstrap(standard) reps(199) rseed(4242) dots(0) safety(off)
matrix B_cp = e(boot_ci_percentile)
matrix B_cb = e(boot_ci_bc)
matrix B_wm = e(boot_wald_mainresults)
matrix B_pp = e(boot_p_percentile)
matrix B_pb = e(boot_p_bc)
scalar B_jt = e(joint_horiz)
di "VER shift b       " %14.4e (el(A_wm,1,1) - el(B_wm,1,1) - 2)
di "VER shift se      " %14.4e (el(A_wm,1,2) - el(B_wm,1,2))
di "VER shift pct_lb  " %14.4e (el(A_cp,1,1) - el(B_cp,1,1) - 2)
di "VER shift pct_ub  " %14.4e (el(A_cp,2,1) - el(B_cp,2,1) - 2)
di "VER shift bc_lb   " %14.4e (el(A_cb,1,1) - el(B_cb,1,1) - 2)
di "VER shift bc_ub   " %14.4e (el(A_cb,2,1) - el(B_cb,2,1) - 2)
di "VER shift wald_lb " %14.4e (el(A_wm,1,5) - el(B_wm,1,5) - 2)
di "VER shift wald_ub " %14.4e (el(A_wm,1,6) - el(B_wm,1,6) - 2)
di "VER A p_pct=" el(A_pp,1,1) " B p_pct=" el(B_pp,1,1) " diff=" %14.4e (el(A_pp,1,1)-el(B_pp,1,1))
di "VER A p_bc =" el(A_pb,1,1) " B p_bc =" el(B_pb,1,1) " diff=" %14.4e (el(A_pb,1,1)-el(B_pb,1,1))
di "VER A jt=" A_jt " B jt=" B_jt " diff=" %14.4e (A_jt-B_jt)
di "VER dual: A covers 2 = " (el(A_cp,1,1)<=2 & el(A_cp,2,1)>=2) "  B covers 0 = " (el(B_cp,1,1)<=0 & el(B_cp,2,1)>=0)
di "VER DONE"
