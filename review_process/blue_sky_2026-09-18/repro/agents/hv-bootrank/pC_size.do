*=======================================================================
* pC_size.do -- the EMPIRICAL size leg, using wooldid's OWN bootstrap
* draws so the answer does not rest on a normality assumption.
*
* One wooldid run at reps(500) gives e(bootdist): 500 replication draws
* of every coefficient.  Recentre them at their own mean (so the null
* "coefficient = 0" is TRUE by construction in the recentred draws),
* then repeatedly:
*     take 1 draw as the "estimate"  b*
*     take B OTHER draws to build the covariance V*
*     W = b*' inv(V*) b* ,  F = W / q
*     reject when Ftail(q, e(df_r), F) < .05
* That is exactly the arithmetic a user's `test' runs after a
* bootstrap(standard) fit with e(boot_N_V) = B.  The rejection rate is
* therefore the ACTUAL size of that hand test at nominal 5%.
*
* SELF-CONTAINED: needs gen_panel.do next to it.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

local pN 40
local pT 12
local pG 20
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank/gen_panel.do"

dis "=== BIGRUN_BEGIN ==="
cap noisily wooldid y unit t ttre, att esfixedbaseperiod ///
    esprelength(5) espostlength(3) cluster(clu) ///
    bootstrap(standard) reps(500) rseed(77)
dis "BIGRUN_RC=" _rc
dis "BIGRUN_COMPLETED=" e(boot_completed)
dis "BIGRUN_N_V=" e(boot_N_V)
dis "BIGRUN_DFR=" e(df_r)
dis "=== BIGRUN_END ==="

mata:
BD  = st_matrix("e(bootdist)")
cs  = st_matrixcolstripe("e(bootdist)")
dfr = st_numscalar("e(df_r)")
printf("BD dims %f x %f, e(df_r)=%f\n", rows(BD), cols(BD), dfr)
// the pre-period block: a TRUE null in this DGP and the block a user
// actually joint-tests
sel = selectindex(strpos(cs[.,2]', "ES_att_pre") :== 1)
printf("q = %f pre columns: \n", length(sel))
cs[sel, 2]'
A = BD[., sel]
A = select(A, rowmissing(A) :== 0)
q = cols(A)
N = rows(A)
printf("usable complete draws N=%f  q=%f\n", N, q)
A = A :- mean(A)                       // recentre: the null is now TRUE

M = 4000
rseed(20260917)
printf("\n%6s %8s %10s %14s\n", "B", "q", "df_r", "actual size")
Bl = (8, 12, 20, 30, 60, 120, 250)
for (bi = 1; bi <= cols(Bl); bi++) {
    B = Bl[bi]
    if (B + 1 > N) continue
    rej = 0
    nok = 0
    for (m = 1; m <= M; m++) {
        p   = order(runiform(N,1), 1)     // random permutation of 1..N
        bs  = A[p[1], .]'
        Vs  = variance(A[p[2::(B+1)], .])
        Vi  = invsym(Vs)
        if (diag0cnt(Vi) > 0) continue     // singular: a user's test would drop constraints
        W   = (bs' * Vi * bs)
        F   = W / q
        nok = nok + 1
        if (Ftail(q, dfr, F) < 0.05) rej = rej + 1
    }
    printf("SIZE B=%6.0f q=%3.0f dfr=%5.0f nonsingular=%6.0f actual_size=%8.4f\n", B, q, dfr, nok, rej/nok)
}
end

dis "PROBE_C_DONE"
