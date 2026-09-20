*=======================================================================
* pF_selection.do -- is the gap between the PRINTED standard error and
* sqrt(diag(e(V))) sampling noise, or SELECTION?
*
* e(V) is the covariance over the listwise complete-case replications
* only (e(boot_N_V) of them).  Those rows are not a random subset: they
* are the replications in which EVERY column was estimable, i.e. the
* resamples that happened to draw the thin clusters.  If that selection
* is informative, the listwise rows are systematically LESS dispersed
* than the full draw set and e(V) is biased DOWNWARD.
*
* Test, entirely inside one run, from e(bootdist):
*   sd_full[k]     = sd over every available row of column k
*   sd_listwise[k] = sd over the listwise-complete rows only
*   null           = sd over 2000 RANDOM subsets of the same size
* If sd_listwise sits in the extreme lower tail of that null, the
* selection is informative and the shrinkage is not chance.
*
* SELF-CONTAINED: needs gen_admin.do next to it.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank/gen_admin.do"

cap noisily wooldid y state t ttre, att esfixedbaseperiod ///
    esprelength(6) espostlength(11) cluster(state) subgroups(reg) ///
    bootstrap(standard) reps(60) rseed(4242)
dis "F_RC=" _rc
dis "F_COMPLETED=" e(boot_completed)
dis "F_N_V=" e(boot_N_V)

mata:
BD   = st_matrix("e(bootdist)")
nm   = st_matrixcolstripe("e(bootdist)")[., 2]
comp = rows(BD)
lw   = rowmissing(BD) :== 0
nlw  = sum(lw)
printf("F_BDROWS=%f  F_NLW=%f\n", comp, nlw)

// pooled ES_att columns plus the two level rows: what a user reads
want = selectindex((strpos(nm, "ES_att_") :== 1) :| (nm :== "main_att") :| (nm :== "pretreatment_att"))
printf("F_NCOLS_EXAMINED=%f\n", length(want))

rseed(20260917)
M = 2000
nlo = 0
ntot = 0
printf("\n%-24s %10s %10s %8s %10s %8s\n", "coefficient", "sd_full", "sd_lw", "ratio", "p_lower", "navail")
for (j = 1; j <= length(want); j++) {
    k = want[j]
    col = BD[., k]
    av  = selectindex(col :< .)
    if (length(av) < 3) continue
    sdf = sqrt(variance(col[av]))
    sel = selectindex(lw)
    sdl = sqrt(variance(col[sel]))
    // null: random subsets of size nlw drawn from the column's AVAILABLE rows
    cnt = 0
    for (m = 1; m <= M; m++) {
        p = order(runiform(length(av), 1), 1)
        sub = col[av[p[1::min((nlw, length(av)))]]]
        if (sqrt(variance(sub)) <= sdl) cnt = cnt + 1
    }
    plow = cnt / M
    ntot = ntot + 1
    if (plow <= 0.05) nlo = nlo + 1
    printf("FSEL %-20s %10.6f %10.6f %8.3f %10.4f %8.0f\n", nm[k], sdf, sdl, sdl/sdf, plow, length(av))
}
printf("\nF_SUMMARY columns=%f  with p_lower<=.05 = %f\n", ntot, nlo)
end

dis "PROBE_F_DONE"
