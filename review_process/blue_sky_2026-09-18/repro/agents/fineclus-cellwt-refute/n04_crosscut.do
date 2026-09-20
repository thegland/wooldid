*! n04 (NEW PROBE, not seeded) — "nothing is screened in its place, deliberately" (ado:6451).
*! The deleted refusals also refused CROSS-CUTTING clusters, not just finer ones.  The
*! commonest cross-cutting cluster in applied work is TIME.  Under the keep-i rule a
*! twice-drawn period arrives as duplicate (i,t) rows, and a never-drawn period disappears
*! from every unit at once - which moves the set of relative years the replication can
*! estimate.  Old ado:6364 refused this r459.  What does it do now?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 3131
set sortseed 3131
set obs 40
gen long unit = _n
gen onset = cond(unit <= 12, 4, cond(unit <= 24, 6, .))
gen double uf = rnormal(0,1)
gen double tau = cond(missing(onset), 0, 1 + rnormal(0,.4))
expand 9
bysort unit (onset): gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double y = 2 + uf + .1*t + tau*D + rnormal(0,1)
keep y unit t onset

dis "===== ANALYTIC cluster(unit) reference"
cap noi wooldid y unit t onset, att att_i cluster(unit)
dis "REF_RC=" _rc
mat Vr = e(V)
dis "ref se_att=" sqrt(Vr[1,1]) "  se_atti=" sqrt(Vr[2,2])

dis "===== BOOT cluster(unit) (nesting regime)"
cap noi wooldid y unit t onset, att att_i cluster(unit) bootstrap(standard) reps(199) rseed(21) dots(0)
dis "BU_RC=" _rc
mat Vu = e(V)
dis "cluster(unit) boot se_att=" sqrt(Vu[1,1]) "  se_atti=" sqrt(Vu[2,2]) "  df_r=" e(df_r)

dis "===== BOOT cluster(t) -- CROSS-CUTTING i.  Refused r459 before 2026-09-16."
cap noi wooldid y unit t onset, att att_i cluster(t) bootstrap(standard) reps(199) rseed(21) dots(0)
dis "BT_RC=" _rc
mat Vt = e(V)
dis "cluster(t) boot se_att=" sqrt(Vt[1,1]) "  se_atti=" sqrt(Vt[2,2]) "  df_r=" e(df_r)
dis "completed=" e(boot_completed) "  failed reps?"
dis "===== DONE"
