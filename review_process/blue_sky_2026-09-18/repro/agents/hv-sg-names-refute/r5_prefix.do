*=======================================================================
* r5_prefix.do -- NEW LEAD found while refuting hv-sgname-01/02.
*
* wooldid's e(b) column names are PREFIX-NESTED across estimands:
*     main_att  <  main_att_i  <  main_att_it  <  main_att_itime
*     pretreatment_att < pretreatment_att_i < ...
*     main_ame  <  main_ame_i  < ...
* Stata resolves a coefficient reference by UNIQUE PREFIX.  So on a fit
* that posts only ONE of the att family, a reference to a DIFFERENT
* member of the family resolves silently to the one that is there.
*
* HYPOTHESIS: `di _b[main_att]' after `wooldid ..., att_i' returns the
* att_i number with rc 0 and no message, even though no att estimand was
* requested and none is posted.
*
* Leg 0 establishes whether this is generic Stata `_b[]' behaviour, so
* the finding can be stated correctly.
* SELF-CONTAINED.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*===================== LEG 0: is _b[] prefix-resolving at all? =========
clear
set obs 500
gen double ax = rnormal()
gen double axx = rnormal()
gen double yy = ax + 2*axx + rnormal()
qui reg yy ax axx
noi dis "L0_exact_ax = " _b[ax]
noi dis "L0_exact_axx = " _b[axx]
cap noi dis "L0_prefix_a = " _b[a]
noi dis "L0_prefix_a_RC=" _rc
drop ax axx
gen double zzzz = rnormal()
gen double qq = rnormal()
replace yy = zzzz + qq + rnormal()
qui reg yy zzzz qq
cap noi dis "L0_uniqueprefix_zz = " _b[zz]
noi dis "L0_uniqueprefix_zz_RC=" _rc
noi dis "L0_exact_zzzz = " _b[zzzz]

*===================== panel =========================================
clear
set obs 192
gen long unit = _n
gen int grp = .
replace grp = ceil(_n/40)            if _n <= 120
replace grp = 3 + ceil((_n-120)/6)   if _n > 120 & _n <= 132
replace grp = 5 + ceil((_n-132)/30)  if _n > 132
gen int ttre = .
replace ttre = 2004 if inrange(grp,1,3)
replace ttre = 2007 if inrange(grp,4,5)
gen double truetau = cond(inrange(grp,1,3), 0.4, cond(inrange(grp,4,5), 2.5, 0))
gen double ai = rnormal()
expand 10
bysort unit: gen int year = 2000 + _n - 1
gen byte post = (ttre < .) & (year >= ttre)
gen double dose = cond(ttre<., 0.5 + 0.5*mod(unit,3), 0)
gen double y = ai + 0.05*(year-2000) + truetau*post + 0.20*rnormal()
compress

*===================== LEG 1: att ALONE (the reference numbers) ========
qui wooldid y grp year ttre, att esprelength(2) espostlength(2) robust
noi dis "L1_main_att = " _b[main_att]
scalar TRUE_ATT = _b[main_att]

*===================== LEG 2: att_i ALONE ==============================
qui wooldid y grp year ttre, att_i esprelength(2) espostlength(2) robust
tempname B2
matrix `B2' = e(b)
local cn : colnames `B2'
noi dis "L2_COLNAMES = `cn'"
noi dis "L2_main_att_i = " _b[main_att_i]
cap noi dis "L2_main_att  = " _b[main_att]
noi dis "L2_main_att_RC=" _rc
cap noi dis "L2_pretreatment_att = " _b[pretreatment_att]
noi dis "L2_pretreatment_att_RC=" _rc
cap noi lincom main_att
noi dis "L2_LINCOM_main_att_RC=" _rc
cap noi test main_att = 0
noi dis "L2_TEST_main_att_RC=" _rc
noi dis "L2_TRUE_ATT_from_LEG1 = " TRUE_ATT
noi dis "L2_RATIO_wrong_over_right = " _b[main_att]/TRUE_ATT
cap noi esttab ., keep(main_att) nostar b(6)
noi dis "L2_ESTTAB_keep_main_att_RC=" _rc

*===================== LEG 3: att_i + att_it (ambiguous prefix) ========
qui wooldid y grp year ttre, att_i att_it esprelength(2) espostlength(2) robust
cap noi dis "L3_main_att = " _b[main_att]
noi dis "L3_main_att_RC=" _rc

*===================== LEG 4: att + att_i (exact match should win) =====
qui wooldid y grp year ttre, att att_i esprelength(2) espostlength(2) robust
cap noi dis "L4_main_att = " _b[main_att]
noi dis "L4_main_att_RC=" _rc
noi dis "L4_main_att_i = " _b[main_att_i]

*===================== LEG 5: contreat, ame vs ame_i ===================
qui wooldid y grp year ttre, att contreat(dose) contreatcontrols(basic) ///
    esprelength(2) espostlength(2) robust
noi dis "L5_ame_from_att_run = " _b[main_ame]
scalar TRUE_AME = _b[main_ame]
qui wooldid y grp year ttre, att_i contreat(dose) contreatcontrols(basic) ///
    esprelength(2) espostlength(2) robust
tempname B5
matrix `B5' = e(b)
local cn5 : colnames `B5'
noi dis "L5_COLNAMES = `cn5'"
noi dis "L5_main_ame_i = " _b[main_ame_i]
cap noi dis "L5_main_ame  = " _b[main_ame]
noi dis "L5_main_ame_RC=" _rc
noi dis "L5_TRUE_AME_from_att_run = " TRUE_AME

*===================== LEG 6: subgroup names, prefix-nested too? =======
gen byte sg = 1 + mod(grp,2)
qui wooldid y grp year ttre, att_i subgroups(sg) esprelength(2) ///
    espostlength(2) robust
cap noi dis "L6_main_SG1_att = " _b[main_SG1_att]
noi dis "L6_main_SG1_att_RC=" _rc
noi dis "L6_main_SG1_att_i = " _b[main_SG1_att_i]

dis "R5_DONE"
