*! LENS-0 (independent reproduction) of PERMSEM-01 -- run B (fixes a non-conformable
*! select() in run A's mata block; identical data and identical seeds, so run A's
*! sentinels reproduce).  Same rolling-administrative-frame design as L0_permrows.do:
*! 30 districts 2010-2019, 12 full-window, 10 entering the frame in 2015, 8 in 2017;
*! three adopt in 2017, all three on full-window districts in the OBSERVED assignment.
*! Run A measured: reps=120 completed=120 FAILED=0, but 23 replications PARTIAL
*! (e(perm_n) min 97, max 120) -- a shape neither earlier repro produced.
*! This run finishes the arithmetic run A could not reach.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771103
set sortseed 771103
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-lens0"

which wooldid

clear
set obs 30
gen long dist = _n
gen int firstyr = 2010
replace firstyr = 2015 if inrange(dist, 13, 22)
replace firstyr = 2017 if inrange(dist, 23, 30)
gen double onset = cond(inlist(dist, 3, 7, 11), 2017, .)
gen double a_i = rnormal()*1.5
expand 10
bysort dist: gen int year = 2010 + _n - 1
drop if year < firstyr
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.20*(year - 2010) + 0.8*post + rnormal()
qui count
dis "LB_NOBS=" r(N)

cap noi wooldid y dist year onset, cluster(dist) att esprelength(3) espostlength(2) ///
    jointtests adversarial esfixedbaseperiod permute(shuffle) permreps(120) ///
    rseed(4041) permsave(L0b_save) replace
dis "LB_RC=" _rc
dis "LB_reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed)
dis "LB_incomplete=|" e(perm_incomplete) "|  addone=|" e(perm_addone) "|"

matrix _PDB = e(permdist_b)
matrix _PD  = e(permdist)
matrix _PN  = e(perm_n)
matrix _JW  = e(perm_jointW)
matrix _JWO = e(perm_jointWobs)
matrix _AW  = e(perm_advW)
matrix _AWO = e(perm_advWobs)
scalar _NREPS = e(perm_reps)
scalar _NCOMP = e(perm_completed)
scalar _NFAIL = e(perm_failed)
dis "LB_DIMS permdist_b=" rowsof(_PDB) "x" colsof(_PDB) " permdist=" rowsof(_PD) "x" colsof(_PD) ///
    " jointW=" rowsof(_JW) "x" colsof(_JW) " advW=" rowsof(_AW) "x" colsof(_AW)
matrix list _PN, format(%9.0f)

* the file
preserve
use L0b_save.dta, clear
qui ds _rep onset_*, not
local coefvars `r(varlist)'
local ncoef : word count `coefvars'
dis "LB_COEFVARS=`coefvars'"
qui gen int nmiss = 0
foreach v of local coefvars {
  qui replace nmiss = nmiss + missing(`v')
}
qui count if _rep >= 1 & nmiss == `ncoef'
dis "LB_FILE_ALLMISS=" r(N)
qui count if _rep >= 1 & nmiss > 0 & nmiss < `ncoef'
dis "LB_FILE_PARTIAL=" r(N)
* which columns vanish on the partial rows
foreach v of local coefvars {
  qui count if _rep >= 1 & missing(`v')
  dis "LB_COLMISS `v' = " r(N)
}
sort _rep
qui drop if _rep == 0
mkmat `coefvars', matrix(_FALL)
restore

mata:
D  = st_matrix("_PDB")
FA = st_matrix("_FALL")
JW = st_matrix("_JW"); JWO = st_matrix("_JWO")
AW = st_matrix("_AW"); AWO = st_matrix("_AWO")
NR = st_numscalar("_NREPS"); NC = st_numscalar("_NCOMP"); NF = st_numscalar("_NFAIL")

// C1: permdist_b vs the file's kept rows
ok = (rowsum(FA :< .) :> 0)
FC = select(FA, ok)
printf("LB_C1 file-kept rows=%g permdist_b rows=%g reps=%g completed=%g failed=%g\n",
       rows(FC), rows(D), NR, NC, NF)
if (rows(FC) == rows(D) & cols(FC) == cols(D)) {
  dd = abs(FC - D)
  v  = vec(dd)
  v  = select(v, v :< .)
  printf("LB_C1 maxabsdiff(permdist_b, file block)=%g over %g finite cells\n",
         (rows(v) > 0 ? max(v) : 0), rows(v))
  printf("LB_C1 missing cells: file=%g permdist_b=%g\n", sum(FC :>= .), sum(D :>= .))
}

// C2/C3: the W matrices against the file's _rep rows
jn = st_matrixcolstripe("_JW")
an = st_matrixcolstripe("_AW")
allmiss = (rowsum(FA :< .) :== 0)
for (k = 1; k <= cols(JW); k++) {
  m = (JW[.,k] :>= .)
  printf("LB_C2 jointW col %s: rows=%g missing=%g | reps=%g completed=%g failed=%g | file allmissing rows=%g\n",
         jn[k,2], rows(JW), sum(m), NR, NC, NF, sum(allmiss));
  printf("LB_C3 jointW %s: rows missing in W but KEPT in the file = %g ; present in W but all-missing in file = %g\n",
         jn[k,2], sum(m :& !allmiss), sum(!m :& allmiss))
  nn   = sum(JW[.,k] :< .)
  exc  = sum(select(JW[.,k], JW[.,k] :< .) :>= JWO[1,k])
  exn  = sum(JW[.,k] :>= JWO[1,k])
  pc   = (1 + exc)/(nn + 1)
  pn   = (1 + exn)/(rows(JW) + 1)
  printf("LB_C4 %s n_correct=%g nex_correct=%g | N_naive=%g nex_naive=%g | overcount=%g e(perm_failed)=%g reps-n=%g\n",
         jn[k,2], nn, exc, rows(JW), exn, exn - exc, NF, NR - nn)
  printf("LB_C4 %s p_correct=%12.8f (%s)  p_naive=%12.8f (%s)\n",
         jn[k,2], pc, strofreal(pc,"%21x"), pn, strofreal(pn,"%21x"))
  printf("LB_C5 %s POSTED=%12.8f (%s)  matches_correct=%g matches_naive=%g\n",
         jn[k,2], st_numscalar("e(" + jn[k,2] + ")"), strofreal(st_numscalar("e(" + jn[k,2] + ")"),"%21x"),
         (st_numscalar("e(" + jn[k,2] + ")") == pc), (st_numscalar("e(" + jn[k,2] + ")") == pn))
}
for (k = 1; k <= cols(AW); k++) {
  m = (AW[.,k] :>= .)
  nn   = sum(AW[.,k] :< .)
  exc  = sum(select(AW[.,k], AW[.,k] :< .) :>= AWO[1,k])
  exn  = sum(AW[.,k] :>= AWO[1,k])
  pc   = (1 + exc)/(nn + 1)
  pn   = (1 + exn)/(rows(AW) + 1)
  printf("LB_ADV %s: rows=%g missing=%g kept-in-file-but-missing-in-W=%g\n",
         an[k,2], rows(AW), sum(m), sum(m :& !allmiss))
  printf("LB_ADV %s n=%g p_correct=%12.8f p_naive=%12.8f overcount=%g POSTED=%12.8f matches_correct=%g\n",
         an[k,2], nn, pc, pn, exn - exc, st_numscalar("e(" + an[k,2] + ")"),
         (st_numscalar("e(" + an[k,2] + ")") == pc))
}
end

dis "LB_DONE"
