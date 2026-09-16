* mc_ri_ext.do - randomization-inference EXTENSIONS with no MC coverage at HEAD (Stage C
* scaffold, written 2026-09-13 from the Stage A audit: mc_drivers F5 as amended by the critic).
* Conventions follow mc_ri2.do: positional args, one CSV row per (rep, arm) appended with
* open/write/close, every surface read BY NAME with a positional fallback recorded in namesok.
*
* Phases (each its own DGP + arm; TAU is the planted post-treatment effect where stated):
*   E1 contreat  permute(shuffle) with contreat(dose) contreatcontrols(basic) on a dta_cont-LIKE
*                DGP: a continuous dose that varies across units WITHIN each cohort, entering
*                the outcome as tau*dose on treated post cells (dta_cont: treated*intensity).
*                Cohort layout = mc_ri.do's R1 (12 cohorts, onset labels {5 x3, 7 x3, . x6}
*                shuffled uniformly; assignment space 18480 > the M1 prereg's ~3000 tie floor)
*                rather than dta_cont's 10 cohorts / 5 treated (space C(10,5) = 252, far below
*                it), with 5 units per cohort as dta_cont has multiple observations per (i,t).
*                tau = 0: SIZE of p for main_att, main_ame, ES_ame_post0/1, pretreatment_att
*                and the C-family joint p e(joint_Choriz), all exact on the lattice.
*   E2 right     R1 DGP, tau = 0, permside(right): p = (1 + #{s_r >= s_obs})/(B + 1); exact size
*                (P(p <= a) = floor(a(B+1))/(B+1)); e(perm_pconv) must read "right".
*   E3 level90   R1 DGP, tau = TAU (default 1), level(90): COVERAGE of tau by the shift interval
*                e(perm_ci) at 90 (the only level() != 95 arm in any MC).
*   E4 strata    mc_ri.do's R4 DGP (16 clusters, 2 strata of 8, 3 treated per stratum at onset
*                5 / 7), tau = 0, permstrata(str) at B draws (R4 was B = 99): exact size.
*   E5 tstat0    R1 DGP, tau = 0, permtstat: SIZE of the studentized p (the tau = 0 companion
*                of mc_ri2 D4, which measured only coverage at tau = 1).
*   E6 fineclus  ADDED 2026-09-16 for ticket T1 (PREREG_FEATURE_FINECLUS_2026-09-16.md F2/P1).
*                The R1 DGP UNCHANGED and tau = 0, but the call is the coarse-i / fine-cluster
*                pairing the help recommends and the engine refused until T1: i is the onset
*                COHORT (3 levels - onset 5, onset 7, never-treated) and cluster() is the unit
*                (12), strictly finer, under permonsetcohorts. Every replication therefore
*                rebuilds i from the onsets it drew. SIZE of p, judged by the SAME lattice rule
*                and the SAME bands as E2/E4/E5, which are its controls in the same campaign:
*                a wrong rebuild that is merely different still produces plausible p-values, so
*                exactness under the sharp null is the readout that detects it and the only one.
*                The assignment space is R1's 18480, unchanged - the rebuild regroups units,
*                it does not touch the generator.
* Recorded per row: ncomp/nfail/npos, namesok, e(perm_pconv), the arm's level, e(perm_p) by
* name for main_att / pretreatment_att / ES_att_post0 / ES_att_post1 / main_ame /
* ES_ame_post0 / ES_ame_post1 (ame columns '.' outside E1), e(perm_p_lo/hi) main_att,
* e(joint_horiz), e(joint_Choriz), e(perm_ci) ll/ul for main_att and main_ame, e(perm_sd)
* main_att, e(b) main_att / main_ame, and the usable draw count for main_att from e(permdist_b).
*
* args: 1=ado 2=shard 3=nreps 4=baseseed 5=outdir 6=B (permreps; production 199)
*       7=TAU (E3's planted effect; default 1) 8=phases (+-separated list, e.g. E1+E2 - a comma
*       would be read by `do` as the start of options; default E1+E2+E3+E4+E5)
* Production cost at B=199: E2-E5 ~1 min per rep (mc_ri2's D1 is 3.5 s at 199 x fast mode on the
* same DGP - E2/E3/E5 are that; E4 slightly more); E1 (contreat + margins per replication)
* several times that - measure in the smoke and shard accordingly.
*
* ACCEPTANCE (to be preregistered before the production run; analyze_mc.py's binomial regions):
*   - E1, E2, E4, E5: rejection COUNT at .05 and .10 inside the central ~95% region of
*     Binomial(M, alpha) for every recorded p (E1: att, ame, joint_Choriz; the lattice makes the
*     nominal size exact at B = 199, so any departure is a real finding);
*   - E3: coverage COUNT of tau by [ll_main, ul_main] inside the Binomial(M, .90) region, and
*     the blank-interval share reported (the shift interval is blank when the null band excludes
*     0, mc_ri2's off-centre case);
*   - integrity: rc = 0, nfail = 0 and ncomp = B on every row, namesok = 1, pconv = "right" on
*     E2 rows, ndraw = B for main_att.
*
* USAGE (2026-09-13). Positional, shared with mc_bands.do and mc_adv_bs12.do; the runner is
* run_stagec.sh, which honours the machine's Stata instance cap:
*     ./run_stagec.sh riext <outdir> <nreps> <B> <tau> <phases|-> <shard>...
* Direct invocation:
*     <stata> -e do mc_ri_ext.do <ado> <shard> <nreps> <seedbase> <outdir> <B> <tau> <phases>
* phases is E1+E2+E3+E4+E5, + or , separated. B is the permutation replication count and
* wooldid's permreps() floor is 20, so B below 20 makes the permutation arms refuse rather than
* run. Writes <outdir>/mc_ri_ext_s<shard>.csv, one row per (rep, arm). SMOKE-RUN 2026-09-13 at
* nreps 1, B 20: rc 0, 5 rows. NOT yet campaigned and there is no analyzer: a scaffold, not
* evidence.
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local SHARD = `2'
local NREPS = `3'
local BASE  = `4'
local OUTDIR `5'
local B = `6'
local TAU `7'
if "`TAU'" == "" local TAU 1
local PHASES `8'
* the default list is UNCHANGED by the T1 addition: E6 is opted into by name, so a rerun of
* the existing campaign reproduces the existing arms and nothing else
if "`PHASES'" == "" local PHASES E1+E2+E3+E4+E5
local PHASES = subinstr(subinstr("`PHASES'", "+", " ", .), ",", " ", .)

global MCRIX_OUT `"`OUTDIR'/mc_ri_ext_s`SHARD'.csv"'
global MCRIX_SHARD = `SHARD'

*----------------------------------------------------------------- DGPs
* R1 (mc_ri2.do's mkdgp): 12 clusters, labels {5 x3, 7 x3, . x6}, 10 periods, tau on post
cap program drop mkR1
program define mkR1
  args seed tau
  clear
  set seed `seed'
  qui set obs 12
  gen i = _n
  gen u = runiform()
  sort u
  gen ttre = cond(_n <= 3, 5, cond(_n <= 6, 7, .))
  sort i
  gen a = rnormal()
  qui expand 10
  bysort i: gen t = _n
  gen y = a + .1*t + rnormal()
  if `tau' != 0 qui replace y = y + `tau'*(t >= ttre & ttre < .)
  drop a u
  sort i t
  ereturn clear
end

* R4 (mc_ri.do): 16 clusters, 2 strata of 8, 3 treated per stratum (onset 5 in stratum 1, 7 in 2)
cap program drop mkR4
program define mkR4
  args seed
  clear
  set seed `seed'
  qui set obs 16
  gen i = _n
  gen str = 1 + (_n > 8)
  gen u = runiform()
  sort str u
  by str: gen rk = _n
  gen ttre = cond(str == 1 & rk <= 3, 5, cond(str == 2 & rk <= 3, 7, .))
  sort i
  gen a = rnormal()
  qui expand 10
  bysort i: gen t = _n
  gen y = a + .1*t + rnormal()
  drop a u rk
  sort i t
  ereturn clear
end

* CONT: R1's cohort layout, 5 units per cohort, unit-level dose varying within cohort,
* effect tau*dose on treated post cells (dta_cont-like: y = ... + treated*intensity)
cap program drop mkCONT
program define mkCONT
  args seed tau
  clear
  set seed `seed'
  qui set obs 12
  gen i = _n
  gen u = runiform()
  sort u
  gen ttre = cond(_n <= 3, 5, cond(_n <= 6, 7, .))
  sort i
  gen a = rnormal()
  qui expand 5
  bysort i: gen unit = _n
  gen dose = 1 + runiform()
  qui expand 10
  bysort i unit: gen t = _n
  gen y = a + .1*t + rnormal()
  if `tau' != 0 qui replace y = y + `tau'*dose*(t >= ttre & ttre < .)
  drop a u
  sort i unit t
  ereturn clear
end

*----------------------------------------------------------------- row writer
cap program drop emitrow
program define emitrow
  args phase arm tau rep dgpseed rseed rc secs level
  local FLD ncomp nfail npos namesok ///
    p_main p_pre p_post0 p_post1 p_ame p_amepost0 p_amepost1 plo_main phi_main ///
    jointp jointCp ll_main ul_main ll_ame ul_ame sd_main b_main b_ame ndraw ///
    permblank n_main n_post1 n_amepost0
  foreach nm of local FLD {
    local `nm' .
  }
  local pconv
  if `rc' == 0 {
    local pconv `e(perm_pconv)'
    cap local ncomp = e(perm_completed)
    cap local nfail = e(perm_failed)
    cap local npos  = e(perm_npossible)
    cap local jointp = e(joint_horiz)
    cap local jointCp = e(joint_Choriz)
    tempname P
    cap matrix `P' = e(perm_p)
    if _rc == 0 {
      local namesok 1
      local i_main  = colnumb(`P', "main_att")
      local i_pre   = colnumb(`P', "pretreatment_att")
      local i_post0 = colnumb(`P', "ES_att_post0")
      local i_post1 = colnumb(`P', "ES_att_post1")
      local i_ame   = colnumb(`P', "main_ame")
      local i_amepost0 = colnumb(`P', "ES_ame_post0")
      local i_amepost1 = colnumb(`P', "ES_ame_post1")
      * the att family is always present: positional fallback 1..4 flags namesok 0
      local kk 0
      foreach z in main pre post0 post1 {
        local ++kk
        if `i_`z'' >= . | `i_`z'' == 0 {
          local i_`z' `kk'
          local namesok 0
        }
      }
      foreach z in main pre post0 post1 ame amepost0 amepost1 {
        if `i_`z'' < . & `i_`z'' > 0 local p_`z' = `P'[1, `i_`z'']
      }
      * CONTRACT CHANGE 2026-09-15. e(perm_n) is a MATRIX of PER-COLUMN draw counts and
      * e(perm_blank) a scalar; both are posted by the package (wooldid.ado:7971, :8031) and
      * neither was ever recorded here. Without them the lattice check - which reads main_att
      * alone - proves nothing about the OTHER columns' denominators: an on-lattice p is
      * reachable at several n below 199, so a column that quietly lost draws is
      * indistinguishable from one that did not. That was the one gap the 2026-09-14 source
      * read could not close on the E1 over-rejection, so the counts are recorded from now on.
      tempname PN
      cap matrix `PN' = e(perm_n)
      if _rc == 0 {
        foreach z in main post1 amepost0 {
          if `i_`z'' < . & `i_`z'' > 0 local n_`z' = `PN'[1, `i_`z'']
        }
      }
      cap local permblank = e(perm_blank)
      tempname PL PH CI SD BB PD
      cap matrix `PL' = e(perm_p_lo)
      if _rc == 0 local plo_main = `PL'[1, `i_main']
      cap matrix `PH' = e(perm_p_hi)
      if _rc == 0 local phi_main = `PH'[1, `i_main']
      cap matrix `CI' = e(perm_ci)
      if _rc == 0 {
        local rll = rownumb(`CI', "ll")
        local rul = rownumb(`CI', "ul")
        if `rll' >= . local rll 1
        if `rul' >= . local rul 2
        local ll_main = `CI'[`rll', `i_main']
        local ul_main = `CI'[`rul', `i_main']
        if `i_ame' < . & `i_ame' > 0 {
          local ll_ame = `CI'[`rll', `i_ame']
          local ul_ame = `CI'[`rul', `i_ame']
        }
      }
      cap matrix `SD' = e(perm_sd)
      if _rc == 0 local sd_main = `SD'[1, `i_main']
      cap matrix `BB' = e(b)
      if _rc == 0 {
        local kb = colnumb(`BB', "main_att")
        if `kb' < . local b_main = `BB'[1, `kb']
        local kb = colnumb(`BB', "main_ame")
        if `kb' < . local b_ame = `BB'[1, `kb']
      }
      cap matrix `PD' = e(permdist_b)
      if _rc == 0 {
        local kd = colnumb(`PD', "main_att")
        if `kd' >= . local kd `i_main'
        mata: st_local("ndraw", strofreal(sum(st_matrix("`PD'")[., `kd'] :< .)))
      }
    }
  }
  local row "`phase',`arm',`tau',${MCRIX_SHARD},`rep',`dgpseed',`rseed',`rc',`secs',`pconv',`level'"
  foreach nm of local FLD {
    local x = strtrim(string(``nm'', "%18.12g"))
    if "`x'" == "" local x .
    local row "`row',`x'"
  }
  tempname fh
  file open `fh' using "${MCRIX_OUT}", write append
  file write `fh' `"`row'"' _n
  file close `fh'
end

*----------------------------------------------------------------- header (once per shard)
cap confirm file "${MCRIX_OUT}"
if _rc {
  tempname fh0
  file open `fh0' using "${MCRIX_OUT}", write replace
  file write `fh0' "phase,arm,tau,shard,rep,dgpseed,rseed,rc,secs,pconv,level," ///
    "ncomp,nfail,npos,namesok," ///
    "p_main,p_pre,p_post0,p_post1,p_ame,p_amepost0,p_amepost1,plo_main,phi_main," ///
    "jointp,jointCp,ll_main,ul_main,ll_ame,ul_ame,sd_main,b_main,b_ame,ndraw," ///
    "permblank,n_main,n_post1,n_amepost0" _n
  file close `fh0'
}

local ESOPT esprelength(2) espostlength(2)

*----------------------------------------------------------------- the loop
foreach ph of local PHASES {
  forval rep = 1/`NREPS' {
    local dgpseed = `BASE' + `SHARD'*100000 + `rep'
    local rseed   = `BASE' + 500000 + `SHARD'*100000 + `rep'
    local level 95
    local tau 0
    timer clear 91
    timer on 91
    if "`ph'" == "E1" {
      local arm contreat
      mkCONT `dgpseed' 0
      cap qui wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) ///
          `ESOPT' jointtests permute(shuffle) permreps(`B') rseed(`rseed') dots(0)
    }
    else if "`ph'" == "E2" {
      local arm right
      mkR1 `dgpseed' 0
      cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
          permute(shuffle) permreps(`B') permside(right) rseed(`rseed') dots(0)
    }
    else if "`ph'" == "E3" {
      local arm level90
      local level 90
      local tau `TAU'
      mkR1 `dgpseed' `TAU'
      cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
          permute(shuffle) permreps(`B') level(90) rseed(`rseed') dots(0)
    }
    else if "`ph'" == "E4" {
      local arm strata
      mkR4 `dgpseed'
      cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
          permute(shuffle) permreps(`B') permstrata(str) rseed(`rseed') dots(0)
    }
    else if "`ph'" == "E5" {
      local arm tstat0
      mkR1 `dgpseed' 0
      cap qui wooldid y i t ttre, att cluster(i) `ESOPT' jointtests ///
          permute(shuffle) permreps(`B') permtstat rseed(`rseed') dots(0)
    }
    else if "`ph'" == "E6" {
      * T1: the same R1 data, the same seeds, the same statistic - only the (i, cluster)
      * pairing moves. icoh is R1's onset grouping (a deterministic function of ttre, so the
      * DGP is untouched); cluster(i) is then strictly finer than i and the run is reachable
      * only under permonsetcohorts
      local arm fineclus
      mkR1 `dgpseed' 0
      qui gen icoh = cond(ttre == 5, 1, cond(ttre == 7, 2, 3))
      cap qui wooldid y icoh t ttre, att cluster(i) `ESOPT' jointtests ///
          permute(shuffle) permreps(`B') permonsetcohorts rseed(`rseed') dots(0)
    }
    else {
      di as error "mc_ri_ext: unknown phase `ph'"
      exit 198
    }
    local rc = _rc
    timer off 91
    qui timer list 91
    if `rep' == 1 & `rc' == 0 {
      local cn : colnames e(perm_p)
      di "MCRIX_PERMP_COLS `ph': `cn'"
    }
    emitrow `ph' `arm' `tau' `rep' `dgpseed' `rseed' `rc' `=r(t91)' `level'
    di "MCRIX_PROGRESS `ph' shard=`SHARD' rep=`rep' rc=`rc' secs=" %6.1f r(t91)
  }
  di "MCRIX_PHASE_DONE `ph' shard=`SHARD'"
}
di "MCRIX_SHARD_DONE `SHARD'"
exit, clear
