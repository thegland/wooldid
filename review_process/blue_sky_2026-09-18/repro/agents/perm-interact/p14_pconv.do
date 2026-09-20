*==============================================================================
* p14: are the permutation p-value CONVENTIONS, and the p<->interval coherence,
* the numbers they are documented to be?  Recomputed BY HAND from e(permdist_b)
* (/ e(permdist_se) under permtstat) and e(perm_comparator_b), for EVERY column.
*   exact enumeration  -> plain shares
*   sampled draws      -> add-one (1+count)/(n+1), except permrawp which is c/n
*   default p = min(1, 2 min(P_lo,P_hi)) ; permabsp = share(|s_r| >= |s_obs|)
*   permside(left)=P_lo, permside(right)=P_hi
* COHERENCE (help): "the interval excludes zero exactly when the p-value falls
*   below the level" - checked for the default and permabsp conventions only
*   (permside keeps a two-sided interval under a one-sided p, by design).
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242
clear
set obs 6
gen long unit = _n
gen double ttre = .
replace ttre = 2003 if inlist(unit,1,2)
replace ttre = 2005 if unit==3
expand 7
bysort unit: gen int t = 1999 + _n
gen double a_i = rnormal()
gen double g_t = 0.3*(t-2002)
gen double y = a_i + g_t + 1.2*(ttre<. & t>=ttre) + rnormal()*0.5
gen byte sg2 = mod(unit,2)
save "p14_data.dta", replace

mata:
void pconv(string scalar tag)
{
    real matrix D, S, C, Cs, P, CI, PLO, PHI
    real colvector s, keep
    string matrix ds
    real scalar k, n, plo, phi, pabs, praw, pexp, nbad, worst, d, exact, ts
    real scalar ncoh, nchk, rej, excl, alpha, lvl
    string scalar conv, worstnm

    exact = st_numscalar("e(perm_exact)")
    conv  = st_global("e(perm_pconv)")
    ts    = (st_global("e(perm_stat)") == "t")
    lvl   = st_numscalar("e(perm_level)")
    alpha = 1 - lvl/100
    D  = st_matrix("e(permdist_b)")
    ds = st_matrixcolstripe("e(permdist_b)")
    C  = st_matrix("e(perm_comparator_b)")
    P  = st_matrix("e(perm_p)")
    CI = st_matrix("e(perm_ci)")
    PLO= st_matrix("e(perm_p_lo)")
    PHI= st_matrix("e(perm_p_hi)")
    if (ts) {
        S  = st_matrix("e(permdist_se)")
        Cs = st_matrix("e(perm_comparator_se)")
        D  = D :/ S
        C  = C :/ Cs
    }
    nbad = 0; worst = 0; worstnm = ""; nchk = 0; ncoh = 0
    for (k=1; k<=cols(D); k++) {
        if (P[1,k] >= .) continue
        if (C[1,k] >= .) continue
        keep = (D[.,k] :< .)
        if (ts) keep = keep :& (D[.,k] :!= .)
        s = select(D[.,k], keep)
        if (rows(s) < 2) continue
        n = rows(s)
        if (exact == 1) {
            plo  = sum(s :<= C[1,k])/n
            phi  = sum(s :>= C[1,k])/n
            pabs = sum(abs(s) :>= abs(C[1,k]))/n
        }
        else {
            plo  = (1 + sum(s :<= C[1,k]))/(n+1)
            phi  = (1 + sum(s :>= C[1,k]))/(n+1)
            pabs = (1 + sum(abs(s) :>= abs(C[1,k])))/(n+1)
        }
        praw = sum(abs(s) :>= abs(C[1,k]))/n
        if      (conv == "absolute") pexp = pabs
        else if (conv == "raw")      pexp = praw
        else if (conv == "left")     pexp = plo
        else if (conv == "right")    pexp = phi
        else                         pexp = min((1, 2*min((plo,phi))))
        nchk = nchk + 1
        d = abs(pexp - P[1,k])
        if (d > worst) {
            worst = d
            worstnm = ds[k,2]
        }
        if (d > 1e-12) nbad = nbad + 1
        /* p <-> interval coherence, default and absolute conventions only */
        if ((conv == "equaltailed" | conv == "absolute") & CI[1,k] < . & CI[2,k] < .) {
            excl = (CI[1,k] > 0 | CI[2,k] < 0)
            rej  = (P[1,k] < alpha)
            if (excl != rej) {
                ncoh = ncoh + 1
                printf("PCOH|%s|%s|p=%9.6f|ci=[%10.6f,%10.6f]|excl=%f rej=%f\n",
                       tag, ds[k,2], P[1,k], CI[1,k], CI[2,k], excl, rej)
            }
        }
    }
    printf("PCONV|%s|conv=%s|exact=%f|stat=%s|ncols=%f|checked=%f|MISMATCH=%f|WORSTABSDIFF=%12.3e|at=%s|COHERENCE_VIOL=%f\n",
           tag, conv, exact, st_global("e(perm_stat)"), cols(D), nchk, nbad, worst, worstnm, ncoh)
    /* and the one-sided shares by name */
    k = 0
    for (n=1; n<=cols(D); n++) if (ds[n,2]=="main_att") k = n
    if (k > 0) {
        keep = (D[.,k] :< .)
        s = select(D[.,k], keep)
        n = rows(s)
        if (exact == 1) {
            plo = sum(s :<= C[1,k])/n
            phi = sum(s :>= C[1,k])/n
        }
        else {
            plo = (1+sum(s :<= C[1,k]))/(n+1)
            phi = (1+sum(s :>= C[1,k]))/(n+1)
        }
        printf("PLOHI|%s|hand_lo=%9.6f rep_lo=%9.6f | hand_hi=%9.6f rep_hi=%9.6f\n",
               tag, plo, PLO[1,k], phi, PHI[1,k])
    }
}
end

cap prog drop pc
prog def pc
  syntax anything(name=tag) [, OPTS(string asis) SAMP(integer 0)]
  if `samp' == 0 local eng permexact
  else local eng permreps(`samp') rseed(99)
  cap noi wooldid y unit t ttre , att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) `eng' dots(0) `opts'
  dis "PCONV|`tag'|RC=" _rc
  if _rc != 0 exit
  mata: pconv("`tag'")
end

dis "########## EXACT ENUMERATION (plain shares) ##########"
pc x_def
pc x_abs   , opts(permabsp)
pc x_raw   , opts(permrawp)
pc x_left  , opts(permside(left))
pc x_right , opts(permside(right))
pc x_ts    , opts(permtstat)
pc x_ts_abs, opts(permtstat permabsp)
pc x_lvl90 , opts(level(90))

dis "########## SAMPLED DRAWS (add-one, except permrawp) ##########"
pc s_def   , samp(80)
pc s_abs   , samp(80) opts(permabsp)
pc s_raw   , samp(80) opts(permrawp)
pc s_left  , samp(80) opts(permside(left))
pc s_right , samp(80) opts(permside(right))
pc s_ts    , samp(80) opts(permtstat)
pc s_sg    , samp(80) opts(subgroups(sg2))
dis "P14_DONE"
