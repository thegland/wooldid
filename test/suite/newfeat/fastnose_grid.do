* fastnose_grid.do — M4 fastnose parity grid (preregistered: test/suite/PREREG_NEWFEAT_2026-08-24.md §M4)
*
* For each of 23 grid configs, wooldid runs FULL and with `fastnose` added, and every
* shared e(b) entry is compared MATCHED BY COLUMN NAME (never position) — plus
* e(gteffects) cell-by-cell when the config includes gteffects. PASS iff for every
* entry |diff| <= 1e-12 OR reldif <= 1e-10. Any missing-vs-nonmissing flip or
* column-name-set mismatch = FAIL naming the config, column, and both values.
* Worst |diff| / reldif per config is displayed even on pass (FGRID info lines).
*
* reldif here is the MAGNITUDE-relative difference |a-b|/max(|a|,|b|) (0 when a==b):
* Stata's own reldif(), |a-b|/(|b|+1), is always <= |a-b|, which would make the
* prereg's "|diff| <= 1e-12 OR" arm redundant; the two-arm rule is only coherent with
* the magnitude-relative form, which also matches the sthlp's characterization of the
* documented fast-path divergence (~1e-15 abs, rel up to ~1e-13 on subgroup pairwise
* contrasts; wooldid.sthlp fastnose entry). Note: an entry that is ~0 in one leg and
* exactly 0 in the other can post worst_rel near 1 while PASSING via the abs arm —
* the FGRID worst_rel figure is informational, not itself a gate.
*
* Probe-established design notes (2026-08-24, scratch probes on this build):
*  - e(gteffects) keeps POPULATED se/p/lb/ub columns under fastnose: those SEs come
*    from the stage-1 regression's coefficient-level VCE, not from margins (which is
*    what fastnose skips), so all 15 columns are compared cell-by-cell.
*  - gteffects x contreat() is refused by design ("gteffects is not available with
*    contreat(): the cohort-by-period interaction coefficients are dose-zero
*    intercepts under a continuous treatment, not treatment effects."), so the
*    continuous-side gteffects config (c_gt) is dta_cont WITHOUT contreat.
*  - poisson x [aw=w] is refused (rc 198), and the poisson + ES-window full leg runs
*    ~5.3 minutes on dta_binary (the very margins cost fastnose exists to skip) —
*    far over this suite's runtime budget — so the second poisson config is
*    poisson + gteffects (~19s full leg) instead.
*
* Refusal contract (also preregistered): fastnose x each of {unconditionalse,
* adversarial, jointtests, pretrendpower(80), impvalues, permtstat, permtci} must be
* refused loudly. Every rc was pinned empirically at 198 on this build and is
* asserted exactly below. Each base call is valid without fastnose, so the fastnose
* gate is the screen that fires (probe-verified message texts):
*  - unconditionalse/impvalues/adversarial/jointtests/pretrendpower: the F9 gates in
*    wooldid.ado (~line 4225 block), each "Error: fastnose is not compatible with ...".
*    pretrendpower() itself requires esfixedbaseperiod (wooldid_ptpcheck screen), so
*    that call carries esfixedbaseperiod — without it the ptpcheck screen fires first
*    (also rc 198, but a different, non-fastnose message).
*  - permtstat/permtci: fires at parse time inside wooldid_permctl (~line 4440),
*    BEFORE the F9 block and before any permute mechanics or RNG use, with the
*    fastnose-specific message "Error: permtstat/permtci randomize on t-statistics,
*    which fastnose does not compute." The call needs permute(shuffle) permreps(20)
*    + cluster() for the combination to be reachable (permtstat/permtci without
*    permute() trip the earlier "require(s) the permute() option" screen instead).
*
* Args: 1 = path to wooldid.ado (file or dir); 2 = suite data dir (contains dta_binary.dta).
* Run:  "<stata>" -e do fastnose_grid.do <adopath> <datadir>
*    -> expect FASTNOSE_GRID_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 97531
set sortseed 97531

local ADO `1'
if `"`ADO'"' == "" local ADO /tmp/wdsuite_build/newfeat/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local DATA `2'
if `"`DATA'"' == "" local DATA /Users/hegland/Documents/Research/wooldid_dev/test/data

global pf 0

**************************** comparison machinery ****************************
* fgcmp: name-matched cell-by-cell comparison of two result matrices.
* Accumulates worst diffs in scalars __fgwabs/__fgwrel (reset per config by fgreset)
* and posts the violation count in __fgviol. At most 10 violating cells are printed
* per comparison (with a suppression note); the wrapper fgchk turns any violation
* into one $pf increment.
mata:
void fgcmp(string scalar cfg, string scalar what, string scalar mfull, string scalar mfast)
{
    real matrix A, B
    string matrix sf, sn
    string colvector fn, gn
    real colvector map, idx
    real scalar k, nr, r, j, a, b, ad, rd, viol, shown, wa, wr, nmiss

    A  = st_matrix(mfull)
    B  = st_matrix(mfast)
    sf = st_matrixcolstripe(mfull)
    sn = st_matrixcolstripe(mfast)
    fn = strtrim(sf[., 1]) :+ ":" :+ strtrim(sf[., 2])
    gn = strtrim(sn[., 1]) :+ ":" :+ strtrim(sn[., 2])
    viol = 0

    if (rows(A) != rows(B)) {
        printf("{err}FGRID %s FAIL: %s row count differs: full=%g fastnose=%g\n",
               cfg, what, rows(A), rows(B))
        st_numscalar("__fgviol", 1)
        return
    }
    if (cols(A) != cols(B)) {
        printf("{err}FGRID %s FAIL: %s column-name sets differ: full has %g cols, fastnose has %g\n",
               cfg, what, cols(A), cols(B))
        nmiss = 0
        for (j = 1; j <= cols(A); j++) {
            if (sum(gn :== fn[j]) == 0) {
                if (nmiss++ < 5) printf("{err}FGRID %s FAIL: %s column %s present in full leg only\n", cfg, what, fn[j])
            }
        }
        nmiss = 0
        for (j = 1; j <= cols(B); j++) {
            if (sum(fn :== gn[j]) == 0) {
                if (nmiss++ < 5) printf("{err}FGRID %s FAIL: %s column %s present in fastnose leg only\n", cfg, what, gn[j])
            }
        }
        st_numscalar("__fgviol", 1)
        return
    }
    k = cols(A)
    if (fn == gn) map = (1::k)      // ordered-identical names: 1:1 name alignment
    else {
        // orders differ: match by name, which needs unique names on both sides
        if (rows(uniqrows(fn)) != k | rows(uniqrows(gn)) != k) {
            printf("{err}FGRID %s FAIL: %s column names not unique and column orders differ; cannot match by name\n",
                   cfg, what)
            st_numscalar("__fgviol", 1)
            return
        }
        map = J(k, 1, 0)
        for (j = 1; j <= k; j++) {
            idx = selectindex(gn :== fn[j])
            if (rows(idx) != 1) {
                printf("{err}FGRID %s FAIL: %s column %s present in full leg but not fastnose leg\n",
                       cfg, what, fn[j])
                viol++
            }
            else map[j] = idx[1]
        }
        if (viol > 0) {
            st_numscalar("__fgviol", viol)
            return
        }
    }

    wa = st_numscalar("__fgwabs")
    wr = st_numscalar("__fgwrel")
    nr = rows(A)
    shown = 0
    for (j = 1; j <= k; j++) {
        for (r = 1; r <= nr; r++) {
            a = A[r, j]
            b = B[r, map[j]]
            if (a >= . & b >= .) continue                      // missing in both: equal
            if ((a >= .) != (b >= .)) {                        // missing-vs-nonmissing flip
                viol++
                if (shown++ < 10) {
                    printf("{err}FGRID %s FAIL: %s col %s row %g: missing flip: full=%21.15g fastnose=%21.15g\n",
                           cfg, what, fn[j], r, a, b)
                }
                continue
            }
            ad = abs(a - b)
            rd = (ad == 0 ? 0 : ad / max((abs(a), abs(b))))
            if (ad > wa) wa = ad
            if (rd > wr) wr = rd
            if (!(ad <= 1e-12 | rd <= 1e-10)) {
                viol++
                if (shown++ < 10) {
                    printf("{err}FGRID %s FAIL: %s col %s row %g: full=%21.15g fastnose=%21.15g absdiff=%10.3e reldiff=%10.3e\n",
                           cfg, what, fn[j], r, a, b, ad, rd)
                }
            }
        }
    }
    if (viol > 10) printf("{err}FGRID %s FAIL: %s: %g further violating cells suppressed\n", cfg, what, viol - 10)
    st_numscalar("__fgwabs", wa)
    st_numscalar("__fgwrel", wr)
    st_numscalar("__fgviol", viol)
}
end

cap prog drop fgreset
prog def fgreset
  scalar __fgwabs = 0
  scalar __fgwrel = 0
end

cap prog drop fgchk
prog def fgchk
  * fgchk <cfg> <eb|gteff> <fullmat> <fastmat>
  args cfg what mfull mfast
  local disp "e(b)"
  if "`what'" == "gteff" local disp "e(gteffects)"
  scalar __fgviol = 0
  mata: fgcmp("`cfg'", "`disp'", "`mfull'", "`mfast'")
  if scalar(__fgviol) > 0 {
    dis as error "FGRID `cfg' FAIL: `disp' comparison: " scalar(__fgviol) " violating entries (see lines above)"
    global pf = $pf + 1
  }
end

cap prog drop fgline
prog def fgline
  args cfg
  dis "FGRID `cfg' worst_abs=" %10.3e scalar(__fgwabs) " worst_rel=" %10.3e scalar(__fgwrel)
end

**************************** grid configs ****************************
* nm = config name; ds = dataset; gt = 1 if e(gteffects) is compared; cm = the wooldid
* call after "wooldid " (the fastnose leg appends " fastnose" at the end of the options).
* dta_cont gets sg = mod(i,2)+1 generated on load (runbattery.do precedent).
local nm1  b_att
local ds1  dta_binary
local gt1  0
local cm1  `"y i t ttre, att cluster(i)"'
local nm2  b_espool
local ds2  dta_binary
local gt2  0
local cm2  `"y i t ttre, att cluster(i) esprelength(3) espostlength(3)"'
local nm3  b_esf
local ds3  dta_binary
local gt3  0
local cm3  `"y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-2)"'
local nm4  b_antic
local ds4  dta_binary
local gt4  0
local cm4  `"y i t ttre, att cluster(i) anticipation(1) esprelength(2) espostlength(2)"'
local nm5  b_aw
local ds5  dta_binary
local gt5  0
local cm5  `"y i t ttre [aw=w], att cluster(i) esprelength(2) espostlength(2)"'
local nm6  b_pw
local ds6  dta_binary
local gt6  0
local cm6  `"y i t ttre [pw=w], att cluster(i) esprelength(2) espostlength(2)"'
local nm7  b_att_attit
local ds7  dta_binary
local gt7  0
local cm7  `"y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2)"'
local nm8  b_atti
local ds8  dta_binary
local gt8  0
local cm8  `"y i t ttre, att_i cluster(i)"'
local nm9  b_sg_noes
local ds9  dta_binary
local gt9  0
local cm9  `"y i t ttre, att cluster(i) subgroups(sg)"'
local nm10 b_sg_es
local ds10 dta_binary
local gt10 0
local cm10 `"y i t ttre, att cluster(i) subgroups(sg) esprelength(2) espostlength(2)"'
local nm11 b_sg_attit
local ds11 dta_binary
local gt11 0
local cm11 `"y i t ttre, att_it cluster(i) subgroups(sg) esprelength(2) espostlength(2)"'
local nm12 b_pois
local ds12 dta_binary
local gt12 0
local cm12 `"ycount i t ttre, att cluster(i) poisson"'
local nm13 b_pois_gt
local ds13 dta_binary
local gt13 1
local cm13 `"ycount i t ttre, att cluster(i) poisson gteffects"'
local nm14 b_attitime
local ds14 dta_binary
local gt14 0
local cm14 `"y i t ttre, att_itime cluster(i) esprelength(2) espostlength(2)"'
local nm15 b_sink
local ds15 dta_binary
local gt15 1
local cm15 `"y i t ttre [aw=w], att att_it cluster(i) subgroups(sg) esprelength(2) espostlength(2) esfixedbaseperiod esrelativeto(-2) anticipation(1) gteffects"'
local nm16 c_basic
local ds16 dta_cont
local gt16 0
local cm16 `"y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2)"'
local nm17 c_uz
local ds17 dta_cont
local gt17 0
local cm17 `"y i t ttre, contreat(intensity_uz) contreatcontrols(untreatedZero) safety(off) cluster(i) esprelength(2) espostlength(2)"'
local nm18 c_lattice
local ds18 dta_cont
local gt18 0
local cm18 `"y i t ttre [aw=w], contreat(intensity) contreatcontrols(basic) lattice(all 3) cluster(i) esprelength(2) espostlength(2)"'
local nm19 c_gt
local ds19 dta_cont
local gt19 1
local cm19 `"y i t ttre, att cluster(i) esprelength(2) espostlength(2) gteffects"'
local nm20 c_esf
local ds20 dta_cont
local gt20 0
local cm20 `"y i t ttre, contreat(intensity) contreatcontrols(basic) esfixedbaseperiod esprelength(3) espostlength(3) cluster(i)"'
local nm21 c_sg
local ds21 dta_cont
local gt21 0
local cm21 `"y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg)"'
local nm22 b_gt_aw
local ds22 dta_binary
local gt22 1
local cm22 `"y i t ttre [aw=w], att att_it cluster(i) esprelength(2) espostlength(2) gteffects"'
* c_lattice_robust (research_2026-09-08/fastnose_precision.md step 3): the c_lattice call with
* robust in place of cluster(i), so no containment screen and no all-zero-VCE refit can fire -
* lattice + contreat + aweight + fastnose end to end at the unrelaxed tolerance
local nm23 c_lattice_robust
local ds23 dta_cont
local gt23 0
local cm23 `"y i t ttre [aw=w], contreat(intensity) contreatcontrols(basic) lattice(all 3) robust esprelength(2) espostlength(2)"'
local NC 23

**************************** parity grid ****************************
forvalues k = 1/`NC' {
  local cfg `nm`k''

  * ---- full leg ----
  use `"`DATA'/`ds`k''.dta"', clear
  if "`ds`k''" == "dta_cont" qui gen sg = mod(i,2)+1
  cap noi wooldid `cm`k''
  if _rc != 0 {
    dis as error "FGRID `cfg' FAIL: full leg rc=" _rc " [wooldid `cm`k'']"
    global pf = $pf + 1
    continue
  }
  tempname BF GF BN GN
  matrix `BF' = e(b)
  if `gt`k'' == 1 matrix `GF' = e(gteffects)

  * ---- fastnose leg (fresh data load) ----
  use `"`DATA'/`ds`k''.dta"', clear
  if "`ds`k''" == "dta_cont" qui gen sg = mod(i,2)+1
  cap noi wooldid `cm`k'' fastnose
  if _rc != 0 {
    dis as error "FGRID `cfg' FAIL: fastnose leg rc=" _rc " [wooldid `cm`k'' fastnose]"
    global pf = $pf + 1
    cap matrix drop `BF'
    cap matrix drop `GF'
    continue
  }
  matrix `BN' = e(b)
  if `gt`k'' == 1 matrix `GN' = e(gteffects)

  * ---- compare ----
  fgreset
  fgchk `cfg' eb `BF' `BN'
  if `gt`k'' == 1 fgchk `cfg' gteff `GF' `GN'
  fgline `cfg'

  cap matrix drop `BF'
  cap matrix drop `BN'
  cap matrix drop `GF'
  cap matrix drop `GN'
}

**************************** refusal contract ****************************
* Every base call below runs rc 0 without fastnose (probed); with fastnose the pinned
* refusal is exactly rc 198 (all seven, pinned empirically on this build). rc 0 =
* not-refused FAIL; any other rc = wrong-screen FAIL.

cap prog drop fgref
prog def fgref
  * fgref <label> <observed rc>
  args lbl rc
  if `rc' == 0 {
    dis as error "FGRID refuse_`lbl' FAIL: fastnose + `lbl' was NOT refused (rc 0)"
    global pf = $pf + 1
  }
  else if `rc' != 198 {
    dis as error "FGRID refuse_`lbl' FAIL: expected rc 198, observed rc=" `rc'
    global pf = $pf + 1
  }
  else dis "FGRID refuse_`lbl' ok: rc=198"
end

use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre [aw=w], att robust unconditionalse esprelength(2) espostlength(2) fastnose
fgref unconditionalse `=_rc'

use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-2) adversarial fastnose
fgref adversarial `=_rc'

use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) jointtests fastnose
fgref jointtests `=_rc'

* esfixedbaseperiod is required by pretrendpower() itself; without it the ptpcheck
* screen fires first (also rc 198 but not the fastnose gate), so it is on the call
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod pretrendpower(80) fastnose
fgref pretrendpower `=_rc'

use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) impvalues fastnose
fgref impvalues `=_rc'

* permtstat/permtci: refusal fires at parse time in wooldid_permctl, before any
* permute mechanics/RNG, with the fastnose-specific message (see header notes)
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permtstat fastnose
fgref permtstat `=_rc'

use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permtci fastnose
fgref permtci `=_rc'

**************************** verdict ****************************
cap scalar drop __fgwabs __fgwrel __fgviol
if $pf == 0 dis "FASTNOSE_GRID_ALL_PASS"
else dis as error "FASTNOSE_GRID_FAILS=$pf"
exit, clear
