*==========================================================================
* REPRO: the simultaneous band's critical value never sees bootdf(), so with
*        the DEFAULT bootdf(clusters) and a modest cluster count the 95%
*        "simultaneous" band is NARROWER than the 95% POINTWISE interval at
*        EVERY horizon -- on bootdisplay(wald), where the pointwise interval
*        is symmetric and the help's documented non-nesting caveat does not
*        apply.
*
*  Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_band_ignores_bootdf.do
*
*  pointwise (bootdisplay(wald)) : b +/- invt(G-1, .975) * sd     (ado:9576, 9707)
*  simultaneous band             : b +/- q_.95(max_h |dev_h|/sd_h) * sd
*                                  (wooldid_bootband, ado:16876 -- lev only, no df)
*  At G = 8 : invt(7,.975) = 2.3646 while the band's critical value is 2.1089.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

clear
set obs 8
gen int state = _n
gen int onset = .
replace onset = 2003 in 1
replace onset = 2003 in 2
replace onset = 2004 in 3
replace onset = 2004 in 4
replace onset = 2005 in 5
replace onset = 2005 in 6
gen double a_i = rnormal()
expand 3
bys state: gen int unit = _n + 1000*state
gen double u_i = rnormal()*0.5
expand 8
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t - onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset

wooldid y state t ttre , cluster(state) bootstrap(standard) reps(999) rseed(77) ///
   bootdisplay(wald) att esprelength(1) espostlength(1) simultaneousbands dots(0)

local df = e(boot_df)
local cc = e(esband_crit_att)
dis _n "e(boot_df)                 = `df'      (bootdf(clusters) is the default with cluster())"
dis    "pointwise critical value   = " %9.6f invt(`df', 0.975)
dis    "band critical value        = " %9.6f `cc'
dis    "band critical value is SMALLER = " (`cc' < invt(`df',0.975))

mat MR = e(wdidmainresults)
mat SB = e(esband_att)
mata:
  MR = st_matrix("MR"); rs = st_matrixrowstripe("MR")
  SB = st_matrix("SB"); cs = st_matrixcolstripe("SB")
  nv = 0
  printf("\n%-8s %14s %14s %14s %14s   %s\n", "horizon","point_lb","point_ub","band_lb","band_ub","verdict")
  for (j = 1; j <= cols(SB); j++) {
    hit = 0
    for (r = 1; r <= rows(rs); r++) if (rs[r,1] == "ES_att" & rs[r,2] == cs[j,2]) hit = r
    if (hit == 0) continue
    plb = MR[hit,5]; pub = MR[hit,6]; blb = SB[1,j]; bub = SB[2,j]
    bad = (blb > plb + 1e-9) | (bub < pub - 1e-9)
    nv = nv + bad
    printf("%-8s %14.8f %14.8f %14.8f %14.8f   %s\n", cs[j,2], plb, pub, blb, bub,
           (bad ? "BAND NARROWER THAN POINTWISE" : "ok"))
  }
  printf("\nHORIZONS_WITH_BAND_INSIDE_POINTWISE = %f\nFAIL = %f\n", nv, (nv > 0))
end
