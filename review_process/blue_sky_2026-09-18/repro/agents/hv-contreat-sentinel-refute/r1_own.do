* =============================================================================
* R1 -- SKEPTIC's own repro, built from the CLAIMS, on MY OWN data.
*
* Claims under test (from the prober):
*   hv-F5-verify        : b=0 + e(V) diag -1e100 on pretreatment ame rows, while
*                         e(V_es_ame) carries EXACTLY 0 on the same horizons.
*   hv-adv-ame-rank     : the adversarial ame statistic drops the pre block
*                         entirely (bitwise equal to the post-only statistic).
*   hv-noVsentinel-mark : no e() scalar records the e(V) sentinel count.
*
* MY DGP (deliberately different from the prober's): hospital-year administrative
* panel.  i = STATE (12), 30 hospitals per state, T = 2008..2019.  Onset 2012 /
* 2014 / 2016, three never-expanding states.  dose = hospital-specific exposure,
* switched on at that state's expansion year (the canonical exposure x post dose).
*
* Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r1_own.do
* =============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333

dis "=== WHICH (must be the pinned build) ==="
which wooldid

* ------------------------------------------------------------------ DGP ----
clear
set obs 12
gen byte state = _n
gen int expyr = .
replace expyr = 2012 if inlist(state,1,2,3,4)
replace expyr = 2014 if inlist(state,5,6,7)
replace expyr = 2016 if inlist(state,8,9)
* states 10,11,12 never expand
expand 30
bys state: gen int hosp = _n
gen long hid = state*100 + hosp
* hospital-level exposure (share of uninsured discharges), lognormal-ish
gen double expo = 0.05 + 0.35*rbeta(2,3)
gen double hfe = rnormal(0,1.2)
expand 12
bys hid: gen int year = 2007 + _n
gen double dose = expo*(expyr<. & year>=expyr)
gen double yfe = 0.03*(year-2007)
gen double y = 4 + hfe + yfe + 0.8*dose + rnormal(0,0.9)
label var y "log discharges"
compress
save "r1_panel.dta", replace
count
tab dose if dose==0 & expyr<. & year<expyr, missing

* ---------------------------------------------------------- dump program ---
cap prog drop dumpe
prog def dumpe
  syntax , tag(string)
  tempname B V
  cap confirm matrix e(b)
  if _rc {
    dis "DUMP_`tag'_NOB=1"
    exit
  }
  matrix `B' = e(b)
  local bn : colfullnames `B'
  local hasV = 0
  cap confirm matrix e(V)
  if !_rc {
    local hasV 1
    matrix `V' = e(V)
  }
  dis "DUMP_`tag'_K=" colsof(`B') "  HASV=" `hasV'
  local nsent 0
  local j 0
  foreach nm of local bn {
    local ++j
    local bj = `B'[1,`j']
    local vj = .
    if `hasV' local vj = `V'[`j',`j']
    local mk ""
    if `hasV' & `vj' < 0 & !missing(`vj') {
      local mk "SENTINEL"
      local ++nsent
    }
    dis "ROW_`tag'|" %-26s "`nm'" "|b=" %12.0g `bj' "|Vd=" %12.0g `vj' "|`mk'"
  }
  dis "NSENT_`tag'=" `nsent'
  foreach f in att ame att_it ame_it {
    cap confirm matrix e(V_es_`f')
    if !_rc {
      tempname E
      matrix `E' = e(V_es_`f')
      local en : colfullnames `E'
      local nz 0
      local jj 0
      foreach e1 of local en {
        local ++jj
        local d = `E'[`jj',`jj']
        if `d' == 0 local ++nz
        dis "VES_`tag'_`f'|" %-10s "`e1'" "|d=" %14.0g `d'
      }
      dis "VESZERO_`tag'_`f'=" `nz'
    }
    else dis "VES_`tag'_`f'_ABSENT=1"
    cap confirm matrix e(b_es_`f')
    if !_rc dis "BES_`tag'_`f'=1"
    else dis "BES_`tag'_`f'=0"
    dis "NDEGEN_`tag'_`f'=" e(es_ndegen_`f')
    dis "NUMPRE_`tag'_`f'=" e(numpre_`f')
  }
  dis "VCVSTATE_`tag'=" e(wdid_vcvstate)
  dis "DFR_`tag'=" e(df_r)
  dis "ADVP_att_`tag'=" e(adversarial_pval_att)
  dis "ADVP_ame_`tag'=" e(adversarial_pval_ame)
end

* ================================ CONFIGS ==================================
foreach CFG in IID ROB CLHOSP {
  use "r1_panel.dta", clear
  local vopt ""
  if "`CFG'" == "ROB"    local vopt customvce(robust)
  if "`CFG'" == "CLHOSP" local vopt cluster(hid)
  dis "############### CFG=`CFG' ###############"
  cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
      contreat(dose) contreatcontrols(untreatedZero) adversarial jointtests `vopt'
  dis "RC_`CFG'=" _rc
  if _rc continue
  dumpe , tag(`CFG')
  if "`CFG'" == "CLHOSP" estimates save "r1_clhosp", replace
  if "`CFG'" == "IID"    estimates save "r1_iid", replace

  * ---- adversarial: use WOOLDID'S OWN kernel on full vs post-only ----------
  cap confirm matrix e(V_es_ame)
  if !_rc {
    matrix BA = e(b_es_ame)
    matrix VA = e(V_es_ame)
    scalar DFR = e(df_r)
    scalar PAME = e(adversarial_pval_ame)
    local kA = colsof(BA)
    dis "ADVCHK_`CFG'_kA=" `kA'
    * relative-period vector implied by the column names pre3 pre2 pre1 post0..
    matrix REL = J(`kA',1,.)
    local cn : colfullnames BA
    local jj 0
    foreach c1 of local cn {
      local ++jj
      local lab = subinstr("`c1'","b:","",.)
      local lab = subinstr("`lab'",":","",.)
      if strpos("`lab'","pre") {
        local num = substr("`lab'", strpos("`lab'","pre")+3, .)
        matrix REL[`jj',1] = -1*real("`num'")
      }
      else if strpos("`lab'","post") {
        local num = substr("`lab'", strpos("`lab'","post")+4, .)
        matrix REL[`jj',1] = real("`num'")
      }
      dis "RELMAP_`CFG'|" %-10s "`lab'" "|rel=" REL[`jj',1]
    }
    * post-only submatrix
    mata: st_local("npre", strofreal(sum(st_matrix("REL") :< 0)))
    local p1 = `npre'+1
    matrix BAP = BA[1, `p1'..`kA']
    matrix VAP = VA[`p1'..`kA', `p1'..`kA']
    matrix RELP = REL[`p1'..`kA', 1]
    cap noisily mata: qf=.; Wf = wooldid_advW(st_matrix("BA")', st_matrix("VA"), st_matrix("REL"), -1, 0, 0, qf); st_numscalar("Wf",Wf); st_numscalar("qf",qf)
    dis "ADVKERNEL_RC_`CFG'=" _rc
    cap noisily mata: qp=.; Wp = wooldid_advW(st_matrix("BAP")', st_matrix("VAP"), st_matrix("RELP"), -1, 0, 0, qp); st_numscalar("Wp",Wp); st_numscalar("qp",qp)
    dis "ADV_`CFG'_Wfull=" %21x Wf "  q=" qf
    dis "ADV_`CFG'_Wpost=" %21x Wp "  q=" qp
    dis "ADV_`CFG'_BITEQ=" (Wf==Wp & qf==qp)
    dis "ADV_`CFG'_p_from_full=" %12.10f Ftail(qf, DFR, Wf/qf)
    dis "ADV_`CFG'_p_from_post=" %12.10f Ftail(qp, DFR, Wp/qp)
    dis "ADV_`CFG'_p_printed="   %12.10f PAME
    * att comparator
    matrix BT = e(b_es_att)
    matrix VT = e(V_es_att)
    cap noisily mata: qt=.; Wt = wooldid_advW(st_matrix("BT")', st_matrix("VT"), st_matrix("REL"), -1, 0, 0, qt); st_numscalar("Wt",Wt); st_numscalar("qt",qt)
    dis "ADV_`CFG'_att_q=" qt "  p=" %12.10f Ftail(qt, DFR, Wt/qt) "  printed=" %12.10f e(adversarial_pval_att)
  }

  * ---- hv-noVsentinel-marker: does ANY e() scalar record the sentinel count?
  if "`CFG'" == "CLHOSP" {
    dis "=== ereturn list, all (CLHOSP) ==="
    ereturn list
  }
}

* ================ contreatcontrols(basic) cross-check ======================
use "r1_panel.dta", clear
dis "############### CFG=BASIC ###############"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(basic) adversarial jointtests customvce(robust)
dis "RC_BASIC=" _rc
if _rc == 0 dumpe , tag(BASIC)

* ============ realism control: dose defined in EVERY period ================
use "r1_panel.dta", clear
replace dose = expo
dis "############### CFG=ALLPERIOD (time-invariant exposure) ###############"
cap noisily wooldid y state year expyr, att esfixedbaseperiod esprelength(3) espostlength(3) ///
    contreat(dose) contreatcontrols(basic) adversarial jointtests customvce(robust)
dis "RC_ALLPER=" _rc
if _rc == 0 dumpe , tag(ALLPER)

dis "R1_DONE=1"
