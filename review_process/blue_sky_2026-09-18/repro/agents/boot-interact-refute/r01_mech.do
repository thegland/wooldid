* r01_mech.do -- SELF-CONTAINED. BOOTX-1 mechanism + reps ladder, my own data.
* Q: is the printed "joint pre-trend test" under bootstrap() literally
*    Ftail(keff, e(df_r), W/keff) with W = b'S^-1 b on the BOOTSTRAP S, and is
*    e(df_r) = G-1 regardless of reps()?  And how far is it from Hotelling F(k,n-k)?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

* ---------- generator: 30 states x 2 counties, 16 periods, NO pre-trend in the DGP
clear
set obs 30
gen state = _n
gen onset = .
replace onset = 7 + mod(_n, 6) if mod(_n,3)==0
expand 2
bysort state: gen cty = _n
gen county = state*10 + cty
expand 16
bysort county: gen t = _n
bysort county (t): gen double ushock = rnormal(0,.3) if _n==1
bysort county (t): replace ushock = ushock[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .12*rel, 0)
gen double y = 2 + .05*t + ushock + te + e
drop rel te

* ---------- reporter: recompute the ptp Wald by hand and give both references
cap prog drop rpt
prog def rpt
  args tag
  tempname b V S bp Si W
  matrix `b' = e(b_es_att)
  matrix `V' = e(V_es_att)
  local ln : colnames `b'
  local k 0
  local idx
  local j 0
  foreach l of local ln {
    local ++j
    if substr("`l'",1,3)=="pre" {
      local ++k
      local idx `idx' `j'
    }
  }
  matrix `S'  = J(`k',`k',0)
  matrix `bp' = J(`k',1,0)
  local r 0
  foreach i1 of local idx {
    local ++r
    matrix `bp'[`r',1] = `b'[1,`i1']
    local c 0
    foreach i2 of local idx {
      local ++c
      matrix `S'[`r',`c'] = `V'[`i1',`i2']
    }
  }
  matrix `S' = (`S' + `S'')/2
  matrix `Si' = invsym(`S')
  local keff = `k' - diag0cnt(`Si')
  matrix `W' = `bp''*`Si'*`bp'
  local Wh = `W'[1,1]
  local dfr = e(df_r)
  local nb  = e(boot_completed)
  if "`e(boot_method)'"=="" local nb = .
  local phand = .
  if `dfr' < . local phand = Ftail(`keff',`dfr',`Wh'/`keff')
  else local phand = chi2tail(`keff',`Wh')
  local phot = .
  if `nb' < . & `nb' > `keff' local phot = Ftail(`keff', `nb'-`keff', (`nb'-`keff')/((`nb'-1)*`keff')*`Wh')
  dis "RES|`tag'|k=`k'|keff=`keff'|W=" %12.6g `Wh' "|dfr=`dfr'|ncomp=`nb'|NV=" e(boot_N_V) ///
      "|p_printed=" %12.6g e(roth_prep_att) "|p_hand=" %12.6g `phand' "|p_hotelling=" %12.6g `phot' ///
      "|Wstored=" %12.6g e(roth_W_att) "|keffstored=" e(roth_keff_att) ///
      "|slope=" %10.6g e(roth_slope_att) "|se=" %10.6g e(roth_slope_se_att) ///
      "|ll=" %10.6g e(roth_slope_ll_att) "|ul=" %10.6g e(roth_slope_ul_att) ///
      "|mdt80=" %10.6g e(roth_power_att)[1,2] "|rankjoint=" %10.6g e(joint_horiz)
end

local OPT "cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) jointtests pretrendpower(80)"

qui cap noi wooldid y county t onset, `OPT'
dis "A_RC=" _rc
rpt ANALYTIC

foreach R in 20 50 100 300 {
  qui cap noi wooldid y county t onset, `OPT' bootstrap(standard) reps(`R') rseed(777)
  dis "B`R'_RC=" _rc
  rpt BOOT`R'
}

* bootdf variants at fixed reps to show e(df_r) is the only driver
foreach D in clusters normal 3 {
  qui cap noi wooldid y county t onset, `OPT' bootstrap(standard) reps(50) rseed(777) bootdf(`D')
  dis "D`D'_RC=" _rc
  rpt BOOTDF_`D'
}
dis "DONE"
