clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
set seed 20260916
set sortseed 20260916
include "$WD/shapelib.doi"
which wooldid

local FN 120
local FT 8
local FSEED 7
include "$WD/mkfix.doi"
save "$WD/fixA.dta", replace

* ---------- A1: analytic, four levels, ES rows, subgroups ----------
foreach L in 90 95 99 87.5 {
  local Ltag = subinstr("`L'",".","p",.)
  cap noi wooldid y cg t ttre, att att_i att_it att_itime cluster(unit) ///
        esprelength(3) espostlength(3) level(`L')
  dis "RC|A1_`Ltag'|" _rc
  if _rc == 0 {
    dis "ELEVEL|A1_`Ltag'|" e(level)
    wdck A1_`Ltag' e(wdidmainresults)
    tempname M`Ltag'
    mat `M`Ltag'' = e(wdidmainresults)
    * implied critical value
    local cn : colnames e(wdidmainresults)
    mat WW = e(wdidmainresults)
    local nr = rowsof(WW)
    forvalues i = 1/`nr' {
      local b = WW[`i',1]
      local s = WW[`i',2]
      local lo = WW[`i',5]
      local hi = WW[`i',6]
      if `s' < . & `s' > 0 & `lo' < . {
        local z = (`hi'-`lo')/(2*`s')
        local zl = (`b'-`lo')/`s'
        local zu = (`hi'-`b')/`s'
        local rn : rownames WW
        local rr : word `i' of `rn'
        dis "CRIT|A1_`Ltag'|`rr'|z=" %12.9f `z' "|zl=" %12.9f `zl' "|zu=" %12.9f `zu'
      }
    }
    dis "ZREF|A1_`Ltag'|norm=" %12.9f invnormal(1-(100-`L')/200) "|t119=" %12.9f invttail(119,(100-`L')/200) "|t3=" %12.9f invttail(3,(100-`L')/200)
  }
}
* nesting: 90 inside 95 inside 99
cap noi wooldid y cg t ttre, att att_i att_it att_itime cluster(unit) esprelength(3) espostlength(3) level(90)
mat A90 = e(wdidmainresults)
cap noi wooldid y cg t ttre, att att_i att_it att_itime cluster(unit) esprelength(3) espostlength(3) level(99)
mat A99 = e(wdidmainresults)
wdnest A_90in99 A90 A99

* ---------- A2: simultaneous bands, analytic engine ----------
foreach L in 90 95 99 {
  cap noi wooldid y cg t ttre, att att_i cluster(unit) esprelength(3) espostlength(3) ///
        level(`L') simultaneousbands
  dis "RC|A2_`L'|" _rc
  if _rc == 0 {
    dis "BANDENG|A2_`L'|" e(esband_engine) "|fams=" e(esband_families) "|lvl=" e(esband_level)
    foreach f in att att_i {
      cap confirm matrix e(esband_`f')
      if !_rc {
        mat BB = e(esband_`f')
        mat BV = e(V_es_`f')
        mat Bb = e(b_es_`f')
        local crit = e(esband_crit_`f')
        dis "BANDCRIT|A2_`L'|`f'|crit=" %12.9f `crit' "|znorm=" %12.9f invnormal(1-(100-`L')/200)
        local H = colsof(BB)
        local bad = 0
        forvalues h = 1/`H' {
          local lo = BB[1,`h']
          local hi = BB[2,`h']
          local bh = Bb[1,`h']
          local sh = sqrt(BV[`h',`h'])
          * pointwise analytic interval at same level
          local plo = `bh' - invnormal(1-(100-`L')/200)*`sh'
          local phi = `bh' + invnormal(1-(100-`L')/200)*`sh'
          local wband = `hi'-`lo'
          local wpt   = `phi'-`plo'
          if `wband' < `wpt' - 1e-10 {
            local ++bad
            dis "VIOL_BANDNARROW|A2_`L'|`f'|h=`h'|wband=" `wband' "|wpt=" `wpt'
          }
          if !(`lo' < `hi') {
            dis "VIOL_BANDLBUB|A2_`L'|`f'|h=`h'|lo=" `lo' "|hi=" `hi'
          }
          if `bh' < `lo' | `bh' > `hi' {
            dis "VIOL_BANDBIN|A2_`L'|`f'|h=`h'|b=" `bh' "|lo=" `lo' "|hi=" `hi'
          }
        }
        dis "BAND|A2_`L'|`f'|H=`H'|narrow=`bad'"
      }
    }
  }
}

* ---------- A3: subgroups + unconditionalse + impvalues ----------
cap noi wooldid y cg t ttre, att att_i cluster(unit) esprelength(2) espostlength(2) subgroups(sg) level(95)
dis "RC|A3_sg|" _rc
if _rc==0 {
  wdck A3_sg_main e(wdidmainresults)
  foreach v in 0 1 {
    wdck A3_sg`v' e(wdidsgresults_`v')
  }
  cap confirm matrix e(pwsgs)
  if !_rc wdck A3_pw e(pwsgs)
}
cap noi wooldid y cg t ttre, att att_i unconditionalse level(95) esprelength(2) espostlength(2)
dis "RC|A3_unc|" _rc
if _rc==0 wdck A3_unc e(wdidmainresults)
cap noi wooldid y cg t ttre, att att_i impvalues level(95) esprelength(2) espostlength(2)
dis "RC|A3_imp|" _rc
if _rc==0 {
  wdck A3_imp e(wdidmainresults)
  dis "IMPDF|" e(df_r)
}
cap noi wooldid y cg t ttre, att att_i impvalues level(90) esprelength(2) espostlength(2)
if _rc==0 mat I90 = e(wdidmainresults)
cap noi wooldid y cg t ttre, att att_i impvalues level(99) esprelength(2) espostlength(2)
if _rc==0 mat I99 = e(wdidmainresults)
wdnest IMP_90in99 I90 I99

* ---------- A4: gteffects cells ----------
cap noi wooldid y cg t ttre, att cluster(unit) gteffects level(95)
dis "RC|A4_gt|" _rc
if _rc==0 {
  mat GT = e(gteffects)
  local cn : colnames GT
  dis "GTCOLS|`cn'"
  local j=0
  foreach c of local cn {
    local ++j
    if substr("`c'",1,3)=="lb_" local gl=`j'
    if substr("`c'",1,3)=="ub_" local gu=`j'
    if "`c'"=="estimate" local ge=`j'
    if "`c'"=="b" local ge=`j'
  }
  dis "GTIDX|e=`ge'|l=`gl'|u=`gu'"
  local nv=0
  forvalues i=1/`=rowsof(GT)' {
    local lo=GT[`i',`gl']
    local hi=GT[`i',`gu']
    local b =GT[`i',`ge']
    if `lo'<. & `hi'<. {
      if !(`lo'<`hi') { local ++nv 
        dis "VIOL_GTLBUB|row`i'|" `lo' "|" `hi' }
      if `b'<. & (`b'<`lo'-1e-11 | `b'>`hi'+1e-11) { local ++nv 
        dis "VIOL_GTBIN|row`i'|b=" `b' "|" `lo' "|" `hi' }
    }
  }
  dis "GT|rows=" rowsof(GT) "|viol=`nv'"
}
dis "PROBEA_DONE"
