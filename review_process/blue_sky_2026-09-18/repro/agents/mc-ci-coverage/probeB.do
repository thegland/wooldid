clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
set seed 20260916
set sortseed 20260916
include "$WD/shapelib.doi"
use "$WD/fixA.dta", clear

* ---------- B0: nesting on the analytic path ----------
cap noi wooldid y cg t ttre, att att_i att_it att_itime cluster(unit) esprelength(3) espostlength(3) level(90)
mat A90 = e(wdidmainresults)
cap noi wooldid y cg t ttre, att att_i att_it att_itime cluster(unit) esprelength(3) espostlength(3) level(99)
mat A99 = e(wdidmainresults)
wdnest A_90in99 A90 A99

* ---------- B1: simultaneous bands, analytic engine ----------
foreach L in 90 95 99 {
  cap noi wooldid y cg t ttre, att att_i cluster(unit) esprelength(3) espostlength(3) level(`L') simultaneousbands
  dis "RC|B1_`L'|" _rc
  if _rc == 0 {
    dis "BANDENG|B1_`L'|" e(esband_engine) "|fams=" e(esband_families) "|ndraw=" e(esband_ndraw)
    foreach f in att att_i {
      cap confirm matrix e(esband_`f')
      if _rc == 0 {
        mat BB = e(esband_`f')
        mat BV = e(V_es_`f')
        mat Bb = e(b_es_`f')
        local crit = e(esband_crit_`f')
        dis "BANDCRIT|B1_`L'|`f'|crit=" %12.9f `crit' "|znorm=" %12.9f invnormal(1-(100-`L')/200) "|tG=" %12.9f invttail(119,(100-`L')/200)
        local H = colsof(BB)
        local bad = 0
        forvalues h = 1/`H' {
          local lo = BB[1,`h']
          local hi = BB[2,`h']
          local bh = Bb[1,`h']
          local sh = sqrt(BV[`h',`h'])
          local wband = `hi'-`lo'
          local wpt   = 2*invttail(119,(100-`L')/200)*`sh'
          if `wband' < `wpt' - 1e-10 {
            local ++bad
            dis "VIOL_BANDNARROW|B1_`L'|`f'|h=`h'|wband=" `wband' "|wpt=" `wpt'
          }
          if !(`lo' < `hi') {
            dis "VIOL_BANDLBUB|B1_`L'|`f'|h=`h'|lo=" `lo' "|hi=" `hi'
          }
          if `bh' < `lo' | `bh' > `hi' {
            dis "VIOL_BANDBIN|B1_`L'|`f'|h=`h'|b=" `bh' "|lo=" `lo' "|hi=" `hi'
          }
        }
        dis "BAND|B1_`L'|`f'|H=`H'|narrow=`bad'"
      }
    }
  }
}

* ---------- B2: subgroups / unconditionalse / impvalues ----------
cap noi wooldid y cg t ttre, att att_i cluster(unit) esprelength(2) espostlength(2) subgroups(sg) level(95)
dis "RC|B2_sg|" _rc
if _rc==0 {
  wdck B2_sg_main e(wdidmainresults)
  wdck B2_sg0 e(wdidsgresults_0)
  wdck B2_sg1 e(wdidsgresults_1)
  wdck B2_pw e(pwsgs)
}
cap noi wooldid y cg t ttre, att att_i unconditionalse level(95) esprelength(2) espostlength(2)
dis "RC|B2_unc|" _rc
if _rc==0 {
  wdck B2_unc e(wdidmainresults)
}
cap noi wooldid y cg t ttre, att att_i impvalues level(95) esprelength(2) espostlength(2)
dis "RC|B2_imp|" _rc
if _rc==0 {
  wdck B2_imp e(wdidmainresults)
  dis "IMPDF|" e(df_r)
  mat list e(wdidmainresults)
}
cap noi wooldid y cg t ttre, att att_i impvalues level(90) esprelength(2) espostlength(2)
if _rc==0 mat I90 = e(wdidmainresults)
cap noi wooldid y cg t ttre, att att_i impvalues level(99) esprelength(2) espostlength(2)
if _rc==0 mat I99 = e(wdidmainresults)
wdnest IMP_90in99 I90 I99
cap noi wooldid y cg t ttre, att att_i impvalues level(87.5) esprelength(2) espostlength(2)
if _rc==0 {
  dis "IMPLAB|" `"`: colnames e(wdidmainresults)'"'
}

* ---------- B3: gteffects cells ----------
cap noi wooldid y cg t ttre, att cluster(unit) gteffects level(95)
dis "RC|B3_gt|" _rc
if _rc==0 {
  mat GT = e(gteffects)
  local cn : colnames GT
  dis "GTCOLS|`cn'"
  local j=0
  local ge=0
  local gl=0
  local gu=0
  foreach c of local cn {
    local ++j
    if substr("`c'",1,3)=="lb_" {
      local gl=`j'
    }
    if substr("`c'",1,3)=="ub_" {
      local gu=`j'
    }
    if "`c'"=="estimate" {
      local ge=`j'
    }
    if "`c'"=="b" & `ge'==0 {
      local ge=`j'
    }
  }
  dis "GTIDX|e=`ge'|l=`gl'|u=`gu'"
  local nv=0
  local ncmp=0
  forvalues i=1/`=rowsof(GT)' {
    local lo=GT[`i',`gl']
    local hi=GT[`i',`gu']
    local b =GT[`i',`ge']
    if `lo'<. & `hi'<. {
      local ++ncmp
      if !(`lo'<`hi') {
        local ++nv
        dis "VIOL_GTLBUB|row`i'|" `lo' "|" `hi'
      }
      if `b'<. {
        if (`b'<`lo'-1e-11 | `b'>`hi'+1e-11) {
          local ++nv
          dis "VIOL_GTBIN|row`i'|b=" `b' "|" `lo' "|" `hi'
        }
      }
    }
  }
  dis "GT|rows=" rowsof(GT) "|ncmp=`ncmp'|viol=`nv'"
}
dis "PROBEB_DONE"
