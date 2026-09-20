clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
include "$WD/shapelib.doi"
use "$WD/fixA.dta", clear

* generic scanner over any 2xK boot CI matrix
cap prog drop ckci
prog def ckci
  args tag mname bname
  cap confirm matrix `mname'
  if _rc {
    dis "CI|`tag'|`mname'|ABSENT"
    exit
  }
  mat CC = `mname'
  mat BO = `bname'
  local cn : colnames CC
  local nb = 0
  local nh = 0
  local no = 0
  local nc = 0
  forvalues k = 1/`=colsof(CC)' {
    local lo = CC[1,`k']
    local hi = CC[2,`k']
    local bb = BO[1,`k']
    local nm : word `k' of `cn'
    if `lo'<. & `hi'<. {
      local ++nc
      if !(`lo'<`hi') {
        local ++nb
        dis "VIOL_CILBUB|`tag'|`nm'|" `lo' "|" `hi'
      }
      if `bb'<. {
        if `bb'<`lo' | `bb'>`hi' {
          local ++no
          dis "VIOL_BOUT|`tag'|`nm'|b=" `bb' "|[" `lo' "," `hi' "]"
        }
      }
    }
    if (`lo'<. & `hi'>=.) | (`lo'>=. & `hi'<.) {
      local ++nh
      dis "VIOL_CIHALF|`tag'|`nm'|lb=" `lo' "|ub=" `hi'
    }
  }
  dis "CI|`tag'|`mname'|K=" colsof(CC) "|ncmp=`nc'|lbub=`nb'|half=`nh'|bout=`no'"
end

foreach L in 90 95 99 {
  cap noi wooldid y cg t ttre, att att_i cluster(unit) subgroups(sg) gteffects ///
      esprelength(2) espostlength(2) bootstrap(standard) reps(99) rseed(31337) ///
      level(`L') bootdisplay(bc) dots(0)
  dis "RC|J_`L'|" _rc
  if _rc==0 {
    dis "JLVL|`L'|e(level)=" e(level) "|completed=" e(boot_completed)
    wdck J_`L'_main e(wdidmainresults)
    wdck J_`L'_sg0  e(wdidsgresults_0)
    wdck J_`L'_sg1  e(wdidsgresults_1)
    wdck J_`L'_pw   e(pwsgs)
    foreach fam in normal percentile bc {
      ckci J_`L'_`fam' e(boot_ci_`fam') e(b)
    }
    * gteffects cells
    mat GT = e(gteffects)
    local cn : colnames GT
    local ge=0
    local gl=0
    local gu=0
    local j=0
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
    }
    local nv=0
    local nc=0
    local nh=0
    forvalues i=1/`=rowsof(GT)' {
      local lo=GT[`i',`gl']
      local hi=GT[`i',`gu']
      local b =GT[`i',`ge']
      if `lo'<. & `hi'<. {
        local ++nc
        if !(`lo'<`hi') {
          local ++nv
          dis "VIOL_GTLBUB|J_`L'|row`i'|" `lo' "|" `hi'
        }
        if `b'<. {
          if `b'<`lo'-1e-11 | `b'>`hi'+1e-11 {
            local ++nv
            dis "VIOL_GTBIN|J_`L'|row`i'|b=" `b' "|[" `lo' "," `hi' "]"
          }
        }
      }
      if (`lo'<. & `hi'>=.) | (`lo'>=. & `hi'<.) {
        local ++nh
        dis "VIOL_GTHALF|J_`L'|row`i'|lb=" `lo' "|ub=" `hi'
      }
    }
    dis "GT|J_`L'|rows=" rowsof(GT) "|ncmp=`nc'|viol=`nv'|half=`nh'|cilab=" `"`: colnames GT'"'
    if `L'==90 {
      mat K90 = e(wdidmainresults)
      mat S90 = e(wdidsgresults_0)
      mat W90 = e(pwsgs)
    }
    if `L'==99 {
      mat K99 = e(wdidmainresults)
      mat S99 = e(wdidsgresults_0)
      mat W99 = e(pwsgs)
    }
  }
}
wdnest J_main_90in99 K90 K99
wdnest J_sg0_90in99  S90 S99
wdnest J_pw_90in99   W90 W99
dis "PROBEJ_DONE"
