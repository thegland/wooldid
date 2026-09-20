* ================================================================
* ptplib.do -- programs shared by the pretrendpower() Monte Carlos.
* `ptpstat' recomputes, from a posted wooldid estimate, EXACTLY the
* quantities wooldid_ptp computes (verified to 1-2 ulp against
* e(roth_*) in p02_algebra.log), parsing the event-study column
* names the way ado:10450-10462 does so it is not tied to one design.
* ================================================================
capture program drop ptpstat
program define ptpstat, rclass
  syntax , [ref(real -1) type(string) p1(real 0.50) p2(real 0.80)]
  if "`type'" == "" local type att
  tempname B V S bp mm Si Q1 Q2 Q3
  capture confirm matrix e(b_es_`type')
  if _rc != 0 {
    return scalar ok = 0
    exit
  }
  capture confirm matrix e(V_es_`type')
  if _rc != 0 {
    return scalar ok = 0
    exit
  }
  matrix `B' = e(b_es_`type')
  matrix `V' = e(V_es_`type')
  local dfr = e(df_r)
  local cn : colnames `B'
  local k 0
  local idx
  local mv
  local j 0
  foreach c of local cn {
    local ++j
    if substr("`c'", 1, 3) == "pre" {
      local ++k
      local idx `idx' `j'
      local mv `mv' `= -real(substr("`c'", 4, .)) - (`ref')'
    }
  }
  if `k' == 0 {
    return scalar ok = 0
    exit
  }
  matrix `S'  = J(`k', `k', 0)
  matrix `bp' = J(`k', 1, 0)
  matrix `mm' = J(`k', 1, 0)
  local r 0
  foreach a of local idx {
    local ++r
    matrix `bp'[`r', 1] = `B'[1, `a']
    local mvv : word `r' of `mv'
    matrix `mm'[`r', 1] = `mvv'
    local c 0
    foreach b of local idx {
      local ++c
      matrix `S'[`r', `c'] = `V'[`a', `b']
    }
  }
  matrix `S' = (`S' + `S'')/2
  matrix `Si' = invsym(`S')
  local keff = `k' - diag0cnt(`Si')
  matrix `Q1' = `mm''*`Si'*`mm'
  matrix `Q2' = `mm''*`Si'*`bp'
  matrix `Q3' = `bp''*`Si'*`bp'
  local c2 = `Q1'[1, 1]
  local dfok = 0
  if `dfr' > 0 & `dfr' < . local dfok = 1
  return scalar ok    = 1
  return scalar k     = `k'
  return scalar keff  = `keff'
  return scalar dfr   = `dfr'
  return scalar c2    = `c2'
  return scalar slope = `Q2'[1,1]/`c2'
  return scalar se    = 1/sqrt(`c2')
  return scalar W     = `Q3'[1,1]
  if `dfok' {
    return scalar p    = Ftail(`keff', `dfr', `Q3'[1,1]/`keff')
    return scalar mdt1 = sqrt(npnF(`keff', `dfr', invF(`keff', `dfr', .95), 1-`p1')/`c2')
    return scalar mdt2 = sqrt(npnF(`keff', `dfr', invF(`keff', `dfr', .95), 1-`p2')/`c2')
    local tc = invttail(`dfr', .025)
    return scalar crit = `tc'
  }
  else {
    return scalar p    = chi2tail(`keff', `Q3'[1,1])
    return scalar mdt1 = sqrt(npnchi2(`keff', invchi2(`keff', .95), 1-`p1')/`c2')
    return scalar mdt2 = sqrt(npnchi2(`keff', invchi2(`keff', .95), 1-`p2')/`c2')
    local tc = invnormal(.975)
    return scalar crit = `tc'
  }
  local nis = 0
  forvalues j = 1/`k' {
    if abs(`bp'[`j',1])/sqrt(`S'[`j',`j']) > `tc' local nis = 1
  }
  return scalar nis = `nis'
  return matrix bp = `bp'
  return matrix S  = `S'
  return matrix mm = `mm'
end

* one Monte-Carlo arm.  `gen' is the name of a program that (re)draws y
* given globals DELTA and REP; `est' is the wooldid option string.
capture program drop runarm
program define runarm
  syntax , delta(real) reps(int) out(string) est(string) [ref(real -1) seed0(int 20260000)]
  file open fh using "`out'", write text replace
  file write fh "rep,rc,k,keff,dfr,slope,se,W,p,mdt50,mdt80,nis,crit,b1,b2,b3,s11,s22,s33,s12,s13,s23" _n
  file close fh
  forvalues r = 1/`reps' {
    set seed `= `seed0' + `r''
    drawy `delta'
    cap qui wooldid $WDCALL , `est'
    local rc = _rc
    local line "`r',`rc'"
    if `rc' == 0 {
      ptpstat , ref(`ref')
      if r(ok) == 1 {
        tempname bp S
        matrix `bp' = r(bp)
        matrix `S'  = r(S)
        local line "`line',`r(k)',`r(keff)',`r(dfr)',`r(slope)',`r(se)',`r(W)',`r(p)',`r(mdt1)',`r(mdt2)',`r(nis)',`r(crit)'"
        local line "`line',`=`bp'[1,1]',`=`bp'[2,1]',`=`bp'[3,1]'"
        local line "`line',`=`S'[1,1]',`=`S'[2,2]',`=`S'[3,3]',`=`S'[1,2]',`=`S'[1,3]',`=`S'[2,3]'"
      }
      else local line "`line',.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,."
    }
    else local line "`line',.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,."
    file open fh using "`out'", write text append
    file write fh "`line'" _n
    file close fh
  }
end
