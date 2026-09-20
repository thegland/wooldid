*--- fuzz assertion library.  include AFTER clear all, BEFORE any wooldid call.
capture program drop wdfz_snap
program define wdfz_snap, rclass
  qui count
  return local nobs = r(N)
  qui ds
  local vl `r(varlist)'
  return local nvar : word count `vl'
  local sb : sortedby
  return local sort "`sb'"
  qui datasignature
  return local sig "`r(datasignature)'"
end

capture program drop wdfz_assert
program define wdfz_assert, rclass
  syntax , ID(string) [ LEVTXT(string) ENG(string) NOBS0(string) NVAR0(string) SORT0(string) SIG0(string) ]
  local V ""
  if "`levtxt'" == "" local levtxt 95

  * ---------- data integrity (runs whatever the rc) ----------
  qui count
  if r(N) != `nobs0' local V "`V' DATA_NOBS(`=r(N)'vs`nobs0')"
  qui ds
  local vl `r(varlist)'
  local nv : word count `vl'
  if `nv' != `nvar0' local V "`V' DATA_NVAR(`nv'vs`nvar0')"
  local sb : sortedby
  if "`sb'" != "`sort0'" local V "`V' DATA_SORT(`sb'/`sort0')"
  qui datasignature
  if "`r(datasignature)'" != "`sig0'" local V "`V' DATA_SIG(`r(datasignature)'/`sig0')"

  * ---------- e() surface ----------
  if "`e(cmd)'" != "wooldid" {
    local V "`V' NO_ECMD(`e(cmd)')"
    return local viol "`V'"
    exit
  }
  capture qui count if e(sample)
  if _rc local V "`V' NO_ESAMPLE"
  else {
    local nsamp = r(N)
    if `nsamp' <= 0 local V "`V' SAMP_ZERO"
    return local nsamp = `nsamp'
  }
  if missing(e(N)) local V "`V' EN_MISS"
  else if e(N) <= 0 local V "`V' EN_ZERO(`=e(N)')"
  return local eN = e(N)
  capture confirm scalar e(N_clust)
  if !_rc {
    if e(N_clust) < 2 local V "`V' NCLUST_LT2(`=e(N_clust)')"
  }

  * ---------- results matrices ----------
  local mats : e(matrices)
  local nres 0
  foreach m of local mats {
    local keep 0
    if substr("`m'",1,15) == "wdidmainresults" local keep 1
    if substr("`m'",1,13) == "wdidsgresults"   local keep 1
    if `keep' == 0 continue
    capture matrix WWW = e(`m')
    if _rc continue
    local nres = `nres' + 1
    local R = rowsof(WWW)
    local C = colsof(WWW)
    if `C' < 6 {
      local V "`V' MATCOLS(`m'=`C')"
      continue
    }
    local cnames : colnames WWW
    local cn5 : word 5 of `cnames'
    if strpos("`cn5'","`levtxt'") == 0 local V "`V' LEVELNAME(`m':`cn5'/req`levtxt')"
    local rnames : rownames WWW
    forvalues r = 1/`R' {
      local rn : word `r' of `rnames'
      local est = WWW[`r',1]
      local se  = WWW[`r',2]
      local tt  = WWW[`r',3]
      local pp  = WWW[`r',4]
      local lb  = WWW[`r',5]
      local ub  = WWW[`r',6]
      if !missing(`se') & `se' < 0                                    local V "`V' SE_NEG(`m':`rn'=`se')"
      if !missing(`se') & `se'==0 & !missing(`est') & `est'!=0        local V "`V' SE_ZERO(`m':`rn')"
      if !missing(`se') & `se'==0 & !missing(`est') & `est'==0        local V "`V' DEGEN00(`m':`rn')"
      if !missing(`pp') & (`pp' < 0 | `pp' > 1)                       local V "`V' P_RANGE(`m':`rn'=`pp')"
      if !missing(`pp') & `pp' == 0                                   local V "`V' P_ZERO(`m':`rn')"
      if !missing(`lb') & !missing(`ub') & `lb' >= `ub'               local V "`V' CI_ORDER(`m':`rn')"
      if !missing(`lb') & !missing(`est') & `est' < `lb' - 1e-8*max(1,abs(`est')) local V "`V' CI_LO(`m':`rn')"
      if !missing(`ub') & !missing(`est') & `est' > `ub' + 1e-8*max(1,abs(`est')) local V "`V' CI_HI(`m':`rn')"
      if !missing(`se') & `se' > 0 & !missing(`est') & missing(`tt')  local V "`V' T_MISS(`m':`rn')"
      if missing(`est') & (!missing(`se') | !missing(`pp'))           local V "`V' EST_MISS(`m':`rn')"
      if !missing(`se') & `se'>0 & !missing(`est') & !missing(`tt') & "`eng'"=="" {
        local timp = `est'/`se'
        if abs(`tt' - `timp') > 1e-6*max(1,abs(`tt'))                 local V "`V' T_INCONS(`m':`rn'=`tt'/`timp')"
      }
      if !missing(`lb') & !missing(`ub') & !missing(`se') & `se'>0 & "`eng'"=="" {
        local hw = (`ub'-`lb')/(2*`se')
        if `hw' < 0.5 | `hw' > 10                                     local V "`V' CIW_ODD(`m':`rn'=`hw')"
      }
    }
  }
  if `nres' == 0 local V "`V' NO_RESULTMAT"

  * ---------- b / V ----------
  capture matrix BBB = e(b)
  local hasb = (_rc==0)
  capture matrix VVV = e(V)
  local hasv = (_rc==0)
  if `hasb' == 0 local V "`V' NO_EB"
  if `hasb' & `hasv' {
    local kb = colsof(BBB)
    local kv = rowsof(VVV)
    local kc = colsof(VVV)
    if `kb' != `kv' | `kv' != `kc' local V "`V' BV_DIM(b`kb'/V`kv'x`kc')"
    else {
      matrix VVT = VVV'
      local md = mreldif(VVV, VVT)
      if missing(`md') local V "`V' V_HASMISS"
      else {
        if `md' > 1e-10 local V "`V' V_ASYM(`md')"
        matrix DGG = vecdiag(VVV)
        local tr = 0
        local negd = 0
        forvalues k = 1/`kv' {
          local d = DGG[1,`k']
          if !missing(`d') {
            local tr = `tr' + abs(`d')
            if `d' < 0 local negd = `negd' + 1
          }
        }
        if `negd' > 0 local V "`V' V_DIAGNEG(`negd')"
        if `tr' == 0  local V "`V' V_ALLZERO(`kv')"
        if `tr' > 0 & `md' <= 1e-10 & `kv' > 1 {
          matrix VSS = (VVV + VVV')/2
          capture matrix symeigen EXX ELL = VSS
          if _rc == 0 {
            local mine = ELL[1,`kv']
            local maxe = ELL[1,1]
            if `mine' < -1e-8*max(1,abs(`maxe')) local V "`V' V_EIGNEG(`mine'/`maxe')"
          }
        }
      }
    }
  }
  return local viol "`V'"
end
