*==========================================================================
* gtlib.do -- self-contained helpers for the gteffects aggregation-identity
*             probes.  include this file, then use:
*                mkpanel , n(..) t(..) onsets("..") [knobs]
*                gtcheck , tag("label")
* Nothing here reads the repo test suite; the panel generator is my own.
*==========================================================================

cap prog drop mkpanel
prog def mkpanel
  syntax , [n(int 60) T(int 10) t0(int 2000) ONSets(string) NEVershare(real 0.35) ///
            YSD(real 0.0) GAPprob(real 0) ENTRYprob(real 0) EXITprob(real 0) ///
            PERcell(int 1) HETero(int 1) SEED(int 20260917) tstep(int 1) AISD(real 1) TSSD(real 0.5) GB1(real 0.4) ///
            TAUlist(string) ANTicip(int 0) TRENDsd(real 0)]
  clear
  if "`onsets'" == "" local onsets "2004 2006 2008"
  set seed `seed'
  set sortseed `seed'
  local nons : word count `onsets'
  qui set obs `n'
  gen long unit = _n
  gen double u = runiform()
  gen byte coh = 0
  qui replace coh = 1 + floor((u - `nevershare')/((1-`nevershare')/`nons')) if u >= `nevershare'
  qui replace coh = `nons' if coh > `nons'
  gen double ttre = .
  forval k = 1/`nons' {
    local ok : word `k' of `onsets'
    qui replace ttre = `ok' if coh == `k'
  }
  gen double ai = rnormal()*`aisd'
  gen double trend_i = rnormal()*`trendsd'
  gen double hmult = 1
  if `hetero' == 1 qui replace hmult = 1 + 0.5*rnormal()
  qui expand `t'
  bys unit: gen int k = _n
  gen int year = `t0' + (k-1)*`tstep'
  drop k
  if `entryprob' > 0 {
    bys unit: gen double ue = runiform() if _n == 1
    bys unit: replace ue = ue[1]
    gen int entdelay = 0
    qui replace entdelay = 1 + floor(runiform()*3) if ue < `entryprob'
    bys unit: replace entdelay = entdelay[1]
    qui drop if year < `t0' + entdelay*`tstep'
    drop ue entdelay
  }
  if `exitprob' > 0 {
    bys unit: gen double ux = runiform() if _n == 1
    bys unit: replace ux = ux[1]
    gen int exlen = 0
    qui replace exlen = 1 + floor(runiform()*3) if ux < `exitprob'
    bys unit: replace exlen = exlen[1]
    qui drop if year > `t0' + (`t'-1-exlen)*`tstep'
    drop ux exlen
  }
  if `gapprob' > 0 {
    gen double ug = runiform()
    qui drop if ug < `gapprob'
    drop ug
  }
  if `percell' > 1 {
    qui expand `percell'
    bys unit year: gen byte sub = _n
  }
  tempvar tsh
  bys year: gen double `tsh' = rnormal() if _n == 1
  bys year: replace `tsh' = `tsh'[1]
  gen double tau = .
  if "`taulist'" != "" {
    forval k = 1/`nons' {
      local tk : word `k' of `taulist'
      qui replace tau = `tk' if coh == `k'
    }
  }
  else {
    qui replace tau = 0.5 + 0.35*coh
  }
  gen byte post = (year >= ttre - `anticip') & ttre < .
  gen double relyr_true = year - ttre
  gen double te = 0
  qui replace te = tau*hmult*(1 + 0.15*relyr_true) if post == 1
  gen double x1 = rnormal()
  gen double y = 2 + ai + `tssd'*`tsh' + trend_i*(year-`t0') + te + `gb1'*x1 + rnormal()*`ysd'
  qui drop u
  gen byte cohp = coh + 1
  sort unit year
end

cap prog drop gtcheck
prog def gtcheck
  syntax , [TAG(string) RTOL(real 1e-10) PRE ES]
  if "`tag'" == "" local tag run
  cap confirm matrix e(gteffects)
  if _rc != 0 {
    dis as result "GT|`tag'|NOTABLE|e(gteffects) absent"
    exit
  }
  mata: gtident("`tag'", `rtol')
  if "`pre'" != "" mata: gtidentpre("`tag'", `rtol')
  if "`es'"  != "" mata: gtidentes("`tag'", `rtol')
end

*--------------------------------------------------------------------------
cap prog drop gtdefcheck
prog def gtdefcheck
  syntax , [TAG(string)]
  if "`tag'" == "" local tag run
  cap confirm matrix e(gteffects)
  if _rc != 0 exit
  preserve
  clear
  qui svmat double e(gteffects), names(col)
  qui keep if relyr >= 0
  cap confirm variable w_att_it
  if _rc == 0 {
    qui sum w_att_it
    local mn = r(mean)
    local sd = r(sd)
    dis as result "GT|`tag'|DEF_att_it|cv=" %10.3e cond(`mn'!=0,`sd'/abs(`mn'),.) " min=" %14.10g r(min) " max=" %14.10g r(max)
  }
  cap confirm variable w_att_i
  if _rc == 0 {
    tempvar s
    qui bys i: egen double `s' = total(w_att_i)
    tempvar tg
    qui bys i: gen byte `tg' = _n==1
    qui sum `s' if `tg'
    dis as result "GT|`tag'|DEF_att_i_percoh|min=" %14.10g r(min) " max=" %14.10g r(max) " n=" r(N)
  }
  cap confirm variable w_att_itime
  if _rc == 0 {
    tempvar s2 cvv tg2
    qui bys i: egen double `s2' = total(w_att_itime)
    qui gen double `cvv' = w_att_itime/`s2'
    qui sum `cvv'
    dis as result "GT|`tag'|DEF_att_itime_shareACROSScoh|min=" %14.10g r(min) " max=" %14.10g r(max)
    qui bys i: gen byte `tg2' = _n==1
    qui sum `s2' if `tg2'
    dis as result "GT|`tag'|DEF_att_itime_percoh|min=" %14.10g r(min) " max=" %14.10g r(max) " n=" r(N)
  }
  restore
end

*--------------------------------------------------------------------------
* gtwincheck : the inwindow column must name exactly the horizons the run
*              reports in e(b) (ES_<est>_pre<k> / ES_<est>_post<k>).
*--------------------------------------------------------------------------
cap prog drop gtwincheck
prog def gtwincheck
  syntax , [TAG(string)]
  if "`tag'" == "" local tag run
  cap confirm matrix e(gteffects)
  if _rc != 0 exit
  local bn : colnames e(b)
  local hz
  foreach v of local bn {
    if regexm("`v'", "^ES_att_post([0-9]+)$") local hz `hz' `=regexs(1)'
    if regexm("`v'", "^ES_att_pre([0-9]+)$")  local hz `hz' `=-1*real(regexs(1))'
  }
  local hz : list sort hz
  preserve
  clear
  qui svmat double e(gteffects), names(col)
  qui levelsof relyr if inwindow == 1, local(iw)
  local iw2
  foreach v of local iw {
    local iw2 `iw2' `=`v''
  }
  local iw2 : list sort iw2
  local miss : list hz - iw2
  local extra : list iw2 - hz
  dis as result "WIN|`tag'|ebhoriz=[`hz']|inwindow=[`iw2']|missing=[`miss']|extra=[`extra']"
  restore
end
