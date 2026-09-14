* se_cov_mc.do - the 500-draw standard-error coverage MC that Part 8 C.6 names as the
* measurement nobody had run. 2026-09-14.
*
* WHAT C.6 LEFT OPEN. Thirty draws at G = 40 showed wooldid's analytic SEs sitting 20-40% below
* the sampling SD on two arms (staggered esfixedbaseperiod with i = the unit id; and the pooled
* esprelength(0) reference in BOTH DGPs), with csdid and did_imputation - which return the SAME
* points to 1e-14 / 1e-7 on the same data with the same cluster(i) - covering correctly in the
* same draws. Because the points agree, the benchmarks isolate the variance path exactly. Thirty
* draws is a lead, not a finding. C.6 specifies the shape needed before anyone rules: 500+ draws,
* arms {esfixedbaseperiod, pooled} x {i = unit, i = cohort} x {common onset, staggered}, G swept,
* both external packages as controls in every cell. This driver is that shape, sharded.
*
* The author's 2026-09-13 ruling removed the help sentence that steered users toward the geometry
* in question ("call it a day; I am skeptical of it"). That closed the HELP half. It did not
* measure the SEs, and a removed sentence is not evidence that a variance formula is right.
*
* DGP: se_mc.do's, scaled in G. G units, T = 12 periods.
*   hom  - common onset t = 7 for the first half of units, tau = 0.5
*   stag - onset 5 for the first quarter, onset 8 for the second quarter, rest never treated;
*          tau = 1.0 (early) or 0.4 (late), plus 0.05 per period since onset
* Cohort SHARES are held fixed as G varies, so the analytic truth in semc_analyze.py is unchanged.
*
* Arms per draw: wooldid esfixedbaseperiod with i = unit (esf/esfeb) and with i = the cohort id
* (esfc/esfceb); wooldid pooled esprelength(0) (pool/pooleb); csdid long2 + estat event (cs);
* did_imputation (di). The two external arms are CONTROLS: a wooldid arm below band while both
* controls are inside band in the SAME draws is the finding; all three below band is the DGP.
*
* args: 1=ado 2=shard 3=ndraws 4=seedbase 5=outdir 6=G 7=kinds (+ separated; default hom+stag)
* Writes <outdir>/se_cov_s<shard>_G<G>.txt in se_mc.do's pipe format, so semc_analyze.py reads it
* with only a filename change. Every row carries its shard, draw and G.
clear all
set more off
set maxvar 20000
set graphics off

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local SHARD = `2'
local NDRAWS = `3'
local BASE  = `4'
local OUTDIR `5'
local G = `6'
local KINDS `7'
if "`KINDS'" == "" local KINDS hom+stag
local KINDS = subinstr(subinstr("`KINDS'", "+", " ", .), ",", " ", .)

file open MC using "`OUTDIR'/se_cov_s`SHARD'_G`G'.txt", write replace

* se_mc.do's mkd, with G a parameter. Cohort shares are FIXED so the truth is G-invariant.
cap program drop mkd
program define mkd
  args kind G
  clear
  qui set obs `G'
  gen double i = _n
  gen double a = rnormal()
  local Q1 = ceil(`G'/4)
  local H  = ceil(`G'/2)
  if "`kind'" == "stag" gen double ttre = cond(i <= `Q1', 5, cond(i <= `H', 8, .))
  else gen double ttre = cond(i <= `H', 7, .)
  qui expand 12
  bysort i: gen double t = _n
  gen double rel = t - ttre
  gen double tau = 0
  if "`kind'" == "hom" qui replace tau = 0.5 if ttre < . & t >= ttre
  if "`kind'" == "stag" qui replace tau = cond(ttre == 5, 1.0, 0.4) + 0.05*rel if ttre < . & t >= ttre
  gen double y = a + 0.2*t + 0.5*cos(t) + tau + rnormal()*0.4
  gen double g = cond(ttre < ., ttre, 0)
  gen double icoh = cond(ttre < ., ttre, i)
end

cap program drop wv
program define wv
  args tag r B V
  local names : colnames `B'
  local k = colsof(`B')
  forval j = 1/`k' {
    local nm : word `j' of `names'
    file write MC "`tag'|`r'|`nm'|" %21.0g (`B'[1,`j']) "|" %21.0g (sqrt(`V'[`j',`j'])) _n
  }
end

timer clear 1
timer on 1
forval d = 1/`NDRAWS' {
  local r = (`SHARD' - 1)*100000 + `d'
  foreach kind of local KINDS {
    * one seed per (shard, draw, kind): shards never share a stream
    local sd = `BASE' + `SHARD'*1000000 + `d'*10 + cond("`kind'" == "hom", 1, 2)
    set seed `sd'
    set sortseed `sd'
    mkd `kind' `G'
    tempfile dd
    qui save `dd'
    local P 5
    local Q 6
    if "`kind'" == "stag" local P 6
    if "`kind'" == "stag" local Q 8
    cap wooldid y i t ttre, att cluster(i) esprelength(`P') espostlength(`Q') esfixedbaseperiod
    if _rc == 0 {
      matrix B = e(b_es_att)
      matrix V = e(V_es_att)
      wv esf_`kind' `r' B V
      matrix B = e(b)
      matrix V = e(V)
      wv esfeb_`kind' `r' B V
    }
    else file write MC "err|esf_`kind'|`r'|" (_rc) _n
    use `dd', clear
    cap wooldid y icoh t ttre, att cluster(i) esprelength(`P') espostlength(`Q') esfixedbaseperiod
    if _rc == 0 {
      matrix B = e(b_es_att)
      matrix V = e(V_es_att)
      wv esfc_`kind' `r' B V
      matrix B = e(b)
      matrix V = e(V)
      wv esfceb_`kind' `r' B V
    }
    else file write MC "err|esfc_`kind'|`r'|" (_rc) _n
    use `dd', clear
    cap csdid y, ivar(i) time(t) gvar(g) long2
    if _rc == 0 {
      cap estat event, post
      if _rc == 0 {
        matrix B = e(b)
        matrix V = e(V)
        wv cs_`kind' `r' B V
      }
      else file write MC "err|csev_`kind'|`r'|" (_rc) _n
    }
    else file write MC "err|cs_`kind'|`r'|" (_rc) _n
    use `dd', clear
    cap wooldid y i t ttre, att cluster(i) esprelength(0) espostlength(`Q')
    if _rc == 0 {
      matrix B = e(b_es_att)
      matrix V = e(V_es_att)
      wv pool_`kind' `r' B V
      matrix B = e(b)
      matrix V = e(V)
      wv pooleb_`kind' `r' B V
    }
    else file write MC "err|pool_`kind'|`r'|" (_rc) _n
    use `dd', clear
    cap wooldid y icoh t ttre, att cluster(i) esprelength(0) espostlength(`Q')
    if _rc == 0 {
      matrix B = e(b)
      matrix V = e(V)
      wv poolceb_`kind' `r' B V
    }
    else file write MC "err|poolc_`kind'|`r'|" (_rc) _n
    use `dd', clear
    cap did_imputation y i t ttre, horizons(0/`=`Q'-1') minn(0) cluster(i)
    if _rc == 0 {
      matrix B = e(b)
      matrix V = e(V)
      wv di_`kind' `r' B V
    }
    else file write MC "err|di_`kind'|`r'|" (_rc) _n
  }
  timer off 1
  qui timer list 1
  file write MC "progress|`r'|" %9.1f (r(t1)) _n
  di "SECOV_PROGRESS shard=`SHARD' G=`G' draw=`d' cum=" %8.1f r(t1)
  timer on 1
}
file write MC "SECOV_DONE `SHARD' `G'" _n
file close MC
di "SECOV_SHARD_DONE `SHARD' `G'"
exit, clear
