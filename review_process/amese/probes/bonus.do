* bonus.do — the author's coarse()-vs-contreat() example, with diagnostics
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"

cap program drop diagn
program define diagn
  args tag
  dis as result "DIAG|`tag'|rc|" _rc
  tempname V W
  cap matrix `V' = e(V)
  if _rc == 0 {
    mata: V = st_matrix(st_local("V")); printf("DIAG|`tag'|eV|max|%12.6g|allzero|%1.0f\n", max(abs(V)), max(abs(V)) == 0)
  }
  else dis "DIAG|`tag'|eV|absent"
  cap matrix `V' = e(V_es_ame)
  if _rc == 0 {
    mata: V = st_matrix(st_local("V")); ev = symeigenvalues(V); printf("DIAG|`tag'|VESame|rank|%2.0f|of|%2.0f\n", sum(ev :> 1e-10*max(ev)), cols(V))
  }
  else dis "DIAG|`tag'|VESame|absent"
  cap matrix `W' = e(wdidmainresults)
  if _rc == 0 {
    local r1 = rownumb(`W', "main:att")
    local r2 = rownumb(`W', "main:ame")
    local s1 .
    local s2 .
    if `r1' < . local s1 = `W'[`r1', 2]
    if `r2' < . local s2 = `W'[`r2', 2]
    dis "DIAG|`tag'|attSE|`s1'|ameSE|`s2'"
  }
end

clear
set seed 2010
qui set obs 30
gen i = _n
gen ttre = cond(mod(_n,2)==0, 6, .)
qui expand 6
bysort i: gen unit = _n
gen dose = exp(0.8*runiform() + 0.2*mod(i,5))
qui expand 12
bysort i unit: gen t = _n
gen y = 0.1*i + 0.05*t + rnormal(0,1000)
sort i unit t
tempfile d
save `d'

use `d', clear
dis as result "===== 1: contreat robust ====="
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) esf espre(2) robust
diagn ct_rb
use `d', clear
dis as result "===== 2: coarse robust ====="
cap noi wooldid y i t ttre, att coarse(c.dose) esf espre(2) robust
diagn co_rb
use `d', clear
dis as result "===== 3: contreat cluster(unit) ====="
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(unit) esf espre(2)
diagn ct_cu
use `d', clear
dis as result "===== 4: coarse cluster(unit) ====="
cap noi wooldid y i t ttre, att coarse(c.dose) cluster(unit) esf espre(2)
diagn co_cu
use `d', clear
dis as result "===== 5: contreat cluster(i) ====="
cap noi wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) cluster(i) esf espre(2) verbose
diagn ct_ci
use `d', clear
dis as result "===== 6: coarse cluster(i) ====="
cap noi wooldid y i t ttre, att coarse(c.dose) cluster(i) esf espre(2) verbose
diagn co_ci
dis "BONUS_DONE"
exit, clear
