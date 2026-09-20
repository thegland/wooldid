*! harvest.do -- copy the permutation surfaces out of e() into fixed names, then
*! run the independent recheck.  Locals expected: htag, habsp, hexact (0/1).
foreach _m in _PD _CM _OB _PK _PDS _CMS {
  cap matrix drop `_m'
}
matrix _PD = e(permdist_b)
matrix _CM = e(perm_comparator_b)
matrix _OB = e(b)
matrix _PK = e(perm_kernel)
cap confirm matrix e(permdist_se)
if _rc == 0 matrix _PDS = e(permdist_se)
cap confirm matrix e(perm_comparator_se)
if _rc == 0 matrix _CMS = e(perm_comparator_se)
local hlvl = e(perm_level)
* ---- stripe surfaces must equal the kernel rows they are documented to be ----
tempname A B
matrix `A' = e(perm_p)
matrix `B' = _PK[3,1...]
mata: st_local("d", strofreal(max(abs(st_matrix("`A'") - st_matrix("`B'"))))) 
dis "`htag' CHK perm_p==kernel3 maxabsdiff=`d'"
matrix `A' = e(perm_ci)
matrix `B' = _PK[7..8,1...]
mata: st_local("d", strofreal(max(abs(st_matrix("`A'") - st_matrix("`B'")))))
dis "`htag' CHK perm_ci==kernel78 maxabsdiff=`d'"
* ---- estout handoff: e(p)/e(ci_l)/e(ci_u) must be the permutation numbers ----
cap confirm matrix e(p)
if _rc == 0 {
  matrix `A' = e(p)
  matrix `B' = e(perm_p)
  mata: st_local("d", strofreal(max(abs(st_matrix("`A'") - st_matrix("`B'")))))
  dis "`htag' CHK e(p)==e(perm_p) maxabsdiff=`d'"
  matrix `A' = e(ci_l) \ e(ci_u)
  matrix `B' = e(perm_ci)
  mata: st_local("d", strofreal(max(abs(st_matrix("`A'") - st_matrix("`B'")))))
  dis "`htag' CHK e(ci_l/u)==e(perm_ci) maxabsdiff=`d'"
}
else dis "`htag' CHK e(p) ABSENT"
* ---- the displayed table must carry the same p and interval ----
matrix `A' = e(wdidmainresults)
local nr = rowsof(`A')
local worst = 0
forvalues r = 1/`nr' {
  if `A'[`r',1] < . {
    local pk = _PK[3,`r']
    local pt = `A'[`r',4]
    if abs(`pk' - `pt') > `worst' local worst = abs(`pk'-`pt')
    local lk = _PK[7,`r']
    local lt = `A'[`r',5]
    if abs(`lk' - `lt') > `worst' & `lk' < . & `lt' < . local worst = abs(`lk'-`lt')
    local uk = _PK[8,`r']
    local ut = `A'[`r',6]
    if abs(`uk' - `ut') > `worst' & `uk' < . & `ut' < . local worst = abs(`uk'-`ut')
  }
}
dis "`htag' CHK displayedtable==kernel maxabsdiff=`worst'"
mata: recheck("`htag'", `hlvl', `habsp', `hexact')
