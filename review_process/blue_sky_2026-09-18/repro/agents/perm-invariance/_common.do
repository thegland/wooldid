*==========================================================================
* _common.do -- shared generator + dumpers for the perm-invariance probes.
* include'd by every probe do-file in this directory.  Self-contained.
*==========================================================================

*--- a realistic staggered panel: 40 counties in 10 states, 10 annual
*    periods, policy adopted AT THE STATE level in 2004 / 2006 / 2008,
*    three never-adopting states.  i = county (unit), cluster = state.
prog def mkpanel
  syntax [, nunit(integer 40) percl(integer 4) nper(integer 10) tzero(integer 2000) sd(integer 90210)]
  clear
  set seed `sd'
  set sortseed `sd'
  qui set obs `nunit'
  qui gen long unit  = _n
  qui gen long state = ceil(unit/`percl')
  qui gen double ai  = rnormal()
  qui gen double ttre = .
  qui su state, meanonly
  local S = r(max)
  * roughly: first 30% of states adopt 2004, next 30% 2006, next 10% 2008,
  * the rest never
  local a1 = ceil(0.3*`S')
  local a2 = ceil(0.6*`S')
  local a3 = ceil(0.7*`S')
  qui replace ttre = 2004 if state <= `a1'
  qui replace ttre = 2006 if state >  `a1' & state <= `a2'
  qui replace ttre = 2008 if state >  `a2' & state <= `a3'
  qui expand `nper'
  qui bysort unit: gen int t = `tzero' - 1 + _n
  qui gen double tfe  = 0.2*(t - `tzero')
  qui gen double eps  = rnormal()*0.5
  qui gen byte   post = (t >= ttre) & !missing(ttre)
  qui gen double tau  = 0.8 + 0.15*(t - ttre)
  qui replace    tau  = 0 if post == 0
  qui gen double y    = ai + tfe + tau*post + eps
  qui gen double w    = 1 + mod(unit, 5)/4
  qui drop tfe eps tau
  qui order unit state t ttre y
  qui compress
end

*--- dump every e() surface that inference is read off, in exact hex
prog def dmp
  args tag
  tempname fh
  file open `fh' using "dump_`tag'.txt", write text replace
  file write `fh' "TAG `tag'" _n
  foreach s in `: e(scalars)' {
    file write `fh' "S `s' " %21x (e(`s')) _n
  }
  foreach m in `: e(macros)' {
    file write `fh' "M `m' `e(`m')'" _n
  }
  foreach mt in `: e(matrices)' {
    local keep 0
    if substr("`mt'",1,4) == "perm" local keep 1
    if "`mt'" == "wdidmainresults"   local keep 1
    if "`mt'" == "b"                 local keep 1
    if strpos("`mt'", "sg")          local keep 1
    if strpos("`mt'", "results")     local keep 1
    if strpos("`mt'", "band")        local keep 1
    if strpos("`mt'", "step")        local keep 1
    if `keep' {
      tempname MM
      matrix `MM' = e(`mt')
      local nr = rowsof(`MM')
      local nc = colsof(`MM')
      local rn : rowfullnames `MM'
      local cn : colfullnames `MM'
      file write `fh' "X `mt' `nr' `nc'" _n
      file write `fh' "XR `mt' `rn'" _n
      file write `fh' "XC `mt' `cn'" _n
      forval r = 1/`nr' {
        forval c = 1/`nc' {
          file write `fh' "  `r' `c' " %21x (`MM'[`r',`c']) _n
        }
      }
    }
  }
  file close `fh'
end

*--- dump the permsave() ASSIGNMENT block, canonicalised: columns ordered
*    by the NUMERIC cluster id parsed out of the onset_<id> name, so the
*    digest is comparable across cluster relabelings
cap mata: mata drop wd_asg()
mata:
void wd_asg(string scalar outf)
{
    string rowvector vn
    string colvector sel
    real colvector ids, p
    real matrix D
    real scalar k, i, j, fh
    string scalar s
    vn  = st_varname(1..st_nvar())
    sel = J(0,1,"")
    ids = J(0,1,.)
    for (k=1; k<=cols(vn); k++) {
        if (substr(vn[k],1,6)=="onset_") {
            sel = sel \ vn[k]
            ids = ids \ strtoreal(substr(vn[k],7,.))
        }
    }
    fh = fopen(outf,"w")
    if (rows(sel)==0) {
        fput(fh,"NO_ASSIGNMENT_BLOCK")
        fclose(fh)
        return
    }
    p   = order(ids,1)
    sel = sel[p]
    fput(fh, "IDS " + invtokens(strofreal(ids[p]',"%20.0g")))
    D = st_data(., invtokens(sel'))
    for (i=1;i<=rows(D);i++) {
        s = ""
        for (j=1;j<=cols(D);j++) {
            s = s + " " + strofreal(D[i,j],"%21x")
        }
        fput(fh, strofreal(i-1) + s)
    }
    fclose(fh)
}
end

prog def asgdmp
  args tag
  preserve
  qui use "ps_`tag'.dta", clear
  qui sort _rep
  mata: wd_asg("asg_`tag'.txt")
  restore
end
