* cmplib.do -- defines program cmpf: compare two bootsave/permsave .dta draw files
cap prog drop cmpf
prog def cmpf
  syntax , f1(string) f2(string) tag(string) [maxrep(integer 0)]
  qui use "`f1'.dta", clear
  if `maxrep' > 0 qui keep if _rep <= `maxrep'
  qui ds
  local v1 `r(varlist)'
  foreach v of local v1 {
    if "`v'" == "_rep" continue
    qui rename `v' q1`v'
  }
  tempfile t1
  qui save `t1'
  qui use "`f2'.dta", clear
  if `maxrep' > 0 qui keep if _rep <= `maxrep'
  qui ds
  local v2 `r(varlist)'
  qui merge 1:1 _rep using `t1', nogen
  local worst = 0
  local worstv ""
  local ncom = 0
  local nw = 0
  local wworst = 0
  local wworstv ""
  foreach v of local v2 {
    if "`v'" == "_rep" continue
    cap confirm variable q1`v'
    if _rc continue
    local ++ncom
    qui gen double __d = reldif(`v', q1`v')
    qui replace __d = 0 if missing(`v') & missing(q1`v')
    qui sum __d, meanonly
    if r(max) > `worst' {
      local worst = r(max)
      local worstv "`v'"
    }
    if substr("`v'",1,2) == "w_" {
      local ++nw
      qui gen double __e = abs(`v' - q1`v')
      qui sum __e, meanonly
      if r(max) > `wworst' {
        local wworst = r(max)
        local wworstv "`v'"
      }
      qui drop __e
    }
    qui drop __d
  }
  dis "CMP|`tag'|ncommon=" `ncom' "|MAXRELDIF=" %12.4e `worst' "|on=" "`worstv'"
  dis "CMP|`tag'|nwdraw=" `nw' "|WMAXABSDIFF=" %12.4e `wworst' "|on=" "`wworstv'"
  cap confirm variable main_att
  if _rc == 0 {
    cap confirm variable q1main_att
    if _rc == 0 {
      qui gen double __d = reldif(main_att, q1main_att)
      qui sum __d, meanonly
      dis "CMP|`tag'|main_att draws MAXRELDIF=" %12.4e r(max)
      list _rep main_att q1main_att in 1/8, noobs
      qui drop __d
    }
  }
  cap confirm variable main_att_i
  if _rc == 0 {
    cap confirm variable q1main_att_i
    if _rc == 0 {
      qui gen double __d = reldif(main_att_i, q1main_att_i)
      qui sum __d, meanonly
      dis "CMP|`tag'|main_att_i draws MAXRELDIF=" %12.4e r(max)
      qui drop __d
    }
  }
end
