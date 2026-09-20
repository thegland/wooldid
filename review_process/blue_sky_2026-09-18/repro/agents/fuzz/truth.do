capture program drop wdfz_truth
program define wdfz_truth, rclass
  syntax , SMP(string) IVAR(string) [ WVAR(string) TAU(string) ]
  if "`wvar'" == "" local wvar 1
  if "`tau'"  == "" local tau truetau
  cap drop __trsmp __trpost __trw __trwit __trwip __trtag __trtpi
  cap drop __trw_att __trw_it __trw_i __trw_itime
  qui gen byte   __trsmp  = (`smp') == 1
  qui gen byte   __trpost = __trsmp & ttre < . & t >= ttre
  qui gen double __trw    = `wvar'
  qui egen double __trwit = total(cond(__trsmp, __trw, 0)), by(`ivar' t)
  qui egen double __trwip = total(cond(__trpost, __trw, 0)), by(`ivar')
  qui egen byte   __trtag = tag(`ivar' t) if __trpost
  qui egen double __trtpi = total(__trtag == 1), by(`ivar')
  qui gen double __trw_att   = __trw
  qui gen double __trw_it    = __trw / __trwit
  qui gen double __trw_i     = __trw / __trwip
  qui gen double __trw_itime = __trw / (__trwit * __trtpi)
  foreach _e in att it i itime {
    local _nm `_e'
    if "`_e'" == "it"    local _nm att_it
    if "`_e'" == "i"     local _nm att_i
    if "`_e'" == "itime" local _nm att_itime
    qui summarize `tau' [aw = __trw_`_e'] if __trpost == 1, meanonly
    return scalar `_nm' = r(mean)
  }
  qui count if __trpost == 1
  return scalar npost = r(N)
  drop __trsmp __trpost __trw __trwit __trwip __trtag __trtpi
  drop __trw_att __trw_it __trw_i __trw_itime
end
