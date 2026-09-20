*==============================================================================
* dumpe.do -- defines prog dumpe: writes every visible e() result to a text file
*             in exact hex (%21x) so two runs can be diffed mechanically.
*==============================================================================
cap prog drop dumpe
prog def dumpe
  syntax , f(string)
  tempname fh MM
  file open `fh' using "`f'", write replace text
  file write `fh' "== SCALARS ==" _n
  local SLIST : e(scalars)
  local MLIST : e(macros)
  local TLIST : e(matrices)
  foreach s of local SLIST {
    local v = e(`s')
    local hx = strofreal(`v',"%21x")
    file write `fh' "S `s' = `hx'" _n
  }
  file write `fh' "== MACROS ==" _n
  foreach l of local MLIST {
    local v `"`e(`l')'"'
    file write `fh' `"M `l' = `v'"' _n
  }
  file write `fh' "== MATRICES ==" _n
  foreach m of local TLIST {
    cap mat drop `MM'
    cap mat `MM' = e(`m')
    if _rc != 0 {
      file write `fh' "X `m' = UNREADABLE" _n
      continue
    }
    local nr = rowsof(`MM')
    local nc = colsof(`MM')
    local rn : rownames `MM'
    local cn : colnames `MM'
    file write `fh' "T `m' `nr' x `nc'" _n
    file write `fh' "  rn `rn'" _n
    file write `fh' "  cn `cn'" _n
    forvalues r = 1/`nr' {
      local line ""
      forvalues c = 1/`nc' {
        local v = `MM'[`r',`c']
        local hx = strofreal(`v',"%21x")
        local line `"`line' `hx'"'
      }
      file write `fh' `"  r`r'`line'"' _n
    }
  }
  file close `fh'
end
