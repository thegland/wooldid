* ---- helper programs (shared by repro_old.do / repro_new.do) ----
cap program drop dumpres
program define dumpres
  args cfg rc dir
  tempname fh
  file open `fh' using "`dir'/`cfg'.csv", write replace
  file write `fh' "rc|`rc'" _n
  qui ds
  file write `fh' `"varlist|`r(varlist)'"' _n
  local escalars : e(scalars)
  foreach s of local escalars {
    local v : di %21x e(`s')
    file write `fh' "scalar|`s'|`v'" _n
  }
  local emacros : e(macros)
  foreach m of local emacros {
    file write `fh' `"macro|`m'|`e(`m')'"' _n
  }
  local emats : e(matrices)
  foreach m of local emats {
    tempname M
    cap matrix `M' = e(`m')
    if _rc != 0 continue
    local R = rowsof(`M')
    local C = colsof(`M')
    local rn : rownames `M'
    local cn : colnames `M'
    file write `fh' "matnames|`m'|`R'x`C'|`rn'|`cn'" _n
    forval r = 1/`R' {
      forval c = 1/`C' {
        local v : di %21x `M'[`r',`c']
        file write `fh' "mat|`m'|`r'|`c'|`v'" _n
      }
    }
    cap matrix drop `M'
  }
  file close `fh'
end
cap program drop pcell
program define pcell
  args tag m r c
  tempname M S
  cap matrix `M' = e(`m')
  if _rc != 0 {
    dis "`tag'=ABSENT(no e(`m'))"
    exit
  }
  if `r' > rowsof(`M') | `c' > colsof(`M') {
    dis "`tag'=ABSENT(out of range)"
    exit
  }
  scalar `S' = `M'[`r',`c']
  dis "`tag'=" %21x `S' " dec=" %14.0g `S'
end
cap program drop pname
program define pname
  args tag name
  tempname B S
  cap matrix `B' = e(b)
  if _rc != 0 {
    dis "`tag'=ABSENT(no e(b))"
    exit
  }
  local names : colnames `B'
  local k 0
  local pos 0
  local n 0
  foreach w of local names {
    local ++k
    if "`w'" == "`name'" {
      local ++n
      if `pos' == 0 local pos `k'
    }
  }
  if `pos' == 0 {
    dis "`tag'=ABSENT count=0 ncols=" colsof(`B')
    exit
  }
  scalar `S' = `B'[1,`pos']
  dis "`tag'=" %21x `S' " dec=" %14.0g `S' " firstpos=`pos' count=`n' ncols=" colsof(`B')
end
cap program drop pshape
program define pshape
  args tag m
  tempname M
  cap matrix `M' = e(`m')
  if _rc != 0 {
    dis "`tag'=ABSENT"
    exit
  }
  local rn : rownames `M'
  dis "`tag'=" rowsof(`M') "x" colsof(`M') " rownames=[`rn']"
end
cap program drop pscal
program define pscal
  args tag s
  local es : e(scalars)
  local em : e(macros)
  if strpos(" `es' ", " `s' ") > 0 {
    tempname S
    scalar `S' = e(`s')
    dis "`tag'=" %21x `S' " dec=" %14.0g `S' " kind=scalar"
    exit
  }
  if strpos(" `em' ", " `s' ") > 0 {
    dis "`tag'=[`e(`s')'] kind=macro"
    exit
  }
  dis "`tag'=ABSENT"
end
cap program drop fresh
program define fresh
  set seed 98765
  set sortseed 98765
  ereturn clear
  cap macro drop WOOLDID_LASTRC WOOLDID_CLEANUP
  cap matrix drop _all
  cap scalar drop _all
  cap estimates drop _all
end
