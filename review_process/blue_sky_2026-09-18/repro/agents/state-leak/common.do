*! common.do -- state-differencing harness for wooldid blue-sky probing
*! defines: wdsnap <tag>   (writes snap_<tag>.txt and dirs_<tag>.log)
*!          wdgen          (builds a realistic staggered panel in memory)
*!          wddumpe <tag>  (dumps e() in %21x)

cap prog drop wdsnap
prog def wdsnap
  args tag
  qui {
    tempname fh
    file open `fh' using "snap_`tag'.txt", write text replace

    file write `fh' "=== C-CLASS ===" _n
    foreach cc in N k changed filename frame level emptycells type varabbrev more rmsg trace matsize maxvar linesize sortlist seed dp scheme adosize maxiter matalnum matastrict fvtrack fvwrap noisily processors {
      local vv
      capture local vv = c(`cc')
      if _rc != 0 {
        capture local vv "`c(`cc')'"
      }
      file write `fh' "[C] `cc' = `vv'" _n
    }
    * rngstate is thousands of characters -- never `dis' it (batch logs wrap at 80)
    local rs "`c(rngstate)'"
    file write `fh' "[C] rngstate_len = " %12.0f (strlen("`rs'")) _n
    file write `fh' "[C] rngstate = `rs'" _n
    * sort rng state, if this Stata exposes one
    local srs
    capture local srs "`c(sortrngstate)'"
    file write `fh' "[C] sortrngstate = `srs'" _n
    local sq
    capture local sq = c(rngseed_mt64s)
    file write `fh' "[C] rngseed_mt64s = `sq'" _n
    file write `fh' "[C] rngstream = `c(rngstream)'" _n
    file write `fh' "[C] rng = `c(rng)'" _n

    file write `fh' "=== DATA ===" _n
    capture confirm variable _n_never_exists
    local nvars = c(k)
    file write `fh' "[D] nobs = " %12.0f (_N) _n
    file write `fh' "[D] nvars = `nvars'" _n
    if `nvars' > 0 {
      foreach v of varlist _all {
        local ty : type `v'
        local fmt : format `v'
        local vl : value label `v'
        local lb
        capture mata: st_local("lb", subinstr(subinstr(subinstr(st_varlabel(st_varindex(st_local("v"))), char(34), "~"), char(96), "^"), char(36), "#"))
        file write `fh' "[VAR] `v' | `ty' | `fmt' | vl=`vl' | lbl=`lb'" _n
      }
      capture datasignature
      file write `fh' "[D] datasignature = `r(datasignature)'" _n
      * row-order fingerprint: __wdrowid is planted by wdgen and must equal _n
      capture confirm variable __wdrowid
      if _rc == 0 {
        capture assert __wdrowid == _n
        file write `fh' "[D] roworder_intact = " %4.0f (_rc == 0) _n
        capture sum __wdrowid, meanonly
      }
    }
    * dataset characteristics
    local chs : char _dta[]
    file write `fh' "[D] dtachars = `chs'" _n
    foreach ch of local chs {
      local cv : char _dta[`ch']
      file write `fh' "[DCHAR] `ch' = " %12.0f (strlen("`macval(cv)'")) _n
    }

    file write `fh' "=== GLOBALS ===" _n
    local gl : all globals
    foreach g of local gl {
      * skip harness-owned names
      if substr("`g'",1,4) == "WDH_" continue
      * several of Stata's own S_* globals hold a COMPOUND-QUOTED value, whose
      * quotes/backticks/dollars break any Stata line the value is substituted into.
      * Fetch AND sanitise inside Mata (char() literals only, so nothing the macro
      * expander can see), then write the sanitised text.
      local gv
      capture mata: st_local("gv", subinstr(subinstr(subinstr(st_global(st_local("g")), char(34), "~"), char(96), "^"), char(36), "#"))
      if _rc != 0 local gv <UNREADABLE>
      file write `fh' "[GLOBAL] `g' = `gv'" _n
    }

    file write `fh' "=== FILES-IN-PWD (see shell ls) ===" _n
    file close `fh'

    * things with no programmatic accessor go to a named log
    capture log close wdsnaplog
    log using "dirs_`tag'.log", text replace name(wdsnaplog)
    noi dis "--MATRIXDIR--"
    capture noisily matrix dir
    noi dis "--SCALARDIR--"
    capture noisily scalar dir
    noi dis "--ESTIMATESDIR--"
    capture noisily estimates dir
    noi dis "--FRAMESDIR--"
    capture noisily frames dir
    noi dis "--MATADESC--"
    capture noisily mata: mata describe
    noi dis "--LABELDIR--"
    capture noisily label dir
    noi dis "--ENDDIRS--"
    log close wdsnaplog
  }
end

* ---------------------------------------------------------------------------
* wdgen: a small, messy-but-plausible staggered panel.
*   60 units in 12 states, t = 2001..2012, three treated cohorts + never-treated
*   carries: variable labels, value labels, non-default formats, a subgroup var,
*            an analytic weight, and a planted row-order fingerprint __wdrowid
cap prog drop wdgen
prog def wdgen
  syntax [, nunits(integer 60) nt(integer 12) seed(integer 20260916)]
  qui {
    clear
    set seed `seed'
    set obs `nunits'
    gen int unit = _n
    gen int state = mod(_n-1, 12) + 1
    gen byte sg = mod(_n-1, 2)
    * onset: 1/4 never treated, rest in 2005 / 2007 / 2009
    gen int ttre = .
    replace ttre = 2005 if mod(_n,4) == 1
    replace ttre = 2007 if mod(_n,4) == 2
    replace ttre = 2009 if mod(_n,4) == 3
    gen double ufe = rnormal()
    expand `nt'
    bys unit: gen int year = 2000 + _n
    gen double tfe = 0.03*(year-2000)
    gen byte post = (year >= ttre) & !missing(ttre)
    gen double eps = rnormal()
    * NB: 0*. is MISSING in Stata, so the dynamic term must be guarded or every
    * never-treated row's outcome goes missing and wooldid (correctly) refuses r459
    gen double y = 1 + ufe + tfe + 0.5*post + 0.2*post*cond(missing(ttre), 0, year-ttre) + eps
    gen double wt = 1 + 0.5*runiform()
    * labels and formats a real dataset would carry
    label var y "log earnings"
    label var unit "county id"
    label var year "calendar year"
    label var ttre "policy onset year"
    label var sg "urban indicator"
    label define sglbl 0 "rural" 1 "urban"
    label values sg sglbl
    format y %9.3f
    format wt %6.4f
    format ttre %ty
    sort unit year
    gen long __wdrowid = _n
    compress
    qui count if missing(y)
    assert r(N) == 0
    qui count if missing(ttre)
    assert r(N) > 0
  }
end

* ---------------------------------------------------------------------------
* wddumpe: dump the parts of e() that carry numbers, in exact hex
cap prog drop wddumpe
prog def wddumpe
  args tag
  qui {
    tempname fh
    file open `fh' using "ereturn_`tag'.txt", write text replace
    file write `fh' "cmd=`e(cmd)' N=`e(N)' rc-ok" _n
    foreach m in b V wdidmainresults b_es_att V_es_att gteffects {
      capture confirm matrix e(`m')
      if _rc == 0 {
        tempname MM
        matrix `MM' = e(`m')
        local nr = rowsof(`MM')
        local nc = colsof(`MM')
        file write `fh' "[MAT] `m' `nr' x `nc'" _n
        forvalues r = 1/`nr' {
          local rn : word `r' of `: rownames `MM''
          forvalues cc = 1/`nc' {
            local cn : word `cc' of `: colnames `MM''
            local vv : di %21x `MM'[`r',`cc']
            file write `fh' "  `m'[`rn',`cn'] = `vv'" _n
          }
        }
      }
      else file write `fh' "[MAT] `m' ABSENT" _n
    }
    foreach s in N level df_r {
      capture local sv : di %21x e(`s')
      file write `fh' "[SCA] `s' = `sv'" _n
    }
    file close `fh'
  }
end

* ---------------------------------------------------------------------------
* wdline: one compact, greppable line per state check (batch logs wrap at 80,
* so each fact gets its own short line)
cap prog drop wdline
prog def wdline
  args tag
  dis "WDL|`tag'|level=" c(level)
  dis "WDL|`tag'|emptycells=" "`c(emptycells)'"
  dis "WDL|`tag'|sortlist=" "`c(sortlist)'"
  dis "WDL|`tag'|N=" c(N) "|k=" c(k)
  dis "WDL|`tag'|frame=" "`c(frame)'"
  dis "WDL|`tag'|changed=" c(changed)
  local chs : char _dta[]
  dis "WDL|`tag'|dtachars=" "`chs'"
  local gl : all globals
  local any 0
  foreach g of local gl {
    if substr("`g'",1,8) == "WOOLDID_" {
      local gv
      capture mata: st_local("gv", subinstr(subinstr(subinstr(st_global(st_local("g")), char(34), "~"), char(96), "^"), char(36), "#"))
      dis "WDL|`tag'|GLOBAL|`g'=" "`gv'"
      local any 1
    }
  }
  if `any' == 0 dis "WDL|`tag'|GLOBAL|none"
  capture confirm variable __wdrowid
  if _rc == 0 {
    capture assert __wdrowid == _n
    dis "WDL|`tag'|roworder_intact=" (_rc == 0)
  }
  * rngstate fingerprint: first/last 24 chars + length (never dis the whole thing)
  local rs "`c(rngstate)'"
  dis "WDL|`tag'|rnglen=" strlen("`rs'")
  dis "WDL|`tag'|rnghead=" substr("`rs'",1,24)
  dis "WDL|`tag'|rngtail=" substr("`rs'",-24,24)
end

* wdstore/wdcmp: stash c(rngstate) and compare later
cap prog drop wdrngeq
prog def wdrngeq
  args tag a b
  dis "WDL|`tag'|rng_equal=" ("${`a'}" == "${`b'}")
end
