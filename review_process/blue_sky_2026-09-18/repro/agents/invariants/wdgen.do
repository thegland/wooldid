* ============================================================================
* wdgen.do -- canonical DGPs + snapshot helper for the blue-sky invariant probes.
* Usage:  do wdgen.do        (defines programs; generates nothing)
*         wdgen_d1           -> unit-level panel  (i == unit)
*         wdgen_d2           -> coarse-cohort panel (many units inside each i-t cell)
*         wdsnap using f.txt , tag(base)   -> dump e() in %21x for exact comparison
* Both DGPs are deterministic given `set seed' / `set sortseed' before the call.
* ============================================================================
cap prog drop wdgen_d1
prog def wdgen_d1
  syntax , [nu(integer 60) nt(integer 10)]
  clear
  qui set obs `nu'
  gen long i = _n
  gen byte arm = mod(_n,4)
  gen double ttre = .
  qui replace ttre = 2003 if arm==1
  qui replace ttre = 2005 if arm==2
  qui replace ttre = 2007 if arm==3
  gen double ai = rnormal()
  gen long clus = ceil(i/3)            /* 3 units per cluster : COARSER than i */
  gen byte sg = mod(_n,2)              /* subgroup, constant within i          */
  qui expand `nt'
  bysort i: gen int t = 1999 + _n
  gen double x  = rnormal()            /* within-cell-varying control          */
  gen double w  = 1 + runiform()       /* sample weight                        */
  gen double rel = t - ttre
  gen double te  = cond(ttre<. & t>=ttre, 1 + 0.35*rel + 0.5*arm, 0)
  gen double y   = 3 + ai + 0.20*(t-1999) + 0.7*x + te + rnormal()
  order i t ttre y x w sg clus arm
end

cap prog drop wdgen_d2
prog def wdgen_d2
  syntax , [nc(integer 8) nj(integer 12) nt(integer 10)]
  clear
  qui set obs `nc'
  gen long i = _n                       /* i is a COARSE cohort                */
  gen double ttre = .
  qui replace ttre = 2003 if inlist(i,1,2)
  qui replace ttre = 2005 if inlist(i,3,4)
  qui replace ttre = 2007 if inlist(i,5,6)
  gen double ai = rnormal()
  qui expand `nj'
  bysort i: gen long j = _n             /* unit inside the cohort              */
  gen long uid = _n                     /* globally unique unit id             */
  gen double aj = rnormal()
  gen byte sg = mod(i,2)                /* constant within i (sgdom requires)  */
  qui expand `nt'
  bysort uid: gen int t = 1999 + _n
  gen double x   = rnormal()
  gen double dose = runiform()*2        /* continuous treatment                */
  gen double w   = 1 + runiform()
  gen double rel = t - ttre
  gen double te  = cond(ttre<. & t>=ttre, 1 + 0.30*rel + 0.8*dose, 0)
  gen double y   = 2 + ai + aj + 0.15*(t-1999) + 0.5*x + te + rnormal()
  order i uid j t ttre y x dose w sg
end

* ---------------------------------------------------------------------------
* wdsnap: append an exact (%21x) dump of the interesting parts of e() to a file.
* Everything a metamorphic comparison needs, in one greppable block per tag.
* ---------------------------------------------------------------------------
cap prog drop wdsnap
prog def wdsnap
  syntax using/ , tag(string) [mats(string)]
  qui file open wdsf using "`using'", write append
  file write wdsf "=== TAG `tag'" _n
  foreach s in N N_clust df_r {
    cap local v = e(`s')
    if _rc==0 file write wdsf "scal `s' = " %21x (`v') _n
  }
  tempname B V
  cap matrix `B' = e(b)
  if _rc==0 {
    local cn : colnames `B'
    local ce : coleq `B'
    forval k=1/`=colsof(`B')' {
      local nm : word `k' of `cn'
      local eq : word `k' of `ce'
      file write wdsf "b[`eq':`nm'] = " %21x (`B'[1,`k'])
      cap matrix `V' = e(V)
      if _rc==0 file write wdsf "  v = " %21x (`V'[`k',`k'])
      file write wdsf _n
    }
  }
  foreach m of local mats {
    cap matrix `B' = e(`m')
    if _rc==0 {
      local rn : rownames `B'
      local cn : colnames `B'
      forval r=1/`=rowsof(`B')' {
      forval k=1/`=colsof(`B')' {
        local rr : word `r' of `rn'
        local kk : word `k' of `cn'
        file write wdsf "`m'[`rr',`kk'] = " %21x (`B'[`r',`k']) _n
      }
      }
    }
    else file write wdsf "`m' ABSENT" _n
  }
  file close wdsf
end
