* p05: cleanmatrices and the Mata kernels -- does `clear mata' inside the worker
* change THIS run's e(V), and does it poison the SESSION for later runs?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid
dis "ADOSIZE=" c(adosize)

cap prog drop wdvst
prog def wdvst
  args tag
  dis "VST|`tag'|vcvstate=" "`e(wdid_vcvstate)'"
  dis "VST|`tag'|vblocks=" "`e(wdid_vblocks)'"
  dis "VST|`tag'|vblockn=" e(wdid_vblockn)
  capture confirm matrix e(V)
  if _rc == 0 {
    tempname VV
    matrix `VV' = e(V)
    * count strictly-nonzero off-diagonal entries
    local n = rowsof(`VV')
    local off = 0
    forvalues a = 1/`n' {
      forvalues b = 1/`n' {
        if `a' != `b' {
          if `VV'[`a',`b'] != 0 & !missing(`VV'[`a',`b']) local ++off
        }
      }
    }
    dis "VST|`tag'|Vdim=" `n' "|nonzero_offdiag=" `off'
  }
  dis "VST|`tag'|--MATA-KERNELS-BEGIN--"
  capture noisily mata: mata describe wooldid_*
  dis "VST|`tag'|mata_describe_rc=" _rc
  dis "VST|`tag'|--MATA-KERNELS-END--"
end

* ============ ARM A: reference -- plain run, with the ES window so vsplice has work
wdgen
capture noisily wooldid y unit year ttre, cluster(state) esprelength(3) espostlength(3)
dis "A_RC=" _rc
wdvst A_plain
wddumpe cm_plain

* ============ ARM B: same run WITH cleanmatrices
wdgen
capture noisily wooldid y unit year ttre, cluster(state) esprelength(3) espostlength(3) cleanmatrices
dis "B_RC=" _rc
wdvst B_cleanmatrices
wddumpe cm_clean

* ============ ARM C: a PLAIN run immediately after the cleanmatrices run.
* If the session's Mata kernels did not come back, this run silently loses the
* PS-62 covariance splice too.
wdgen
capture noisily wooldid y unit year ttre, cluster(state) esprelength(3) espostlength(3)
dis "C_RC=" _rc
wdvst C_plain_after_clean
wddumpe cm_after

* ============ ARM D: does a mata-dependent feature still work after cleanmatrices?
wdgen
capture noisily wooldid y unit year ttre, cluster(state) esprelength(3) espostlength(3) gteffects
dis "D_RC=" _rc
capture confirm matrix e(gteffects)
dis "D_GTEFFECTS_PRESENT=" (_rc == 0)

* ============ ARM E: user's own Mata object destroyed by cleanmatrices
wdgen
mata: myprecious = 42
capture mata: mata describe myprecious
dis "E_MYPRECIOUS_BEFORE_RC=" _rc
matrix USERMAT = J(2,2,1)
capture noisily wooldid y unit year ttre, cluster(state) cleanmatrices
dis "E_RC=" _rc
capture mata: mata describe myprecious
dis "E_MYPRECIOUS_AFTER_RC=" _rc
capture confirm matrix USERMAT
dis "E_USERMAT_AFTER_RC=" _rc

* ============ ARM F: adosize raised -- does the ado still get re-read?
set adosize 10000
dis "F_ADOSIZE=" c(adosize)
wdgen
capture noisily wooldid y unit year ttre, cluster(state) esprelength(3) espostlength(3) cleanmatrices
dis "F_RC=" _rc
wdvst F_clean_bigado
wdgen
capture noisily wooldid y unit year ttre, cluster(state) esprelength(3) espostlength(3)
dis "F2_RC=" _rc
wdvst F2_plain_after_clean_bigado
wddumpe cm_after_bigado

dis "P05_DONE"
