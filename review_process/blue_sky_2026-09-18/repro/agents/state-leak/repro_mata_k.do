*! SELF-CONTAINED REPRO -- state-leak agent, finding 3 (UNCERTAIN: attribution open)
*!
*! A single plain wooldid run leaves three globals in the user's MATA workspace --
*! __00000L, __00000M and, notably, `k' -- and two scalars in the user's Stata
*! scalar namespace (__gtools_gfile_ghdfeabs, __gtools_gfile_gregvcov).
*! The `__gtools_' pair is self-attributing.  This file asks the question that
*! matters: does a user's own Mata variable named `k' survive a wooldid run?
*! It also asks who creates it -- reghdfe alone, gtools alone, or wooldid.
*!
*! Run:  cd <this directory>
*!       "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_mata_k.do
*!       then read repro_mata_k.log
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop mkpanel
prog def mkpanel
  qui {
    clear
    set obs 60
    gen int unit  = _n
    gen int state = mod(_n-1, 12) + 1
    gen int ttre = .
    replace ttre = 2005 if mod(_n,4) == 1
    replace ttre = 2007 if mod(_n,4) == 2
    replace ttre = 2009 if mod(_n,4) == 3
    gen double ufe = rnormal()
    expand 12
    bys unit: gen int year = 2000 + _n
    gen byte post = (year >= ttre) & !missing(ttre)
    gen double y = 1 + ufe + 0.03*(year-2000) + 0.5*post + rnormal()
  }
end

* ---------------- 1: the consequence question
mkpanel
mata: k = 12345
mata: myvec = (1, 2, 3)
dis "BEFORE: user's mata k and myvec"
capture noisily mata: k
capture noisily mata: myvec
capture noisily wooldid y unit year ttre, cluster(state) att
dis "WOOLDID_RC=" _rc
dis "AFTER: is the user's mata k still 12345?"
capture noisily mata: k
dis "K_READ_RC=" _rc
capture noisily mata: myvec
dis "MYVEC_READ_RC=" _rc
dis "--- full mata workspace after the run ---"
capture noisily mata: mata describe
dis "--- scalar namespace after the run ---"
capture noisily scalar dir

* ---------------- 2: attribution -- reghdfe alone, in a session with no wooldid
clear all
mkpanel
dis "--- mata workspace in a FRESH session ---"
capture noisily mata: mata describe
qui reghdfe y post, absorb(unit year)
dis "--- after reghdfe ALONE (no wooldid in this session) ---"
capture noisily mata: mata describe
capture noisily scalar dir

* ---------------- 3: attribution -- gtools alone
clear all
mkpanel
capture noisily gegen double gm = mean(y), by(unit)
dis "--- after gtools ALONE ---"
capture noisily mata: mata describe
capture noisily scalar dir

* ==================================================================
* MEASURED (2026-09-17, pinned build c1d8979), from the full dirs snapshot in
* p16_snapdiff.do arm s1 -- `diff dirs_s1_pre.log dirs_s1_post.log':
*   --SCALARDIR--  (empty before)
*     __gtools_gfile_ghdfeabs =         65
*     __gtools_gfile_gregvcov =         65
*   --MATADESC--   (empty before)
*                8   real scalar                 __00000L
*                8   real scalar                 __00000M
*                8   real scalar                 k
* wooldid.ado's ONLY top-level Mata assignments are to `wooldid_Mtmp'
* (ado:6689, 8121, 8791, 8795), which are namespaced and do not fire on a plain
* run, so `k' is very likely a dependency's.  This file settles it.
dis "REPRO_DONE"
