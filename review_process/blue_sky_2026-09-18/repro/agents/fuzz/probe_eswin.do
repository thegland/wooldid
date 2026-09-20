* Why does the headline att move when the requested ES window changes?
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
foreach ds in canon unbal {
  foreach w in "" "esprelength(1) espostlength(1)" "esprelength(2) espostlength(3)" "esprelength(3) espostlength(3)" "espostlength(3)" {
    use "`D'/ds_`ds'.dta", clear
    ereturn clear
    capture qui wooldid y unit t ttre, att cluster(region) `w' summarystats
    if _rc==0 {
      matrix M = e(wdidmainresults)
      dis "ESWIN|`ds'|`w'|att=" M[1,1] "|se=" M[1,2] "|N=" e(N) "|refN=" e(escomp_refN) "|trtN=" e(escomp_trtN)
    }
    else dis "ESWIN|`ds'|`w'|RC=`=_rc'"
  }
}
dis "ESWIN_DONE"
