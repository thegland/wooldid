* p9 -- EXACT IDENTITY along the DATA axis: the same spec on 20 named realistic panel
* shapes, analytic vs fastnose, every e() member dumped in %21x.  The documented
* envelope is abs ~1e-15 / rel ~1e-13 (sthlp:295).  Anything bigger is a defect; a row
* or name set that differs at all is a defect.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fastnose"
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
run "`WD'/wdump.do"
which wooldid
cap erase "`WD'/p9dump.txt"

local shapes canon unbal rotate admin tinytreated rcs biennial holeyear fewG5 pw_survey w_spread3 wzero miss_y sg_gaps sg_tiny clustonset onsetgap hetero_unbal longT manyG

local sn 0
foreach SH of local shapes {
  local ++sn
  local SHAPE `SH'
  include "`GENDIR'/gen_shapes.do"
  dis "########## SHAPE `SH' (`sn') : `SHAPEDESC'   N=" _N
  qui save "`WD'/p9panel.dta", replace
  foreach mode in AN FN {
    local ex
    if "`mode'" == "FN" local ex fastnose
    use "`WD'/p9panel.dta", clear
    ereturn clear
    cap noi qui wooldid y unit t ttre , att att_it att_i att_itime esprelength(3) espostlength(3) `ex'
    dis "SH_`SH'_`mode'_RC=" _rc
    if _rc == 0 wdump , tag("`SH'_`mode'") file("`WD'/p9dump.txt")
  }
  * reset the dispatcher's locals for the next shape
  local SHAPEDESC
}
dis "DONE_P9"
