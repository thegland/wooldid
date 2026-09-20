clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local OUT    "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* ---- ex_canon : zero-noise balanced staggered, effects heterogeneous WITHIN cohort
local gseed 811001
local gysd 0
local gtauhet 0.5
local gN 60
local gT 8
local gonsets 2003 2005
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_canon.dta", replace
dis "EXBUILT canon " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* ---- ex_unbal : zero noise, late entry / early exit / holes
local gseed 811002
local gysd 0
local gtauhet 0.5
local gN 60
local gT 9
local gonsets 2003 2005
local gentry 0.35
local gexit 0.30
local ggap 0.08
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_unbal.dta", replace
dis "EXBUILT unbal " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* ---- ex_wts : zero noise with unit-level weights spanning 2 orders
local gseed 811003
local gysd 0
local gtauhet 0.5
local gN 60
local gT 8
local gonsets 2003 2005
local gwtype unit
local gwspread 2
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_wts.dta", replace
dis "EXBUILT wts " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* ---- ex_rcs : zero noise, 3 rows per unit-period cell (repeated cross-section)
local gseed 811004
local gysd 0
local gtauhet 0.5
local gN 30
local gT 7
local gpercell 3
local gonsets 2003 2005
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_rcs.dta", replace
dis "EXBUILT rcs " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* ---- ex_onecoh : zero noise, SINGLE treated cohort, balanced -> all 4 estimands equal
local gseed 811005
local gysd 0
local gtauhet 0.5
local gN 60
local gT 8
local gonsets 2004
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_onecoh.dta", replace
dis "EXBUILT onecoh " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* ---- ex_bien : zero noise, biennial grid
local gseed 811006
local gysd 0
local gtauhet 0.5
local gN 60
local gT 8
local gtstep 2
local gonsets 2004 2008
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_bien.dta", replace
dis "EXBUILT bien " _N
