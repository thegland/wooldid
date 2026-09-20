clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local GENDIR "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local OUT    "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"

* --- ex_noai: NO unit FE, zero noise -> i=cohort is exact too
local gseed 822001
local gysd 0
local gaisd 0
local gtauhet 0.5
local gN 60
local gT 8
local gonsets 2003 2005
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_noai.dta", replace
dis "EX2BUILT noai " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* --- ex_noai_unbal: no unit FE, unbalanced
local gseed 822002
local gysd 0
local gaisd 0
local gtauhet 0.5
local gN 60
local gT 9
local gonsets 2003 2005
local gentry 0.35
local gexit 0.30
local ggap 0.08
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_noai_unbal.dta", replace
dis "EX2BUILT noai_unbal " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* --- ex_sub: repeated cross-sections with households FINER than unit (T1/T2/T3 surface)
local gseed 822003
local gysd 0
local gtauhet 0.5
local gN 36
local gT 7
local gpercell 4
local gsub 3
local gonsets 2003 2005
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_sub.dta", replace
dis "EX2BUILT sub " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* --- ex_wzero: zero noise, row weights spanning 3 orders with 10% exact zeros
local gseed 822004
local gysd 0
local gtauhet 0.5
local gN 60
local gT 8
local gonsets 2003 2005
local gwtype row
local gwspread 3
local gwzero 0.10
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_wzero.dta", replace
dis "EX2BUILT wzero " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* --- ex_ymiss: zero noise with 12% missing outcomes (realistic attrition)
local gseed 822005
local gysd 0
local gtauhet 0.5
local gN 60
local gT 8
local gonsets 2003 2005
local gymiss 0.12
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_ymiss.dta", replace
dis "EX2BUILT ymiss " _N

foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
* --- ex_many: 5 cohorts, uneven sizes, longer T (a harder weighting problem)
local gseed 822006
local gysd 0
local gtauhet 0.6
local gN 90
local gT 12
local gshares 1 3 1 5 2
local gonsets 2003 2004 2006 2008 2009
include "`GENDIR'/gen_core.do"
qui save "`OUT'/ds_ex_many.dta", replace
dis "EX2BUILT many " _N
