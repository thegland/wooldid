*==========================================================================
* gen_shapes.do -- named data shapes for the wooldid blue-sky review.
*
* USAGE:
*     local GENDIR "/private/tmp/.../scratchpad/bluesky/recon/gen"
*     local SHAPE  admin
*     include "`GENDIR'/gen_shapes.do"
*
* The dispatcher CLEARS every gen_core knob first, so shapes never
* contaminate each other.  After the include, `SHAPEDESC' holds a
* one-line description and the data are in memory.
*
* Set `SHAPE' to "list" to print the catalog and stop.
*
* Every shape is seeded (set seed + set sortseed inside gen_core).
*
* OVERRIDING A SHAPE'S KNOBS.  The dispatcher clears every knob first, so
* setting `gysd' before the include does NOT survive.  Pass overrides as
* name/value PAIRS in `GOVR', which are applied after the shape has set
* its own knobs and just before the engine runs:
*
*     local SHAPE unbal
*     local GOVR  "gysd 0 gT 20 gtauhet 0.4"      <-- exact-recovery, longer, heterogeneous
*     include "`GENDIR'/gen_shapes.do"
*
* `gseedovr' is the one-off shorthand for GOVR "gseed <n>" (replicate draws).
*==========================================================================

if "`GENDIR'" == "" {
  di as error "gen_shapes: set local GENDIR to the directory holding gen_core.do"
  error 198
}

if "`SHAPE'" == "list" {
  di as text "REALISTIC: canon canon_exact unbal rotate admin tinytreated rcs rcs_big"
  di as text "           biennial quinq holeyear longT shortT fewG5 fewG10 manyG"
  di as text "           pw_survey w_spread3 wzero miss_y miss_x sg_gaps sg_tiny sg_big"
  di as text "           dose_spike dose_skew dose_lattice dose_cellconst"
  di as text "           clustonset firstperiod afterend bothends pretrend antic"
  di as text "           poisson poissep nevergrp datesm datesq datesd onepre onsetgap"
  di as text "           hetero_unit hetero_unbal"
  di as text "EDGE:      ed_ttrevary ed_ttremix ed_dup ed_singleton ed_zerovar ed_outlier"
  di as text "           ed_badw ed_collin ed_ctrlconst ed_ctrlistreat ed_never1 ed_nopre"
  di as text "           ed_ifflag ed_float ed_tscale ed_onecell ed_tsingle"
  exit
}

* ---- clear every knob so shapes do not contaminate each other ----
foreach _k in gseed gN gcounts gnevern gshares gnever gnevgrp gT gt0 gtstep ///
              gonsets gdropt gentry gentrymax gexit gexitmax ggap gpercell ///
              gsub gysd gtsd gaisd gtaumode gtaus gtauhet gantic gtrendsd gb1 gb2 gdose  ///
              gwtype gwspread gwzero gwneg gwmiss gymiss gxmiss gpois {
  local `_k'
}
local SHAPEDESC

*=========================== REALISTIC ====================================

if "`SHAPE'" == "canon" {
  local SHAPEDESC "reference: balanced staggered panel, 120 units, T=10, 3 cohorts, 40% never"
}
else if "`SHAPE'" == "canon_exact" {
  local gysd 0
  local gtsd 0.5
  local SHAPEDESC "canon with ZERO noise: wooldid must reproduce the design truth exactly"
}
else if "`SHAPE'" == "unbal" {
  local gentry 0.35
  local gexit  0.30
  local ggap   0.08
  local gT     12
  local SHAPEDESC "unbalanced: late entry, early exit and interior holes"
}
else if "`SHAPE'" == "rotate" {
  local gT       12
  local gentry   0.85
  local gentrymax 8
  local gexit    0.85
  local gexitmax 8
  local gN       300
  local SHAPEDESC "rotating panel: most units observed only 2-4 periods"
}
else if "`SHAPE'" == "admin" {
  local gonsets "2004 2005 2006 2007 2008"
  local gcounts "260 3 2 4 2"
  local gnevern 700
  local gT      12
  local gentry  0.12
  local gexit   0.10
  local SHAPEDESC "administrative: 700 never-treated, one 260-unit cohort, four 2-4 unit cohorts"
}
else if "`SHAPE'" == "tinytreated" {
  local gonsets "2005 2007"
  local gcounts "3 3"
  local gnevern 400
  local gT      10
  local SHAPEDESC "huge never-treated group (400), six treated units in two cohorts"
}
else if "`SHAPE'" == "rcs" {
  local gpercell 8
  local gsub     3
  local gN       60
  local gT       10
  local gwtype   row
  local gwspread 2
  local SHAPEDESC "repeated cross-sections: 8 rows per unit-period cell, row-level weights"
}
else if "`SHAPE'" == "rcs_big" {
  local gpercell 20
  local gsub     4
  local gN       40
  local gT       8
  local gwtype   survey
  local gwspread 3
  local SHAPEDESC "heavy repeated cross-sections: 20 rows per cell, survey weights over 3 decades"
}
else if "`SHAPE'" == "biennial" {
  local gtstep  2
  local gt0     2000
  local gT      10
  local gonsets "2005 2009 2013"
  local SHAPEDESC "biennial panel (even years) with odd-year onsets: NO unit is observed at h=0"
}
else if "`SHAPE'" == "quinq" {
  local gtstep  5
  local gt0     1990
  local gT      7
  local gonsets "2001 2007 2012"
  local SHAPEDESC "quinquennial waves with off-wave onsets (census-style)"
}
else if "`SHAPE'" == "holeyear" {
  local gT     12
  local gdropt "2005 2006"
  local SHAPEDESC "two calendar years absent from the entire panel"
}
else if "`SHAPE'" == "longT" {
  local gT      40
  local gt0     1980
  local gonsets "1995 2002 2010"
  local gN      80
  local SHAPEDESC "long panel, T=40"
}
else if "`SHAPE'" == "onepre" {
  local gT      6
  local gt0     2000
  local gonsets "2001 2002"
  local SHAPEDESC "every treated unit has exactly ONE pre-period, and it is its first observed period"
}
else if "`SHAPE'" == "onsetgap" {
  local gT      12
  local gonsets "2004 2006 2008"
  local SHAPEDESC "the onset period itself is missing from the panel for one cohort (no h=0 row)"
}
else if "`SHAPE'" == "shortT" {
  local gT      3
  local gt0     2000
  local gonsets "2002"
  local gN      80
  local SHAPEDESC "short panel, T=3, single cohort treated in the last period"
}
else if "`SHAPE'" == "fewG5" {
  local gonsets "2004 2006"
  local gcounts "2 1"
  local gnevern 2
  local gT      10
  local SHAPEDESC "five units total: few-cluster inference at its limit"
}
else if "`SHAPE'" == "fewG10" {
  local gonsets "2004 2006"
  local gcounts "3 3"
  local gnevern 4
  local gT      10
  local SHAPEDESC "ten units: few-cluster inference"
}
else if "`SHAPE'" == "manyG" {
  local gN  500
  local gT  5
  local gonsets "2002 2003"
  local SHAPEDESC "500 units, T=5: many clusters, short panel"
}
else if "`SHAPE'" == "pw_survey" {
  local gwtype   survey
  local gwspread 3
  local gN       150
  local SHAPEDESC "survey pweights spanning three orders of magnitude"
}
else if "`SHAPE'" == "w_spread3" {
  local gwtype   unit
  local gwspread 3
  local SHAPEDESC "unit-level weights spanning three orders of magnitude"
}
else if "`SHAPE'" == "wzero" {
  local gwtype unit
  local gwzero 0.10
  local SHAPEDESC "10% of rows carry a weight of exactly zero (dropped silently by the markout)"
}
else if "`SHAPE'" == "miss_y" {
  local gymiss 0.15
  local SHAPEDESC "15% of outcomes missing at random"
}
else if "`SHAPE'" == "miss_x" {
  local gxmiss 0.15
  local SHAPEDESC "15% of x1 missing: controls(x1) shrinks the sample, no-controls run does not"
}
else if "`SHAPE'" == "sg_gaps" {
  local SHAPEDESC "subgroup codes 1/3/7/99 with value labels (use sggap)"
}
else if "`SHAPE'" == "sg_tiny" {
  local gN 200
  local SHAPEDESC "subgroup level 2 holds a handful of treated units only (use sgtiny)"
}
else if "`SHAPE'" == "sg_big" {
  local SHAPEDESC "9-digit administrative subgroup codes (use sgbig): e() name overflow surface"
}
else if "`SHAPE'" == "dose_spike" {
  local gdose spike0
  local SHAPEDESC "continuous treatment with a 40% spike at dose zero"
}
else if "`SHAPE'" == "dose_skew" {
  local gdose skew
  local SHAPEDESC "continuous treatment, lognormal with a heavy right tail"
}
else if "`SHAPE'" == "dose_lattice" {
  local gdose lattice
  local SHAPEDESC "continuous treatment on an integer lattice 0..5"
}
else if "`SHAPE'" == "dose_cellconst" {
  local gdose cellconst
  local SHAPEDESC "dose CONSTANT within a cohort-period cell (state-policy design; known G7 pathology)"
}
else if "`SHAPE'" == "clustonset" {
  local gonsets "2005 2005 2005 2009"
  local gcounts "80 3 2 2"
  local gnevern 60
  local gT      12
  local SHAPEDESC "onsets piled on two dates; one enormous cohort and three tiny ones"
}
else if "`SHAPE'" == "firstperiod" {
  local gonsets "2000 2004 2006"
  local gT      10
  local SHAPEDESC "one cohort treated in the FIRST sample period (no pre-period)"
}
else if "`SHAPE'" == "afterend" {
  local gonsets "2004 2006 2012"
  local gT      10
  local SHAPEDESC "one cohort whose onset is after the sample ends (not-yet-treated controls)"
}
else if "`SHAPE'" == "bothends" {
  local gonsets "2000 2004 2012"
  local gT      10
  local SHAPEDESC "one cohort treated at the first period, one after the sample ends"
}
else if "`SHAPE'" == "hetero_unit" {
  local gtauhet 0.40
  local SHAPEDESC "effects heterogeneous WITHIN cohort: att_i at i=unit differs from att_i at i=cohort"
}
else if "`SHAPE'" == "hetero_unbal" {
  local gtauhet 0.40
  local gentry  0.35
  local gexit   0.30
  local gwtype  unit
  local gwspread 2
  local gT      12
  local SHAPEDESC "within-cohort effect heterogeneity + unbalanced + weights: all four estimands differ"
}
else if "`SHAPE'" == "pretrend" {
  local gtrendsd 0.15
  local SHAPEDESC "treated units trend up before onset: parallel trends VIOLATED by construction"
}
else if "`SHAPE'" == "antic" {
  local gantic 2
  local SHAPEDESC "true effects begin 2 periods BEFORE the recorded onset"
}
else if "`SHAPE'" == "poisson" {
  local gpois 1
  local SHAPEDESC "count outcomes ycount / ycount0 for the poisson path"
}
else if "`SHAPE'" == "poissep" {
  local gpois 1
  local SHAPEDESC "count outcomes with an all-zero cohort (separation)"
}
else if "`SHAPE'" == "nevergrp" {
  local gnevgrp 6
  local SHAPEDESC "never-treated split across 6 cohort codes: usable with i=cohort"
}
else if "`SHAPE'" == "datesm" | "`SHAPE'" == "datesq" | "`SHAPE'" == "datesd" {
  local gT      36
  local gt0     0
  local gonsets "12 20 28"
  local SHAPEDESC "real Stata date variables as t and ttre"
}

*============================= EDGE ======================================

else if substr("`SHAPE'", 1, 3) == "ed_" {
  local gN 80
  local gT 10
  if "`SHAPE'" == "ed_never1"  local gnevern 1
  if "`SHAPE'" == "ed_never1"  local gcounts "20 20 20"
  if "`SHAPE'" == "ed_nopre"   local gonsets "2000 2000 2000"
  if "`SHAPE'" == "ed_tsingle" local gN 40
  local SHAPEDESC "edge shape `SHAPE'"
}
else {
  di as error "gen_shapes: unknown SHAPE `SHAPE'  (try SHAPE=list)"
  error 198
}

if "`gseedovr'" != "" local gseed `gseedovr'

* user overrides, applied last
if `"`GOVR'"' != "" {
  local _gnv : word count `GOVR'
  forval _gq = 1(2)`_gnv' {
    local _gk : word `_gq' of `GOVR'
    local _gq2 = `_gq' + 1
    local _gv : word `_gq2' of `GOVR'
    local `_gk' `_gv'
    di as text "  GOVR: `_gk' = `_gv'"
  }
}

include "`GENDIR'/gen_core.do"

*===================== post-generation shaping ===========================

if "`SHAPE'" == "poissep" {
  * the first treated cohort never records a positive count
  qui replace ycount  = 0 if cohort == 1
  qui replace ycount0 = 0 if cohort == 1
}

if "`SHAPE'" == "datesm" {
  * monthly Stata dates: t = months since 2000m1
  qui replace t    = ym(2000, 1) + t
  qui replace ttre = ym(2000, 1) + ttre
  format t ttre %tm
}
if "`SHAPE'" == "datesq" {
  qui replace t    = yq(2000, 1) + t
  qui replace ttre = yq(2000, 1) + ttre
  format t ttre %tq
}
if "`SHAPE'" == "datesd" {
  * DAILY dates carrying MONTHLY data: t jumps ~30 per period
  qui recast long t
  qui replace t    = mdy(1, 1, 2000) + t * 30
  qui replace ttre = mdy(1, 1, 2000) + ttre * 30
  format t ttre %td
}

if "`SHAPE'" == "onsetgap" {
  * cohort 1 is never observed in its own onset period
  qui drop if cohort == 1 & t == 2004
}
if "`SHAPE'" == "ed_ttrevary" {
  * onset varies WITHIN a unit: two onsets recorded for the same id
  qui replace ttre = ttre + 2 if cohort == 1 & t >= 2005
  qui replace h    = t - ttre
}
if "`SHAPE'" == "ed_ttremix" {
  * the classic data error: onset recorded only from the onset period onward
  qui replace ttre = . if h < 0
  qui replace h    = t - ttre
}
if "`SHAPE'" == "ed_dup" {
  * a merge accident: every row of cohort 2 appears twice
  qui expand 2 if cohort == 2
}
if "`SHAPE'" == "ed_singleton" {
  * one unit observed in a single period
  bysort unit (t): gen byte _keep1 = (_n == 1)
  qui drop if unit == 1 & _keep1 == 0
  drop _keep1
}
if "`SHAPE'" == "ed_tsingle" {
  * a calendar period in which exactly one unit is observed
  qui drop if t == 2009 & unit != 2
}
if "`SHAPE'" == "ed_onecell" {
  * exactly one unit per cohort: one observation per cohort-period cell
  bysort cohort (unit): gen byte _k1 = (unit == unit[1])
  qui keep if _k1 == 1
  drop _k1
}
if "`SHAPE'" == "ed_zerovar" {
  * y identically constant inside cohort 1, and ync constant everywhere
  qui replace y   = 4 if cohort == 1
  qui replace ync = 4
}
if "`SHAPE'" == "ed_outlier" {
  qui replace y   = y   + 1e12 if mod(unit, 37) == 0 & t == 2006
  qui replace ync = ync + 1e12 if mod(unit, 37) == 0 & t == 2006
}
if "`SHAPE'" == "ed_badw" {
  qui replace w = 0  if mod(unit, 11) == 0
  qui replace w = -1 if mod(unit, 13) == 0
  qui replace w = .  if mod(unit, 17) == 0
}
if "`SHAPE'" == "ed_collin" {
  * xcollin = 2*x1 + 3 is already in the data; add a cohort-constant control
  gen double xcohort = cohort
}
if "`SHAPE'" == "ed_ctrlconst" {
  gen double xcohort = cohort
  label var xcohort "control constant within cohort"
}
if "`SHAPE'" == "ed_ctrlistreat" {
  gen byte xpost = post
  label var xpost "control IS the treatment indicator"
}
if "`SHAPE'" == "ed_never1" {
  * exactly one never-treated unit survives
  qui levelsof unit if ttre >= ., local(_nv)
  local _keepnv : word 1 of `_nv'
  qui drop if ttre >= . & unit != `_keepnv'
}
if "`SHAPE'" == "ed_nopre" {
  * every treated cohort starts at the first sample period
  qui replace ttre = 2000 if ttre < .
  qui replace h    = t - ttre
  qui replace post = (ttre < . & t >= ttre)
}
if "`SHAPE'" == "ed_ifflag" {
  * a flag variable whose values are 0/1/2 -- "if flag" keeps ONLY the 1s
  gen byte flag = 1
  qui replace flag = 2 if mod(unit, 5) == 0
  qui replace flag = 0 if mod(unit, 7) == 0
  label var flag "values 0/1/2 - an if-flag restriction keeps only the 1s"
}
if "`SHAPE'" == "ed_float" {
  * float storage of the numeric payload: the fcol_sg precision class
  qui recast float y
  qui recast float ync
  qui recast float w
  qui recast float dose
  qui recast float x1
}
if "`SHAPE'" == "ed_tscale" {
  * t is a DAILY date but ttre was recorded as a MONTHLY date: scale mismatch
  qui recast long t
  qui replace t = mdy(1, 1, 2000) + (t - 2000) * 365
  format t %td
}

di as text "SHAPE " as result "`SHAPE'" as text ": `SHAPEDESC'"
