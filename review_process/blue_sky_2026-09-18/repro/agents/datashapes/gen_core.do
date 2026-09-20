*==========================================================================
* gen_core.do -- parameterized staggered-DiD data generator for the
*                wooldid blue-sky review.  Built from scratch; uses NO
*                repo data and no repo generator.
*
* USAGE (locals in, data out; `include' runs in YOUR scope so the knobs
*        you set just before the include are the ones that are used):
*
*     local gN   200
*     local gT   12
*     local gonsets "2004 2006 2008"
*     include "<GENDIR>/gen_core.do"
*
* EVERY knob is optional.  Defaults are applied below with the
* "if empty then default" idiom, so an unset local means "default".
* CLEAR YOUR KNOBS between datasets (locals persist): the shape
* dispatcher gen_shapes.do does that for you.
*
* ------------------------- KNOBS ----------------------------------------
* gseed     master seed (set seed AND set sortseed)          [20260916]
* gN        number of units when gcounts is not given        [120]
* gcounts   EXACT unit count per treated cohort, e.g. "300 2 2 5"
*           (overrides gN/gnever/gshares)                    []
* gnevern   EXACT never-treated unit count (with gcounts)    []
* gshares   relative cohort sizes, e.g. "1 1 4"              [all 1]
* gnever    never-treated share of gN                        [0.40]
* gnevgrp   number of distinct never-treated `cohort' codes  [1]
* gT        number of calendar periods                       [10]
* gt0       first calendar period                            [2000]
* gtstep    spacing between observed periods (1=annual,
*           2=biennial, 5=quinquennial)                      [1]
* gonsets   onset periods, one per treated cohort            ["2004 2006 2008"]
* gdropt    calendar periods to delete from the whole panel  []
* gentry    P(unit enters late)                              [0]
* gentrymax max entry delay in periods                       [3]
* gexit     P(unit exits early)                              [0]
* gexitmax  max early-exit length                            [3]
* ggap      P(an interior unit-period row is missing)        [0]
* gpercell  rows per unit-period cell (1=panel, >1=RCS)      [1]
* gsub      households per cell (finer-than-i cluster) when
*           gpercell>1                                       [2]
* gysd      sd of the idiosyncratic error.  gysd==0 gives
*           EXACT recovery of the design truth               [0.30]
* gtsd      sd of the common period shock                    [0.50]
* gaisd     sd of the UNIT fixed effect a_i.  Set it to 0 to
*           make the DGP exact at i=cohort as well as at
*           i=unit (with unit effects present, aggregating to
*           cohorts on an unbalanced panel is an omitted-
*           variable problem of the USER's making, not a bug) [1]
* gtaumode  const | hetero | zero | ramp                     [hetero]
* gtaus     explicit per-cohort effect sizes (overrides)     []
* gtauhet   sd of a UNIT-level multiplier on the true effect;
*           >0 makes effects heterogeneous WITHIN a cohort, so
*           att_i under i=unit stops equalling att_i under
*           i=cohort.  0 leaves the RNG stream untouched.      [0]
* gantic    true effects start this many periods early       [0]
* gtrendsd  sd of unit-specific linear trends; >0 BREAKS
*           parallel trends (treated units trend up)         [0]
* gb1 gb2   coefficients on x1, x2 in y                      [0.5, -0.25]
* gdose     none | unitconst | cellconst | spike0 | skew |
*           lattice                                          [none]
* gwtype    none | unit | row | survey                       [none]
* gwspread  orders of magnitude the weights span             [1]
* gwzero    share of rows with weight exactly 0              [0]
* gwneg     share of rows with a NEGATIVE weight             [0]
* gwmiss    share of rows with a missing weight              [0]
* gymiss    share of rows with a missing y                   [0]
* gxmiss    share of rows with a missing x1                  [0]
* gpois     1 => also build count outcomes ycount, ycount0   [0]
*
* ------------------------- WHAT YOU GET ---------------------------------
* unit     long   unit id, 1..N                (finest panel id)
* cohort   long   cohort id: 1..K treated, K+1.. never-treated
* hh       long   household id INSIDE unit (finer than unit; only
*                 meaningful when gpercell>1, else == unit)
* region   long   coarse id CROSSING cohorts (5 regions)
* t        int    calendar period
* ttre     double onset period, missing for never-treated
* h        int    t - ttre (relative period), missing for never-treated
* post     byte   1 on rows the design treats as treated
* y        double outcome INCLUDING the x1,x2 channel
* ync      double outcome with NO covariate channel (use this for
*                 exact-recovery probes without controls())
* truetau  double the TRUE treatment effect on that row (0 if untreated)
* x1 x2    double covariates (x1 varies within cell, x2 unit-correlated)
* xcell    double covariate CONSTANT within a unit-period cell
* xcollin  double = 2*x1 + 3 (exactly collinear with x1)
* dose     double continuous treatment (0 pre-onset / never-treated)
* w        double sampling weight (all 1 unless gwtype)
* fwt      int    integer frequency-style weight (1..6), always built
* sgall sg2 sggap sgtiny sgbig sgmiss sgvary : subgroup menu (see below)
*                 sgall is constant 1 -- a one-level partition of the sample
* ycount ycount0 : count outcomes when gpois==1
*
* EXACTNESS HARNESS: with gysd==0 and gtrendsd==0 the DGP is exactly
* two-way fixed effects plus truetau, so a correctly specified wooldid
* run must reproduce the design truth to machine precision.  Use
* gen_truth.do to compute the four estimand truths.
*==========================================================================

* ---------------- defaults ----------------
if "`gseed'"    == "" local gseed    20260916
if "`gN'"       == "" local gN       120
if "`gnever'"   == "" local gnever   0.40
if "`gnevgrp'"  == "" local gnevgrp  1
if "`gT'"       == "" local gT       10
if "`gt0'"      == "" local gt0      2000
if "`gtstep'"   == "" local gtstep   1
if `"`gonsets'"' == "" local gonsets 2004 2006 2008
if "`gentry'"   == "" local gentry   0
if "`gentrymax'"== "" local gentrymax 3
if "`gexit'"    == "" local gexit    0
if "`gexitmax'" == "" local gexitmax 3
if "`ggap'"     == "" local ggap     0
if "`gpercell'" == "" local gpercell 1
if "`gsub'"     == "" local gsub     2
if "`gysd'"     == "" local gysd     0.30
if "`gtsd'"     == "" local gtsd     0.50
if "`gaisd'"    == "" local gaisd    1
if "`gtaumode'" == "" local gtaumode hetero
if "`gantic'"   == "" local gantic   0
if "`gtauhet'"  == "" local gtauhet  0
if "`gtrendsd'" == "" local gtrendsd 0
if "`gb1'"      == "" local gb1      0.5
if "`gb2'"      == "" local gb2      -0.25
if "`gdose'"    == "" local gdose    none
if "`gwtype'"   == "" local gwtype   none
if "`gwspread'" == "" local gwspread 1
if "`gwzero'"   == "" local gwzero   0
if "`gwneg'"    == "" local gwneg    0
if "`gwmiss'"   == "" local gwmiss   0
if "`gymiss'"   == "" local gymiss   0
if "`gxmiss'"   == "" local gxmiss   0
if "`gpois'"    == "" local gpois    0

clear
set seed `gseed'
set sortseed `gseed'

local K : word count `gonsets'

* ---------------- cohort sizes ----------------
if "`gcounts'" != "" {
  local KC : word count `gcounts'
  if `KC' != `K' {
    di as error "gen_core: gcounts has `KC' entries but gonsets has `K'"
    error 198
  }
  local ntre 0
  forval g = 1/`K' {
    local c`g' : word `g' of `gcounts'
    local ntre = `ntre' + `c`g''
  }
  if "`gnevern'" == "" local gnevern = max(1, round(`ntre' * `gnever' / (1 - `gnever')))
  local nnev `gnevern'
  local gN = `ntre' + `nnev'
}
else {
  if "`gshares'" == "" {
    local gshares
    forval g = 1/`K' {
      local gshares `gshares' 1
    }
  }
  local tot 0
  forval g = 1/`K' {
    local s`g' : word `g' of `gshares'
    local tot = `tot' + `s`g''
  }
  local ntre = round(`gN' * (1 - `gnever'))
  local acc 0
  forval g = 1/`K' {
    local c`g' = max(1, round(`ntre' * `s`g'' / `tot'))
    local acc = `acc' + `c`g''
  }
  local c`K' = max(1, `c`K'' + (`ntre' - `acc'))
  local ntre 0
  forval g = 1/`K' {
    local ntre = `ntre' + `c`g''
  }
  local nnev = max(1, `gN' - `ntre')
  local gN = `ntre' + `nnev'
}

* ---------------- unit-level frame ----------------
set obs `gN'
gen long unit = _n
gen double _u = runiform()
sort _u unit
gen long _rk = _n

gen long   cohort = .
gen double ttre   = .
local lo 1
forval g = 1/`K' {
  local hi = `lo' + `c`g'' - 1
  local on : word `g' of `gonsets'
  qui replace cohort = `g'   if inrange(_rk, `lo', `hi')
  qui replace ttre   = `on'  if inrange(_rk, `lo', `hi')
  local lo = `hi' + 1
}
qui replace cohort = `K' + 1 + mod(_rk - `lo', `gnevgrp') if cohort >= .
drop _u _rk
sort unit

* region: a COARSE grouping that crosses cohorts (for cluster() above i)
gen long region = mod(unit - 1, 5) + 1

* unit-level heterogeneity
gen double a_i     = rnormal() * `gaisd'
gen double trend_i = 0
if `gtrendsd' > 0 qui replace trend_i = `gtrendsd' * (0.5 * (ttre < .) + 0.5 * rnormal())

* unit-level dose draw (used by the unitconst / spike0 / skew forms)
gen double _dosedraw = exp(rnormal() * 0.6)
if "`gdose'" == "spike0"  qui replace _dosedraw = cond(runiform() < 0.40, 0, exp(rnormal() * 0.6))
if "`gdose'" == "skew"    qui replace _dosedraw = exp(rnormal() * 1.6)
if "`gdose'" == "lattice" qui replace _dosedraw = floor(runiform() * 6)

* unit-level effect heterogeneity (no draw at all when gtauhet == 0, so the
* RNG stream of every other shape is unchanged)
gen double _tauhet = 1
if `gtauhet' > 0 qui replace _tauhet = 1 + rnormal() * `gtauhet'

* unit-level entry / exit
gen int _ent = 0
gen int _ext = 0
if `gentry' > 0 qui replace _ent = ceil(runiform() * `gentrymax') if runiform() < `gentry'
if `gexit'  > 0 qui replace _ext = ceil(runiform() * `gexitmax')  if runiform() < `gexit'

* unit-level subgroup menu
gen byte sgall = 1
gen byte sg2   = (runiform() < 0.5)
gen int  sggap = 1
qui replace sggap = 3  if runiform() < 0.66
qui replace sggap = 7  if runiform() < 0.33
qui replace sggap = 99 if runiform() < 0.10
cap label drop sglbl
label define sglbl 1 "North" 3 "South" 7 "East" 99 "Outlying territories and possessions"
label values sggap sglbl
* sgtiny: level 2 holds a handful of TREATED units only
gen byte sgtiny = 1
qui replace sgtiny = 2 if ttre < . & runiform() < 0.03
* sgbig: 9-10 digit administrative codes (name-construction overflow surface)
gen double sgbig = 480000000 + 1000 * mod(unit, 7)
* sgmiss: a quarter of units carry no subgroup
gen byte sgmiss = sg2
qui replace sgmiss = . if runiform() < 0.25

* unit-level weight base
gen double _wbase = 1
if "`gwtype'" != "none" {
  local wsd = `gwspread' * ln(10) / 4
  qui replace _wbase = exp(rnormal() * `wsd')
}

* ---------------- expand to the panel ----------------
expand `gT'
bysort unit: gen int t = `gt0' + (_n - 1) * `gtstep'

* delete whole calendar periods from the panel
if "`gdropt'" != "" {
  foreach dt of local gdropt {
    qui drop if t == `dt'
  }
}

* entry / exit
qui drop if t <  `gt0' + _ent * `gtstep'
qui drop if t >  `gt0' + (`gT' - 1 - _ext) * `gtstep'

* interior holes (never the first or last row of a unit)
if `ggap' > 0 {
  bysort unit (t): gen byte _mid = (_n > 1 & _n < _N)
  qui drop if _mid == 1 & runiform() < `ggap'
  drop _mid
}

* ---------------- repeated cross-sections ----------------
if `gpercell' > 1 {
  expand `gpercell'
  bysort unit t: gen int _j = _n
  gen long hh = unit * 100 + mod(_j - 1, `gsub') + 1
  drop _j
}
else {
  gen long hh = unit
}

* ---------------- period effects, covariates ----------------
* _rid pins the row order so every by-group draw below is reproducible
sort unit t hh
gen long _rid = _n

sort t _rid
by t: gen double _gt = rnormal() * `gtsd' if _n == 1
by t: replace _gt = _gt[1]

sort unit t _rid
by unit t: gen double _xc = rnormal() if _n == 1
by unit t: replace _xc = _xc[1]
gen double xcell = _xc
drop _xc

gen double x1 = rnormal()
gen double x2 = 0.3 * a_i + 0.5 * rnormal()
gen double xcollin = 2 * x1 + 3

* ---------------- treatment timing and true effects ----------------
gen int  h    = t - ttre
gen byte post = (ttre < . & t >= ttre - `gantic')

gen double dose = 0
if "`gdose'" != "none" {
  qui replace dose = _dosedraw if post == 1
  if "`gdose'" == "cellconst" {
    * dose constant within a cohort-period cell: the state-policy design
    qui replace dose = 0
    qui replace dose = round(0.5 + 0.5 * cohort + 0.25 * h, 0.25) if post == 1
    qui replace dose = 0 if dose < 0
  }
}

gen double truetau = 0
forval g = 1/`K' {
  if "`gtaus'" != "" local tb`g' : word `g' of `gtaus'
  if "`gtaus'" == "" local tb`g' = 0.5 + 0.5 * `g'
  local ts`g' = 0.10 * (`g' - 2)
}
if "`gtaumode'" == "zero" {
  * placebo DGP: no effect anywhere (for size / coverage probes)
}
else {
  forval g = 1/`K' {
    if "`gtaumode'" == "const" qui replace truetau = 1                                    if post == 1 & cohort == `g'
    if "`gtaumode'" == "ramp"  qui replace truetau = `tb`g'' * (1 + max(0, h))            if post == 1 & cohort == `g'
    if "`gtaumode'" == "hetero" qui replace truetau = `tb`g'' + `ts`g'' * max(0, h)       if post == 1 & cohort == `g'
  }
  if "`gdose'" != "none" {
    qui replace truetau = truetau * dose
  }
  qui replace truetau = truetau * _tauhet
}

* ---------------- outcomes ----------------
gen double ync = a_i + trend_i * (t - `gt0') + _gt + truetau + rnormal() * `gysd'
gen double y   = ync + `gb1' * x1 + `gb2' * x2

if `gpois' == 1 {
  gen double _mu   = 0.4 * a_i + 0.2 * _gt + 0.3 * truetau + 1.2
  gen int ycount   = rpoisson(exp(_mu))
  gen int ycount0  = rpoisson(exp(_mu - 2.5))
  drop _mu
}

* ---------------- weights ----------------
gen double w = 1
if "`gwtype'" == "unit"   qui replace w = _wbase
if "`gwtype'" == "row"    qui replace w = _wbase * exp(rnormal() * 0.3)
if "`gwtype'" == "survey" qui replace w = round(_wbase * 250 * exp(rnormal() * 0.3))
if "`gwtype'" == "survey" qui replace w = 1 if w < 1
gen int fwt = 1 + floor(runiform() * 6)

if `gwzero' > 0 qui replace w = 0  if runiform() < `gwzero'
if `gwneg'  > 0 qui replace w = -w if runiform() < `gwneg'
if `gwmiss' > 0 qui replace w = .  if runiform() < `gwmiss'

* ---------------- injected missings ----------------
if `gymiss' > 0 qui replace y  = . if runiform() < `gymiss'
if `gymiss' > 0 qui replace ync = . if y >= .
if `gxmiss' > 0 qui replace x1 = . if runiform() < `gxmiss'

* sgvary: a subgroup that CHANGES within a unit over time (constancy screen)
gen byte sgvary = sg2
qui replace sgvary = 1 - sgvary if t > `gt0' + (`gT' / 2) * `gtstep'

drop _dosedraw _tauhet _ent _ext _wbase _gt _rid

* ---------------- labels and order ----------------
label var unit    "unit id (finest panel id)"
label var cohort  "cohort id (treated 1..K, then never-treated)"
label var hh      "household id inside unit (finer than unit)"
label var region  "coarse region id, crosses cohorts"
label var t       "calendar period"
label var ttre    "treatment onset period (. = never treated)"
label var h       "relative period t - ttre"
label var post    "design-treated row"
label var y       "outcome (with covariate channel)"
label var ync     "outcome (no covariate channel)"
label var truetau "TRUE treatment effect on this row"
label var dose    "continuous treatment dose"
label var w       "sampling weight"
label var fwt     "integer frequency-style weight"

order unit cohort hh region t ttre h post y ync truetau dose x1 x2 xcell xcollin w fwt
sort unit t hh
qui compress unit cohort hh region t fwt

* ---------------- report ----------------
qui count
local _nobs = r(N)
tempvar _tagu
qui bysort unit: gen byte `_tagu' = (_n == 1)
qui count if `_tagu' == 1
local _nun = r(N)
drop `_tagu'
sort unit t hh
di as text "gen_core: N=" as result `_nobs' as text " rows, " as result `_nun' as text " units, T=" as result `gT' ///
   as text ", cohorts=" as result `K' as text ", seed=" as result `gseed'
