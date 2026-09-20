*=======================================================================
* gen_panel.do -- self-contained staggered-DiD panel generator for the
*                 hv-bootrank probes.  NO repo data, NO repo generator.
*
* Locals honoured (set them before `include'):
*   pN      units                      [200]
*   pT      periods                    [12]
*   pt0     first period               [2000]
*   ponsets onset list                 [2004 2006 2008]
*   pG      number of clusters (clu)   [40]   clusters NEST units
*   pnever  never-treated share        [0.40]
*   pysd    idiosyncratic sd           [0.30]
*   pseed   seed                       [77201]
*   ptau    per-cohort effect list     [0.5 0.8 1.1]
*   psmall  units in the LAST cohort (overrides equal split) [.]
*
* Variables out: unit t ttre clu y ync truetau post h
* DGP: y = a_unit + f_t + tau_cohort*post + e,  NO pre-trends
*      (so every pre-period coefficient has a TRUE value of exactly 0).
*=======================================================================
if "`pN'"      == "" local pN      200
if "`pT'"      == "" local pT      12
if "`pt0'"     == "" local pt0     2000
if `"`ponsets'"' == "" local ponsets 2004 2006 2008
if "`pG'"      == "" local pG      40
if "`pnever'"  == "" local pnever  0.40
if "`pysd'"    == "" local pysd    0.30
if "`pseed'"   == "" local pseed   77201
if `"`ptau'"'  == "" local ptau    0.5 0.8 1.1

clear
set seed `pseed'
set sortseed `pseed'

local K : word count `ponsets'
qui set obs `pN'
gen long unit = _n

* never-treated block at the end
local nnev = round(`pN' * `pnever')
local ntre = `pN' - `nnev'

gen long cohidx = .
if "`psmall'" != "" {
    * last cohort deliberately small
    local nrest = `ntre' - `psmall'
    local per   = floor(`nrest' / (`K' - 1))
    local done  = 0
    forvalues k = 1/`=`K'-1' {
        local lo = `done' + 1
        local hi = `done' + `per'
        if `k' == `K' - 1 local hi = `nrest'
        qui replace cohidx = `k' in `lo'/`hi'
        local done = `hi'
    }
    qui replace cohidx = `K' in `=`nrest'+1'/`ntre'
}
else {
    local per = floor(`ntre' / `K')
    forvalues k = 1/`K' {
        local lo = (`k'-1)*`per' + 1
        local hi = `k'*`per'
        if `k' == `K' local hi = `ntre'
        qui replace cohidx = `k' in `lo'/`hi'
    }
}

gen double ttre = .
forvalues k = 1/`K' {
    local ok : word `k' of `ponsets'
    qui replace ttre = `ok' if cohidx == `k'
}

* clusters NEST units: contiguous blocks of units
gen long clu = ceil(unit * `pG' / `pN')
qui replace clu = `pG' if clu > `pG'

gen double a_unit = rnormal()

expand `pT'
bysort unit: gen int t = `pt0' + _n - 1

* period shocks, common to all units
tempname fmat
mata: st_matrix("`fmat'", rnormal(`pT', 1, 0, .5))
gen double f_t = .
forvalues j = 1/`pT' {
    qui replace f_t = `fmat'[`j',1] if t == `pt0' + `j' - 1
}

gen int h    = t - ttre
gen byte post = (h >= 0) & !missing(h)

gen double truetau = 0
forvalues k = 1/`K' {
    local tk : word `k' of `ptau'
    qui replace truetau = `tk' if cohidx == `k' & post == 1
}

gen double ync = a_unit + f_t + truetau
gen double y   = ync + rnormal(0, `pysd')

keep unit t ttre clu y ync truetau post h
order unit t ttre clu y ync truetau post h
sort unit t
qui compress
qui count
dis as text "gen_panel: rows=" as result r(N) as text "  units=`pN'  T=`pT'  G=`pG'  onsets=`ponsets'"
