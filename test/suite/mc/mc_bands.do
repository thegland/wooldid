* mc_bands.do - COVERAGE of the simultaneous sup-t band (simultaneousbands, PS-45) on the
* analytic-simulated and the bootstrap-studentized engines (Stage C scaffold, written 2026-09-13
* from the Stage A audit: known_findings KF-5/KF-6, mc_drivers F5 as amended by the critic).
*
* THE QUESTION. The band shipped 2026-09-10 after every MC had run; test_esband.do pins its
* structure only, and no MC has used a cluster strictly BELOW i - the geometry the help
* recommends and the only one on which the analytic e(V_es_*) (which the analytic band is
* simulated from) is the right object. This driver measures, per rep and engine, whether EVERY
* reported horizon's truth lies inside the band, against the pointwise intervals on the same run.
* Geometry: AB1's shape (advb_grid.json) - G cohorts of U units, T periods, cluster(unitid)
* (m = U clusters per treated cell), even cohorts treated at ONSET = ceil(T/2), iid errors.
* Planted effect: a CONSTANT tau on every treated post period (tau*(t >= ttre)), so with
* esfixedbaseperiod every reported horizon's truth is known EXACTLY: 0 on pre horizons, tau on
* post horizons (the band's row 3, relyr, says which is which).
* Arms (each a separate wooldid run on the same dataset):
*   an  att cluster(unitid) esprelength(4) espostlength(4) esfixedbaseperiod simultaneousbands
*       -> e(esband_engine) = analytic-simulated, band from e(V_es_att) on a fixed internal seed
*   bs  + bootstrap(standard) reps(B)            -> engine bootstrap (coefficient-scale band from
*       the replications; the only bootstrap engine admissible BELOW i)
*   bt  att cluster(i) ... bootstrap(standard) reps(B) boottstat -> engine bootstrap-studentized.
*       boottstat with cluster() REQUIRES cluster() to nest i (wooldid_bootctl, r(459): each
*       replication is estimated on a physical expansion of the resampled clusters), so this
*       arm runs at the cohort level on the same data - 15 treated cohorts of 6 units resampled
*       whole. It is the studentized engine's only admissible level on this DGP, and it is NOT
*       the SI-A1 containment case for the bootstrap (the resampled unit is the cohort).
* Recorded per row: the band matrix e(esband_att) read BY ROW POSITION (row 1 lb, row 2 ub,
* row 3 relyr; wooldid_esplots reads it the same way), e(esband_crit_att), e(esband_n_att)
* (bootstrap engines: listwise-complete replications), band_absent (rc 0 but no e(esband_att):
* the bootstrap kernel REFUSES the whole family on one degenerate horizon), band_all (every
* horizon's truth inside the band), band_hits, pw_all / pw_hits / pw_n from the pointwise
* lb/ub of e(wdidmainresults)'s ES_att rows (the interval of THAT engine: analytic t, the
* studentized interval under boottstat, the percentile interval under bs), and per-horizon hit
* flags for relyr -5..-2 and 0..3 (bd* = band, pw* = pointwise), plus ncomp/nfail.
*
* One CSV row per (rep, arm), appended as it completes.
* args: 1=ado 2=shard 3=nreps 4=baseseed 5=outdir 6=B (reps; production 199) 7=TAU (default 1)
*       8=arms (+-separated list, e.g. an+bs+bt - a comma would be read by `do` as the start
*       of options; default an+bs+bt) 9=G 10=U 11=T (defaults 30 6 12)
* Cost measured in the 2-rep smoke at B=19 (2026-09-13) and scaled to B=199: an ~2 s, bs and bt
* several minutes per rep (a full estimation per replication): shard 8 ways.
*
* ACCEPTANCE (to be preregistered before the production run):
*   - per engine (an and bs below i; bt at the cohort level), the band_all COUNT inside the
*     central ~95% binomial region of Binomial(M, .95) (analyze_mc.py's binom_region; M = 400
*     -> [371, 389]) - the band's nominal claim;
*   - pw_all share BELOW the band_all share (the multiplicity gap the band exists to close; with
*     8 horizons the pointwise family-wise coverage is well under .95 unless the horizons are
*     nearly perfectly correlated) and each per-horizon pointwise hit count inside Binomial(M, .95);
*   - integrity: rc = 0 and band_absent = 0 on every row (every horizon has a usable variance on
*     this geometry), ncomp = B on the bootstrap arms, nh = 8 (the AB1 window).
*
* USAGE (2026-09-13). Positional, shared with mc_ri_ext.do and mc_adv_bs12.do; the runner is
* run_stagec.sh, which honours the machine's Stata instance cap:
*     ./run_stagec.sh bands <outdir> <nreps> <B> <tau> <arms|-> <shard>...
* Direct invocation, if you need to vary the geometry (args 9-11 = G U T, default 30 6 12):
*     <stata> -e do mc_bands.do <ado> <shard> <nreps> <seedbase> <outdir> <B> <tau> <arms> [G U T]
* arms is an+bs+bt (analytic / bootstrap-percentile / bootstrap-studentized), + or , separated.
* B is the bootstrap replication count. Writes <outdir>/mc_bands_s<shard>.csv, one row per
* (rep, arm). SMOKE-RUN 2026-09-13 at nreps 1, B 20: rc 0, 3 rows, all three arms rc 0.
* NOT yet campaigned and there is no analyzer: this file is a scaffold, not evidence.
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local SHARD = `2'
local NREPS = `3'
local BASE  = `4'
local OUTDIR `5'
local B = `6'
local TAU `7'
if "`TAU'" == "" local TAU 1
local ARMS `8'
if "`ARMS'" == "" local ARMS an+bs+bt
local ARMS = subinstr(subinstr("`ARMS'", "+", " ", .), ",", " ", .)
local G `9'
local U `10'
local T `11'
if "`G'" == "" local G 30
if "`U'" == "" local U 6
if "`T'" == "" local T 12
local ONSET = ceil(`T'/2)

global MCBANDS_OUT `"`OUTDIR'/mc_bands_s`SHARD'.csv"'
global MCBANDS_SHARD = `SHARD'
global MCBANDS_TAU = `TAU'

*----------------------------------------------------------------- DGP (AB1 shape + constant tau)
cap program drop mkdgp
program define mkdgp
  args seed tau G U T onset
  clear
  set seed `seed'
  qui set obs `G'
  gen i = _n
  gen ttre = cond(mod(_n, 2) == 0, `onset', .)
  gen a = rnormal()
  qui expand `U'
  bysort i: gen unit = _n
  gen unitid = (i - 1)*`U' + unit
  qui expand `T'
  bysort i unitid: gen t = _n
  gen e = rnormal()
  gen y = a + .05*t + `tau'*(t >= ttre & ttre < .) + e
  drop a e
  sort i unitid t
  ereturn clear
end

*----------------------------------------------------------------- row writer
cap program drop emitrow
program define emitrow
  args arm rep seed rc secs
  local HZ m5 m4 m3 m2 0 1 2 3
  local FLD crit nband nh band_absent band_all band_hits pw_all pw_hits pw_n ncomp nfail
  foreach h of local HZ {
    local FLD `FLD' bd`h'
  }
  foreach h of local HZ {
    local FLD `FLD' pw`h'
  }
  foreach nm of local FLD {
    local `nm' .
  }
  local engine
  if `rc' == 0 {
    local engine `e(esband_engine)'
    cap local crit = e(esband_crit_att)
    cap local nband = e(esband_n_att)
    if "`arm'" != "an" {
      cap local ncomp = e(boot_completed)
      cap local nfail = e(boot_failed)
    }
    * ---- the band, by row position (row 1 lb, row 2 ub, row 3 relyr)
    tempname SB
    cap matrix `SB' = e(esband_att)
    if _rc != 0 | rowsof(`SB') < 3 {
      local band_absent 1
    }
    else {
      local band_absent 0
      local nh = colsof(`SB')
      local band_hits 0
      local band_all 1
      forval c = 1/`nh' {
        local rel = `SB'[3, `c']
        local truth = cond(`rel' < 0, 0, ${MCBANDS_TAU})
        local hit = (`SB'[1, `c'] - 1e-9 <= `truth' & `truth' <= `SB'[2, `c'] + 1e-9)
        local band_hits = `band_hits' + `hit'
        if !`hit' local band_all 0
        local tag = cond(`rel' < 0, "m" + string(-`rel'), string(`rel'))
        cap local bd`tag' = `hit'
      }
    }
    * ---- the pointwise intervals of the same run: e(wdidmainresults) rows ES_att:*, cols 5/6
    tempname MR
    cap matrix `MR' = e(wdidmainresults)
    if _rc == 0 {
      local pw_hits 0
      local pw_n 0
      local pw_all 1
      local rnames : rownames `MR'
      local reqs : roweq `MR'
      local nr = rowsof(`MR')
      local nc = colsof(`MR')
      forval r = 1/`nr' {
        local eq : word `r' of `reqs'
        if "`eq'" != "ES_att" continue
        local nm : word `r' of `rnames'
        local rel .
        if `nc' >= 7 local rel = `MR'[`r', 7]
        if `rel' >= . {
          if substr("`nm'", 1, 3) == "pre" local rel = -real(substr("`nm'", 4, .))
          if substr("`nm'", 1, 4) == "post" local rel = real(substr("`nm'", 5, .))
        }
        if `rel' >= . continue
        local truth = cond(`rel' < 0, 0, ${MCBANDS_TAU})
        local lo = `MR'[`r', 5]
        local hi = `MR'[`r', 6]
        if `lo' >= . | `hi' >= . continue
        local hit = (`lo' - 1e-9 <= `truth' & `truth' <= `hi' + 1e-9)
        local ++pw_n
        local pw_hits = `pw_hits' + `hit'
        if !`hit' local pw_all 0
        local tag = cond(`rel' < 0, "m" + string(-`rel'), string(`rel'))
        cap local pw`tag' = `hit'
      }
    }
  }
  local row "`arm',${MCBANDS_SHARD},`rep',`seed',${MCBANDS_TAU},`rc',`secs',`engine'"
  foreach nm of local FLD {
    local x = strtrim(string(``nm'', "%18.12g"))
    if "`x'" == "" local x .
    local row "`row',`x'"
  }
  tempname fh
  file open `fh' using "${MCBANDS_OUT}", write append
  file write `fh' `"`row'"' _n
  file close `fh'
end

*----------------------------------------------------------------- header (once per shard)
cap confirm file "${MCBANDS_OUT}"
if _rc {
  tempname fh0
  file open `fh0' using "${MCBANDS_OUT}", write replace
  file write `fh0' "arm,shard,rep,seed,tau,rc,secs,engine," ///
    "crit,nband,nh,band_absent,band_all,band_hits,pw_all,pw_hits,pw_n,ncomp,nfail," ///
    "bdm5,bdm4,bdm3,bdm2,bd0,bd1,bd2,bd3,pwm5,pwm4,pwm3,pwm2,pw0,pw1,pw2,pw3" _n
  file close `fh0'
}

local ESOPT esprelength(4) espostlength(4) esfixedbaseperiod simultaneousbands

*----------------------------------------------------------------- the loop
forval rep = 1/`NREPS' {
  local seed  = `BASE' + `SHARD'*10000 + `rep'
  local rseed = `BASE' + 500000 + `SHARD'*10000 + `rep'
  foreach arm of local ARMS {
    timer clear 91
    timer on 91
    mkdgp `seed' `TAU' `G' `U' `T' `ONSET'
    if "`arm'" == "an" {
      cap qui wooldid y i t ttre, att cluster(unitid) `ESOPT'
    }
    else if "`arm'" == "bt" {
      * cluster(i): boottstat refuses a cluster() finer than i (see the header)
      cap qui wooldid y i t ttre, att cluster(i) `ESOPT' ///
          bootstrap(standard) reps(`B') boottstat rseed(`rseed') dots(0)
    }
    else if "`arm'" == "bs" {
      cap qui wooldid y i t ttre, att cluster(unitid) `ESOPT' ///
          bootstrap(standard) reps(`B') rseed(`=`rseed' + 1000000') dots(0)
    }
    else {
      di as error "mc_bands: unknown arm `arm'"
      exit 198
    }
    local rc = _rc
    timer off 91
    qui timer list 91
    emitrow `arm' `rep' `seed' `rc' `=r(t91)'
    di "MCBANDS_PROGRESS rep=`rep' arm=`arm' rc=`rc' secs=" %6.1f r(t91)
  }
}
di "MCBANDS_SHARD_DONE `SHARD'"
exit, clear
