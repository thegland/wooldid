* mc_amese.do — MC-truth driver for the AME-SE diagnosis grid (one shard of one case).
* Per rep: fresh DGP draw, one wooldid run, record the DISPLAYED ame estimates and SEs
* (e(wdidmainresults) rows: main:ame, pretreatment:ame, ES_ame:pre5..post3) plus the
* V_es_ame eigen range. The analyzer compares cross-rep SDs of the estimates (truth)
* against mean reported SEs, per the preregistered ratio bands.
* args: 1=adopath 2=case (e.g. d1_basic_ci) 3=shard 4=nreps 5=baseseed 6=outdir
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
do "/Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh"
local CASE `2'
local SHARD = `3'
local NREPS = `4'
local BASE = `5'
local OUTDIR `6'

tokenize "`CASE'", parse("_")
local DC `1'
local CC `3'
local VC `5'
amese_opts `CC' `VC'
local OPTS `r(opts)'

local ROWS "main:ame pretreatment:ame ES_ame:pre5 ES_ame:pre4 ES_ame:pre3 ES_ame:pre2 ES_ame:post0 ES_ame:post1 ES_ame:post2 ES_ame:post3"

tempname fh
file open `fh' using `"`OUTDIR'/`CASE'_s`SHARD'.csv"', write replace
local hdr "case,shard,rep,seed,rc"
forval j = 1/10 {
  local hdr "`hdr',b`j',s`j'"
}
file write `fh' "`hdr',evmin,evmax,ndegen" _n

forval rep = 1/`NREPS' {
  local seed = `BASE' + `SHARD'*10000 + `rep'
  set seed `seed'
  amese_mkdata `DC'
  cap qui wooldid y i t ttre, `OPTS'
  local rc = _rc
  forval j = 1/10 {
    local b`j' .
    local s`j' .
  }
  local evmin .
  local evmax .
  local nd .
  if `rc' == 0 {
    local nd = e(es_ndegen_ame)
    tempname W V
    cap matrix `W' = e(wdidmainresults)
    if _rc == 0 {
      local j 0
      foreach rn of local ROWS {
        local ++j
        local r = rownumb(`W', "`rn'")
        if `r' < . {
          local b`j' = `W'[`r', 1]
          local s`j' = `W'[`r', 2]
        }
      }
    }
    cap matrix `V' = e(V_es_ame)
    if _rc == 0 {
      mata: ev = symeigenvalues(st_matrix(st_local("V"))); st_local("evmin", strofreal(min(ev), "%16.0g")); st_local("evmax", strofreal(max(ev), "%16.0g"))
    }
    cap matrix drop `W' `V'
  }
  local row "`CASE',`SHARD',`rep',`seed',`rc'"
  forval j = 1/10 {
    local vb : di %16.0g `b`j''
    local vs : di %16.0g `s`j''
    local row "`row',`=trim("`vb'")',`=trim("`vs'")'"
  }
  local ve1 = trim("`evmin'")
  local ve2 = trim("`evmax'")
  local vn : di %12.0g `nd'
  file write `fh' "`row',`ve1',`ve2',`=trim("`vn'")'" _n
}
file close `fh'
dis "MC_AMESE_SHARD_DONE `CASE' `SHARD'"
exit, clear
