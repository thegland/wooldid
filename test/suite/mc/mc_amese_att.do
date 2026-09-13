* mc_amese_att.do — att-side supplement to the AME-SE grid: same designs, records the
* ATT quantities (main:att, pretreatment:att, ES_att pre5..post3) so the att-side SE
* correctness under each VCE is measured rather than assumed.
* INVENTORY NOTE (2026-09-13, critic C2): run_amese.sh runs mc_amese.do only, and no analyzer
* reads this driver's att_<case>_s*.csv (analyze_amese.py's `{case}_s*.csv` glob does not
* match it). Archived 2026-08-25 twin; reads only att rows of e(wdidmainresults), so it is
* unaffected by amescrub.
* args: 1=adopath 2=case 3=shard 4=nreps 5=baseseed 6=outdir
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
local ROWS "main:att pretreatment:att ES_att:pre5 ES_att:pre4 ES_att:pre3 ES_att:pre2 ES_att:post0 ES_att:post1 ES_att:post2 ES_att:post3"

tempname fh
file open `fh' using `"`OUTDIR'/att_`CASE'_s`SHARD'.csv"', write replace
local hdr "case,shard,rep,seed,rc"
forval j = 1/10 {
  local hdr "`hdr',b`j',s`j'"
}
file write `fh' "`hdr'" _n
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
  if `rc' == 0 {
    tempname W
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
    cap matrix drop `W'
  }
  local row "`CASE',`SHARD',`rep',`seed',`rc'"
  forval j = 1/10 {
    local vb : di %16.0g `b`j''
    local vs : di %16.0g `s`j''
    local row "`row',`=trim("`vb'")',`=trim("`vs'")'"
  }
  file write `fh' "`row'" _n
}
file close `fh'
dis "MC_AMESE_ATT_SHARD_DONE `CASE' `SHARD'"
exit, clear
