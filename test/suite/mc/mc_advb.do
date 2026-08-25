* mc_advb.do — adversarial line-test SIZE BOUNDARY driver (AB campaign;
* PREREG_ADVBOUND_2026-08-25). One shard of one config; one CSV row per rep.
*
* The question: the adversarial test is severely oversized under cohort-level
* clustering (A1/A2/A6 20.5-23.2%; D2 48.4% even with 8-unit cohorts) and correctly
* sized in exactly one measured configuration (A4 att 5.75% under robust with 6-unit
* cells). This driver maps the boundary by taking the panel geometry, the VCE, the
* error structure, the onset pattern and the planted-alternative size as ARGUMENTS.
*
* The null is TRUE by construction wherever PWR==0: y carries no treatment effect and
* no dose effect, so the true event-study profile is flat zero, which lies on a line.
* PWR>0 plants tau(h)=PWR*h^2 on treated units (the power arms; null FALSE by design).
*
* args: 1=adopath 2=config 3=shard 4=nreps 5=baseseed 6=outdir
*       7=G (cohorts) 8=U (units per cohort) 9=T (periods)
*       10=vce token: robust | cluster(i) | cluster(unitid) | cluster(cl2)
*       11=errtype: iid | ar1 | cohortshock
*       12=mode: binary | contreat
*       13=STAG (0 = common onset; 1 = onsets alternate ONSET+1/ONSET)
*       14=SCS  (cohort-period shock sd; used only when errtype==cohortshock)
*       15=PWR  (planted quadratic deviation; 0 = null true)
*       16,17=extra wooldid options (e.g. altmethod_adversarial, safety(off))
*
* NOTE the unit identifier: `unitid' is GLOBALLY unique ((i-1)*U + unit). The
* within-cohort counter 1..U repeats across cohorts, so clustering on it would
* silently cluster on U meaningless cross-cohort groups - and nothing would error.
* `cl2' gives exactly 2 sub-clusters per cohort (the m=2 rung of the ladder).
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local CFG `2'
local SHARD = `3'
local NREPS = `4'
local BASE = `5'
local OUTDIR `6'
local G = `7'
local U = `8'
local T = `9'
local VCE `10'
local ERR `11'
local MODE `12'
local STAG = `13'
local SCS = `14'
local PWR = `15'
local EXTRA `16' `17'

local ONSET = ceil(`T'/2)
local CONTOPT
if "`MODE'" == "contreat" local CONTOPT contreat(dose) contreatcontrols(basic)

* the reported horizon set is fixed by the geometry (verified: espre=5, espl=4, the
* relyr -1 reference is remapped out), so 8 att coefficients and q=6
local ROWS "ES_att:pre5 ES_att:pre4 ES_att:pre3 ES_att:pre2 ES_att:post0 ES_att:post1 ES_att:post2 ES_att:post3"
local AROWS "ES_ame:pre5 ES_ame:pre4 ES_ame:pre3 ES_ame:pre2 ES_ame:post0 ES_ame:post1 ES_ame:post2 ES_ame:post3"

tempname fh
file open `fh' using `"`OUTDIR'/`CFG'_s`SHARD'.csv"', write replace
local hdr "config,shard,rep,seed,rc,padv,pjoint,pame,pjointC,n,dfr"
forval j = 1/8 {
  local hdr "`hdr',b`j'"
}
forval j = 1/8 {
  local hdr "`hdr',s`j'"
}
if "`MODE'" == "contreat" {
  forval j = 1/8 {
    local hdr "`hdr',ab`j'"
  }
  forval j = 1/8 {
    local hdr "`hdr',as`j'"
  }
}
file write `fh' "`hdr'" _n

* shard-1 only: the realized V_es lower triangle, so a FAIL can be decomposed into a
* variance ratio by direction rather than argued about
tempname vfh
local VDUMP = (`SHARD' == 1)
if `VDUMP' {
  file open `vfh' using `"`OUTDIR'/`CFG'_V_s1.csv"', write replace
  file write `vfh' "config,rep,which,r,c,v" _n
}

local clcheck 0

forval rep = 1/`NREPS' {
  local seed = `BASE' + `SHARD'*10000 + `rep'
  set seed `seed'
  clear
  qui set obs `G'
  gen i = _n
  * even cohorts treated, odd never-treated (the comparison pool in every period)
  gen ttre = cond(mod(_n, 2) == 0, `ONSET' + `STAG'*mod(_n/2, 2), .)
  * cohort-level intercept: absorbed by the cohort FE, kept for comparability (inert)
  gen a = rnormal()
  qui expand `U'
  bysort i: gen unit = _n
  gen unitid = (i - 1)*`U' + unit
  gen cl2 = (i - 1)*2 + cond(unit <= ceil(`U'/2), 1, 2)
  * dose varies WITHIN each cohort-period cell (across units), so the ame is
  * identified and the SI-A4 void screen cannot fire
  gen dose = exp(0.8*runiform() + 0.2*mod(i, 5))
  qui expand `T'
  bysort i unitid: gen t = _n

  if "`ERR'" == "ar1" {
    bysort unitid (t): gen e = rnormal()
    bysort unitid (t): replace e = .7*e[_n-1] + sqrt(1 - .49)*rnormal() if _n > 1
  }
  else gen e = rnormal()
  * cohort-by-period common shock: correlated ACROSS units within a cohort-period,
  * which neither unit-level clustering nor robust can capture
  if "`ERR'" == "cohortshock" {
    tempvar cs
    bysort i t: gen `cs' = rnormal() if _n == 1
    bysort i t (`cs'): replace `cs' = `cs'[1]
    qui replace e = e + `SCS'*`cs'
  }

  gen y = a + .05*t + `PWR'*cond(ttre < ., (t - ttre)^2, 0) + e
  drop a e
  sort i unitid t

  * one-time assertion that the cluster variable spans the intended number of
  * clusters per treated cell - a silently-wrong cluster variable is exactly the
  * failure mode this driver was written to prevent
  if `clcheck' == 0 & "`VCE'" != "robust" {
    local cvar = subinstr(subinstr("`VCE'", "cluster(", "", .), ")", "", .)
    tempvar rel ctag cn
    qui gen `rel' = t - ttre if ttre < .
    qui gegen `ctag' = tag(i `rel' `cvar') if ttre < . & `rel' >= -5 & `rel' <= 3
    qui gegen `cn' = total(`ctag') if ttre < . & `rel' >= -5 & `rel' <= 3, by(i `rel')
    qui sum `cn'
    dis as result "CLCHECK `CFG' clustervar=`cvar' clusters-per-treated-cell min=" r(min) " max=" r(max)
    drop `rel' `ctag' `cn'
    local clcheck 1
  }

  ereturn clear
  cap qui wooldid y i t ttre, att `VCE' `CONTOPT' esprelength(4) espostlength(4) ///
      esfixedbaseperiod adversarial jointtests `EXTRA'
  local rc = _rc
  foreach nm in padv pjoint pame pjointC n dfr {
    local `nm' .
  }
  forval j = 1/8 {
    local b`j' .
    local s`j' .
    local ab`j' .
    local as`j' .
  }
  if `rc' == 0 {
    local padv = e(adversarial_pval_att)
    local pjoint = e(joint_horiz)
    local pame = e(adversarial_pval_ame)
    local pjointC = e(joint_Choriz)
    cap local n = e(N)
    cap local dfr = e(df_r)
    tempname W
    * e(wdidmainresults), NOT e(b_es_att)/e(V_es_att): those post under cap and can
    * silently fail to post, which would read as missing data rather than an error
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
      if "`MODE'" == "contreat" {
        local j 0
        foreach rn of local AROWS {
          local ++j
          local r = rownumb(`W', "`rn'")
          if `r' < . {
            local ab`j' = `W'[`r', 1]
            local as`j' = `W'[`r', 2]
          }
        }
      }
    }
    cap matrix drop `W'
    if `VDUMP' {
      foreach vk in att ame {
        tempname VM
        cap matrix `VM' = e(V_es_`vk')
        if _rc == 0 {
          local vr = rowsof(`VM')
          forval r = 1/`vr' {
            forval c = 1/`r' {
              local vv : di %18.0g `VM'[`r', `c']
              file write `vfh' "`CFG',`rep',`vk',`r',`c',`=trim("`vv'")'" _n
            }
          }
        }
        cap matrix drop `VM'
      }
    }
  }
  local row "`CFG',`SHARD',`rep',`seed',`rc'"
  foreach nm in padv pjoint pame pjointC n dfr {
    local v : di %18.0g ``nm''
    local row "`row',`=trim("`v'")'"
  }
  forval j = 1/8 {
    local v : di %18.0g `b`j''
    local row "`row',`=trim("`v'")'"
  }
  forval j = 1/8 {
    local v : di %18.0g `s`j''
    local row "`row',`=trim("`v'")'"
  }
  if "`MODE'" == "contreat" {
    forval j = 1/8 {
      local v : di %18.0g `ab`j''
      local row "`row',`=trim("`v'")'"
    }
    forval j = 1/8 {
      local v : di %18.0g `as`j''
      local row "`row',`=trim("`v'")'"
    }
  }
  file write `fh' "`row'" _n
}
file close `fh'
if `VDUMP' file close `vfh'
dis "MC_ADVB_SHARD_DONE `CFG' `SHARD'"
exit, clear
