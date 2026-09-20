*=============================================================================
* p12 — SELF-CONTAINED REPRO (the mechanism).
*
* CLAIM.  After bootstrap(), wooldid posts a bootstrap covariance (e(V), e(V_es_*))
* but keeps e(df_r) on the ANALYTIC scale (bootdf(clusters) => G-1, here 39).  Every
* MULTI-DIMENSIONAL Wald statistic built on that pair is then referred to F(k, 39) as
* if the covariance were known, when in fact it is a sample covariance estimated from
* n = e(boot_completed) replications.  The correct reference is Hotelling's T-squared:
*        (n - k)/((n - 1) k) * W  ~  F(k, n - k).
* wooldid instead reports  Ftail(k, e(df_r), W/k), inflating the F statistic by
* (n - 1)/(n - k) and shrinking the p-value without bound as n -> k.
*
* TWO CONSUMERS ARE AFFECTED, and they are the SAME NUMBER:
*   (1) pretrendpower()'s "joint pre-trend test" line and e(roth_prep_att);
*   (2) a hand `test` over the pre-period coefficients, which sthlp:1340 invites.
* wooldid's OWN jointtests are NOT affected - they are rank tests on the replication
* draws - so one run prints two contradictory bootstrap pre-trend p-values.
*
* The Roth rank guard does NOT catch this: it fires only when invsym drops a dimension,
* i.e. when completed <= k.  At completed = k + 1 it is silent.
*
* Run: stata-mp -e do p12_hotelling.do   (needs the pinned build on the adopath)
*=============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

clear
set obs 40
gen int state = _n
gen int onset = .
replace onset = 9  if state <= 10
replace onset = 11 if state > 10 & state <= 20
replace onset = 13 if state > 20 & state <= 30
gen double sfe = rnormal()
expand 20
bysort state: gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double rely = t - onset
gen double te = cond(D, 0.30 + 0.05*rely, 0)
gen double y = 1 + sfe + 0.05*t + te + rnormal()*0.8

local PB cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) jointtests pretrendpower(80)

dis "===================== ANALYTIC ====================="
qui wooldid y state t onset, `PB'
matrix B0 = e(b)
local bn : colnames B0
local tcpre
foreach nm of local bn {
  if strpos("`nm'","ES_att_pre")>0 local tcpre `tcpre' (`nm' = 0)
}
dis "ANA  df_r          = " e(df_r)
dis "ANA  roth_prep_att = " e(roth_prep_att) "   (W=" e(roth_W_att) ", keff=" e(roth_keff_att) ")"
dis "ANA  joint_horiz   = " e(joint_horiz)
qui test `tcpre'
dis "ANA  hand test p   = " r(p) "  df=" r(df) " drop=" r(drop)

foreach R in 7 8 10 15 25 50 1000 {
  dis "===================== BOOTSTRAP reps(`R') ====================="
  qui wooldid y state t onset, `PB' bootstrap(standard) reps(`R') rseed(4242) dots(0)
  local n   = e(boot_completed)
  local k   = e(roth_keff_att)
  local W   = e(roth_W_att)
  local dfr = e(df_r)
  local rp  = e(roth_prep_att)
  local jh  = e(joint_horiz)
  qui test `tcpre'
  local hp  = r(p)
  local hdf = r(df)
  local hdrop = r(drop)
  * the Hotelling-corrected reference for the SAME statistic
  local Fh = .
  local ph = .
  if `n' > `k' {
    local Fh = (`n' - `k')/((`n' - 1)*`k') * `W'
    local ph = Ftail(`k', `n' - `k', `Fh')
  }
  dis "R`R'  n=" `n' "  k=" `k' "  df_r=" `dfr' "  W=" `W'
  dis "R`R'  wooldid roth_prep (F(k,df_r))   = " `rp'
  dis "R`R'  hand `test` p (same statistic)   = " `hp' "   df=" `hdf' " drop=" `hdrop'
  dis "R`R'  wooldid jointtests joint_horiz   = " `jh'
  dis "R`R'  HOTELLING-CORRECTED p F(k,n-k)   = " `ph' "   (inflation factor (n-1)/(n-k) = " (`n'-1)/(`n'-`k') ")"
  dis "R`R'  roth slope = " e(roth_slope_att) "  se = " e(roth_slope_se_att) "  CI = " e(roth_slope_ll_att) " .. " e(roth_slope_ul_att)
  cap noi matrix list e(roth_power_att), format(%12.6g)
}
dis "P12_DONE"
