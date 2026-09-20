clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
set linesize 200
capture noisily {
  use "main.dta", clear
  qui wooldid y i t ttre, att
  local lev_def = e(b)[1,1]
  qui wooldid y i t ttre, att customvce(robust)
  local lev_exp = e(b)[1,1]
  qui wooldid y i t ttre, att unconditionalse robust
  local lev_unc = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity
  local sem_def = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity customvce(robust)
  local sem_exp = e(b)[1,1]
  qui wooldid y i t ttre, att semielasticity unconditionalse robust
  local sem_unc = e(b)[1,1]
  dis "B_LEV_DEF = " %21x `lev_def'
  dis "B_LEV_EXP = " %21x `lev_exp'
  dis "B_LEV_UNC = " %21x `lev_unc'
  dis "B_LEV_RATIO = " %20.12f `lev_def'/`lev_exp'
  dis "B_SEM_DEF = " %20.12f `sem_def'
  dis "B_SEM_EXP = " %20.12f `sem_exp'
  dis "B_SEM_UNC = " %20.12f `sem_unc'
  dis "B_SEM_RATIO = " %20.12f `sem_def'/`sem_exp'
  dis "B_SEM_PCTERR = " %20.6f 100*(`sem_def'/`sem_exp'-1)
}
dis "B_RC=" _rc
* G: dispersion sweep -- does the DEFAULT answer track the grand mean while the
* correct one tracks each unit's own level?  Widen alpha dispersion, hold the
* grand mean and the level effect fixed.
foreach s in 0 60 120 200 {
  capture noisily {
    use "main.dta", clear
    qui su alpha, meanonly
    local am = r(mean)
    qui su alpha
    qui replace alpha = `am' + (alpha-`am')*(`s'/r(sd))
    qui replace y = alpha + gamma + seceff + 25*post + 10*invnormal(uniform())
    qui count if y <= 0
    local bad = r(N)
    qui wooldid y i t ttre, att semielasticity
    local gd = e(b)[1,1]
    qui wooldid y i t ttre, att semielasticity customvce(robust)
    local ge = e(b)[1,1]
    dis "G_SD=" `s' "  nonpos=" `bad' "  def=" %16.9f `gd' "  exp=" %16.9f `ge' "  ratio=" %12.6f `gd'/`ge'
  }
  dis "G_RC_`s'=" _rc
}
dis "BLOCKB_DONE"
