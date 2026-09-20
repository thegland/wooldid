* p08: two DOCUMENTED quantitative claims about the Dirichlet engine.
*  (1) sthlp:532 "the Bayesian bootstrap ... shares [the pairs cluster bootstrap's]
*      first two moments"  => se(bayes,equal) ~= se(standard) on the same data.
*  (2) sthlp:548 bootalpha(size) "shrinks the cluster-weight dispersion by roughly
*      the square root of the mean cluster size" => se(size)/se(equal) ~= sqrt(C/N).
*  (3) F001 (ado:14289-14312, fixed 2026-09-13): weighted unbalanced replications
*      must stay CENTRED on the point estimate. |mean(bootdist)-b| / sd(bootdist).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

local R 400
* dataset: UNWEIGHTED and WEIGHTED arms, C=30 clusters, unbalanced
foreach arm in unw wtd {
  local wt ""
  if "`arm'" == "wtd" local wt "[aw=bigpop]"
  foreach spec in "standard ." "bayes equal" "bayes size" {
    local md : word 1 of `spec'
    local al : word 2 of `spec'
    local alopt ""
    if "`al'" != "." local alopt bootalpha(`al')
    use p04_panel, clear
    cap noi wooldid y st year onset `wt', att cluster(st) bootstrap(`md') `alopt' reps(`R') rseed(808) dots(0)
    local rc = _rc
    matrix BB = e(b)
    matrix VV = e(V)
    matrix BD = e(bootdist)
    local b0 = BB[1,1]
    local se = sqrt(VV[1,1])
    mata {
      d = st_matrix("BD")[.,1]
      d = select(d, d :< .)
      printf("[`arm' `md' `al'] rc=`rc' n=%f  b=%12.8f  se(e(V))=%12.8f  sd(draws)=%12.8f  mean(draws)-b = %12.8f  => %6.3f boot-sd\n", rows(d), `b0', `se', sqrt(variance(d)), mean(d)-`b0', (mean(d)-`b0')/sqrt(variance(d)))
      st_local("sd_`arm'_`md'_`al'", strofreal(sqrt(variance(d)), "%21x"))
    }
  }
}
qui use p04_panel, clear
qui count
local NN = r(N)
qui bysort st: gen byte f1 = _n==1
qui count if f1
local CC = r(N)
dis "N=`NN' C=`CC'  sqrt(C/N)=" sqrt(`CC'/`NN')
foreach arm in unw wtd {
  local a = `sd_`arm'_bayes_equal'
  local b = `sd_`arm'_standard_.'
  local c = `sd_`arm'_bayes_size'
  dis "[`arm'] sd(bayes equal)/sd(standard) = " %6.4f `a'/`b'  "   (claim ~1)"
  dis "[`arm'] sd(bayes size)/sd(bayes equal) = " %6.4f `c'/`a' "   (claim ~ sqrt(C/N) = " %6.4f sqrt(`CC'/`NN') ")"
}
dis "PROBE_DONE"
