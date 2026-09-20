* p21: why does the bayes percentile interval for att miss more often than the
* standard one when their standard errors differ by only 5%?
clear all
set more off
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
use p12_res, clear
gen long ord = _n
bysort arm draw (ord): gen byte k = _n
gen str12 est2 = est
replace est2 = cond(k==1, "att", "att_i") if arm=="bayes"
keep if (arm=="bayes" & est2=="att") | arm=="standard"
gen double wp = hi_p - lo_p
gen byte covp = (lo_p <= 1) & (hi_p >= 1)
gen byte misslo = (lo_p > 1)
gen byte misshi = (hi_p < 1)
table arm, stat(mean wp covp misslo misshi) stat(n b)
* paired: same 30 datasets
keep draw arm covp wp b lo_p hi_p
reshape wide covp wp b lo_p hi_p, i(draw) j(arm) string
gen byte disc = covpbayes != covpstandard
qui count if disc & covpstandard==1
dis "standard covers, bayes misses : " r(N)
qui count if disc & covpbayes==1
dis "bayes covers, standard misses : " r(N)
qui su wpbayes
local a = r(mean)
qui su wpstandard
dis "mean percentile width  bayes=" `a' "  standard=" r(mean) "  ratio=" `a'/r(mean)
qui corr bbayes bstandard
dis "point estimates identical across arms: " (r(rho) > 0.999999)
list draw bbayes lo_pbayes hi_pbayes lo_pstandard hi_pstandard if disc, noobs abbrev(12)
dis "PROBE_DONE"
