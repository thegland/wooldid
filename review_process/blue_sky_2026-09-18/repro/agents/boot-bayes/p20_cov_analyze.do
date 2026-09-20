* p20: re-analyse p12_res. The str8 `est' variable truncated "main_att_i" to
* "main_att", pooling the bayes arm's two estimands; postfile order recovers them.
* Also rebuild the analytic reference interval from b and se (the original probe
* posted the wrong two columns of e(wdidmainresults)).
clear all
set more off
set processors 2
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"
use p12_res, clear
gen long ord = _n
bysort arm draw (ord): gen byte k = _n
gen str12 est2 = est
replace est2 = cond(k==1, "att", "att_i") if arm=="bayes"
replace est2 = cond(k==1, "att", "att_i") if arm=="analytic"
local TRUE 1.0
gen byte covp = (lo_p <= `TRUE') & (hi_p >= `TRUE') if !missing(lo_p)
gen byte covw = (lo_w <= `TRUE') & (hi_w >= `TRUE') if !missing(lo_w)
* analytic: rebuild a normal interval from b, se
replace covw = (b - 1.959964*se <= `TRUE') & (b + 1.959964*se >= `TRUE') if arm=="analytic"
gen double widw = hi_w - lo_w
replace widw = 2*1.959964*se if arm=="analytic"
collapse (sum) ncovp=covp ncovw=covw (mean) mse=se mwid=widw (sd) sdb=b (count) n=b, by(arm est2)
gen double ratio = mse/sdb
list arm est2 n ncovp ncovw mse sdb ratio mwid, noobs sepby(arm) abbrev(12)
dis "ncovp/ncovw are COUNTS covered out of n, nominal 95%"
dis "PROBE_DONE"
