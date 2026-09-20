* r10b_attrit.do -- SELF-CONTAINED.  NEW probe (mine): how many draws actually sit
* behind e(V_es_att) on a messy-but-plausible panel (few treated clusters, long pre
* window), and does e() record that count anywhere?  wooldid posts e(boot_completed)
* and e(boot_N_V) (listwise over ALL columns) but no per-family V_es count, and
* wooldid_ptp's reference distribution F(keff, e(df_r)) does not depend on it.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 40
gen state = _n
gen onset = .
replace onset = 12 + mod(_n,4) if inlist(_n, 3, 9, 14, 21, 27, 33)
expand 22
bysort state: gen t = _n
bysort state (t): gen double ueff = rnormal(0,.4) if _n==1
bysort state (t): replace ueff = ueff[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .6 + .08*rel, 0)
gen double y = 2 + .04*t + ueff + te + e
drop rel te
tab onset, m

local OPT "cluster(state) esfixedbaseperiod esprelength(9) espostlength(6) jointtests pretrendpower(80)"
cap noi wooldid y state t onset, `OPT'
dis "B_A_RC=" _rc
dis "BRES|ANALYTIC|keff=" e(roth_keff_att) "|npre=" e(roth_npre_att) "|prep=" %12.8g e(roth_prep_att) "|slope=" %12.8g e(roth_slope_att) "|se=" %12.8g e(roth_slope_se_att) "|mdt=" %12.8g e(roth_power_att)[1,2]
cap noi wooldid y state t onset, `OPT' bootstrap(standard) reps(200) rseed(313)
dis "B_200_RC=" _rc
dis "BRES|BOOT200|keff=" e(roth_keff_att) "|npre=" e(roth_npre_att) "|prep=" %12.8g e(roth_prep_att) "|slope=" %12.8g e(roth_slope_att) "|se=" %12.8g e(roth_slope_se_att) "|mdt=" %12.8g e(roth_power_att)[1,2] "|completed=" e(boot_completed) "|NV=" e(boot_N_V) "|rankjoint=" %12.8g e(joint_horiz) "|dfr=" e(df_r)
cap confirm matrix e(bootdist)
dis "BOOTDIST_RC=" _rc
matrix BD = e(bootdist)
mata:
  BD = st_matrix("BD")
  nm = st_matrixcolstripe("BD")[.,2]
  keep = J(0,1,.)
  for (b=1; b<=rows(nm); b++) if (substr(nm[b],1,7)=="ES_att_") keep = keep \ b
  C = BD[., keep]
  printf("VESN| reps rows = %f ; ES_att columns = %f ; LISTWISE COMPLETE = %f\n", rows(BD), rows(keep), sum(rowmissing(C):==0))
  keep2 = J(0,1,.)
  for (b=1; b<=rows(nm); b++) if (substr(nm[b],1,10)=="ES_att_pre") keep2 = keep2 \ b
  C2 = BD[., keep2]
  printf("VESN| pre columns = %f ; listwise complete over PRE only = %f\n", rows(keep2), sum(rowmissing(C2):==0))
  for (b=1; b<=rows(keep); b++) printf("VESCOL| %s available %f\n", nm[keep[b]], sum(BD[.,keep[b]] :< .))
end
dis "DONE"
