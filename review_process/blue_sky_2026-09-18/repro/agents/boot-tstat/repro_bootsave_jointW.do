*==========================================================================
* REPRO: bootsave()'s _rep = 0 row for the jointW_* (and advW*) columns does
*        NOT hold the observed statistic.  It holds an observed event-study
*        COEFFICIENT, matched by column POSITION into a differently-shaped
*        matrix.
*
*  Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_bootsave_jointW.do
*
*  wooldid.sthlp (bootsave() entry): "Under boottstat the replication
*  joint-test statistics are added as jointW_name and the adversarial
*  statistic as advW, WITH THE OBSERVED STATISTICS ON THE _rep = 0 ROW."
*  The same entry sells the file as self-sufficient: "so the draws can be
*  recentred from the file alone".
*
*  Cause (wooldid.ado:14915-14923, wooldid_bootsavemk):
*      JO = st_matrix(jobs)
*      for (k = 1; k <= cols(JW); k++)  ... JO[1, k] ...
*  `JW' (= e(boot_jointW)) has ONE COLUMN PER FAMILY (att, att_it, ...).
*  `jobs' is built by wooldid_bootjobs (ado:14666) and has one column per
*  (family, PRE-HORIZON) pair holding e(b_es_<fam>)'s observed coefficient.
*  The two column sets are unrelated, and JO is indexed by JW's k.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"

clear
set obs 20
gen int state = _n
gen int onset = .
forval s = 1/15 {
  local w = mod(`s'-1, 3) + 1
  local o : word `w' of 2003 2004 2005
  qui replace onset = `o' in `s'
}
gen double a_i = rnormal()
expand 3
bys state: gen int unit = _n + 1000*state
gen double u_i = rnormal()*0.5
expand 8
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1 + 0.2*(t - onset)
qui replace tau = 0 if !post
gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
gen int ttre = onset

wooldid y state t ttre , cluster(state) bootstrap(standard) reps(199) ///
   boottstat rseed(88) att esprelength(3) espostlength(3) jointtests adversarial ///
   bootsave(dj) replace dots(0)

scalar P     = e(joint_horiz)
mat JW       = e(boot_jointW)
mat BES      = e(b_es_att)
dis _n "reported bootstrapped joint pre-trend p, e(joint_horiz) = " %12.9f P
mat list BES

preserve
use dj, clear
ds jointW_* advW*
qui sum jointW_att if _rep == 0
scalar DUMP0 = r(mean)
cap confirm variable advW
if _rc == 0 {
  qui sum advW if _rep == 0
  scalar ADV0 = r(mean)
  qui sum advW if _rep > 0
  dis "advW  _rep==0 = " %14.8f ADV0 "   replication min = " %14.8f r(min) "  max = " %14.8f r(max)
}
else {
  cap confirm variable advW1
  if _rc == 0 {
    qui sum advW1 if _rep == 0
    scalar ADV0 = r(mean)
    qui sum advW1 if _rep > 0
    dis "advW1 _rep==0 = " %14.8f ADV0 "   replication min = " %14.8f r(min) "  max = " %14.8f r(max)
  }
  else scalar ADV0 = .
}
restore
dis "ADV0_IS_NEGATIVE = " (ADV0 < 0)

mata:
  JW = st_matrix("JW"); BES = st_matrix("BES"); bcs = st_matrixcolstripe("BES")
  p  = st_numscalar("P")
  d0 = st_numscalar("DUMP0")
  w  = select(JW[.,1], JW[.,1] :< .)
  n  = rows(w)
  // the engine's rule:  p = (1 + #{W_r >= W_obs}) / (n + 1)
  kge = p*(n+1) - 1
  sw  = sort(w, -1)                       // descending
  printf("\nreplications with a usable W_r          n = %f\n", n)
  printf("implied #{W_r >= W_obs} from the printed p  = %f\n", kge)
  if (kge >= 1 & kge <= n-1)
     printf("so W_obs must lie in (%12.8f , %12.8f]\n", sw[kge+1], sw[kge])
  printf("bootsave _rep==0 value of jointW_att       = %12.8f\n", d0)
  printf("min / max of the replication W_r           = %12.8f / %12.8f\n", min(w), max(w))
  // what the dump value actually is:
  printf("\nobserved e(b_es_att) columns:\n")
  for (j = 1; j <= cols(BES); j++)
     printf("   %-10s %14.8f%s\n", bcs[j,2], BES[1,j],
            (abs(BES[1,j] - d0) < 1e-12 ? "   <== THIS is what the dump's _rep==0 jointW_att holds" : ""))
  // recompute the p the way the help says the file supports
  pmine = (1 + sum(w :>= d0))/(n+1)
  printf("\np recomputed from the file as documented   = %12.9f\n", pmine)
  printf("p reported by wooldid                      = %12.9f\n", p)
  printf("FAIL = %f\n", (abs(pmine - p) > 1e-10))
end
