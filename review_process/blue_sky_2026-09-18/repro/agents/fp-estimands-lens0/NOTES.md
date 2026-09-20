# LENS 0 (independent reproduction) -- FP-SEMIELAS-FE
VERDICT: CONFIRMED, high, realistic. Reproduced from the claim alone on my own DGP.

DGP (mine, third independent design): 60 counties x 14 yrs, alpha_i ~ N(380,120),
calendar trend 6*(t-2008), time-varying 3-level sector (+60/-40), 3 adoption waves
(2012/2015/2018), 15 never-treated, level effect 25, sigma 10. y in 92..815, no y<=0.
Files: lens0_semielas.do/.log, lens0_blockB.do/.log, main.dta

## B  default vs customvce(robust) vs unconditionalse robust (plain att)
LEV  def +1.98d5e4734bd40X+004  exp/unc +1.98d5e4734bd21X+004   ratio 1.000000000000
SEM  def 0.050517350256  exp 0.051389837368  unc 0.051389837368  ratio 0.983022 (-1.70%)
  -> unconditionalse == customvce(robust) to the last digit; both != default.

## B2  the ONLY message printed on a default semielasticity run
cell-collapse caution, text ends "the run is slower, and no answer changes."
No warning about the denominator. (Body of that caution DOES say eydx "divides by
the prediction" -- but never that the prediction omits the absorbed FE.)

## C  homogeneous-level control (y = 400 + 25*post + noise): ratio 0.998492  -> gap collapses

## D  mechanism WITHOUT wooldid (reghdfe y post, absorb(i t)  vs  reg y i.i i.t post)
dydx  24.852437119968 both (identical to 12dp; b to ~1 ulp)
eydx  hdfe 0.053477596449   full 0.058589306995
hand oracles mean(b/xb): 0.053477596277 / 0.058589306995  -> reproduce margins exactly
corr(xb_hdfe, xb_full) = 0.206088;  SD 12.04 vs 127.96;  both means = 465.0333 = mean(y)
_cons = 455.7137 = mean(y) - b*mean(post).  Hand arithmetic
  0.625*(24.85244/455.714) + 0.375*(24.85244/480.566) = 0.0534790 == the default answer.
  => the default denominator is a SAMPLE-LEVEL CONSTANT, not any observation's prediction.

## E  *** NOVEL, not run by either prior agent ***  WITHIN THE DEFAULT ENGINE
Same reghdfe engine, same column space, two parameterizations of the sector effect.
reghdfe sanity: r2 = 0.99524061134268 for BOTH; ranks 3 vs 1; dummies not dropped.
 LEVEL   fe(sector) +1.98d5e4734bd51X+004 == controls(sec2 sec3) +1.98d5e4734bd51X+004  (BIT-IDENTICAL)
 SEMIEL  fe(sector) 0.054316723370  vs  controls(sec2 sec3) 0.054698927107  ratio 0.993013
 correct (controls + customvce(robust)) = 0.056368573561
=> the default engine disagrees WITH ITSELF under a pure reparameterization. Kills the
   "customvce is a different estimator" escape hatch.

## F  reach: esprelength(3) espostlength(3), def vs exp
 main_att         0.049335171 / 0.056769517  ratio 0.869
 pretreatment_att -0.005245170 / 0.000098570  SIGN FLIP
 ES_att_pre3      -0.002450753 / -0.003778026  0.649
 ES_att_pre2      -0.006479679 / +0.009506747  SIGN FLIP (both non-trivial vs att ~0.05)
 ES_att_pre1      -0.006805078 / -0.005433012  1.253
 ES_att_post0/1/2  0.941 / 0.708 / 0.921

## G  dispersion sweep (alpha SD rescaled, grand mean and effect held fixed)
 SD=0    def 0.050572343  exp 0.049430433  ratio 1.0231   <-- TIME FE ALONE suffice
 SD=60   0.044322340 / 0.043812296  ratio 1.0116
 SD=120  0.053990316 / 0.055898368  ratio 0.9659
 SD=200  0.048460868 / 0.061045129  ratio 0.7939
=> error is UNBOUNDED and TWO-SIDED. Two components: (1) the subpop's mean fitted level
   differs from the grand mean (any calendar trend does this; ratio>1), (2) Jensen curvature
   E[1/yhat] > 1/E[yhat] (unit-level dispersion; ratio<1). They can offset, which is why my
   main design shows only -1.7% while the prober saw +65%. No fixed factor, no sign guarantee.

## Code facts verified on the pinned build
ado:1476 reghdfe call has no `d`;  ado:1462/1468 ppmlhdfe DOES pass `d`.
`grep -n "xbd"` over wooldid.ado: ZERO hits -> nothing downstream restores the FE.
ado:1846-1849 poismarginstext is EMPTY on the linear path -> margins uses reghdfe's
default predict(xb).  Not in known_bydesign §1-6 (only BS-6, the poisson xb-scale remap).
