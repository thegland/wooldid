* ---------------------------------------------------------------------------
* gen_ct.do -- shared continuous-treatment generator for the fp-contreat probes
* Creates an EXACT-FIT panel nested in wooldid's continuous-treatment model:
*
*   y = a_i + g_t + delta_term(dose) + 1{t>=ttre}*( kappa_it + beta_it*dose )
*
* with NO residual, so a correctly specified fit recovers beta_it exactly and
* the true post-treatment AME is the (weight-)mean of beta_it over treated
* post rows.  `dnuis' selects which nuisance dose term is planted:
*     none  basic  byTime  byCohort  byTreated  byTtre
* args:  1 = dnuis (default basic)
* leaves: i (cohort group), unit, t, ttre, dose, y, betatrue, w
* ---------------------------------------------------------------------------
args dnuis
if "`dnuis'" == "" local dnuis basic

clear
local NG   = 4          /* cohort groups                                   */
local NU   = 8          /* units per group                                 */
local NT   = 8          /* periods                                         */
set obs `=`NG'*`NU'*`NT''
gen long unit = ceil(_n/`NT')
gen int  t    = _n - (unit-1)*`NT'
gen byte i    = ceil(unit/`NU')
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 6 if i == 2
* i = 3, 4 never treated

* dose: strictly positive, varies WITHIN every (i,t) cell
gen double dose = 1 + 2*runiform()

* sampling weight, non-degenerate
gen double w = 0.5 + runiform()

* cell-specific true intercept and slope, post-treatment only
gen double kappatrue = 0
gen double betatrue  = 0
qui replace kappatrue = 0.3 + 0.11*i + 0.07*(t-ttre) if ttre < . & t >= ttre
qui replace betatrue  = 0.8 + 0.25*i - 0.13*(t-ttre) if ttre < . & t >= ttre

* nuisance dose function g() -- the thing contreatcontrols() must partial out
gen double gnuis = 0
if "`dnuis'" == "basic"     qui replace gnuis = 0.9*dose
if "`dnuis'" == "byTime"    qui replace gnuis = (0.4 + 0.15*t)*dose
if "`dnuis'" == "byCohort"  qui replace gnuis = (0.4 + 0.30*i)*dose
if "`dnuis'" == "byTreated" qui replace gnuis = 0.7*(ttre < .)*dose
if "`dnuis'" == "byTtre"    qui replace gnuis = (0.4 + 0.2*cond(ttre<.,ttre,0))*dose

gen double y = 2 + 0.6*i + 0.35*t + gnuis + kappatrue + betatrue*dose
