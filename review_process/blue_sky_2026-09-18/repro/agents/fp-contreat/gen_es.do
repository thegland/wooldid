* exact-fit generator with a NONTRIVIAL pre-treatment dose slope inside the
* esprelength(2) window, for the event-study ame probes.
* leaves: i unit t ttre dose w relyr betatrue kappatrue y
clear
set obs 384            /* 6 groups x 8 units x 8 periods */
gen long unit = ceil(_n/8)
gen int  t    = _n - (unit-1)*8
gen byte i    = ceil(unit/8)
gen double ttre = .
qui replace ttre = 4 if i == 1
qui replace ttre = 5 if i == 2
qui replace ttre = 6 if i == 3
* i = 4,5,6 never treated
gen double relyr = t - ttre
gen double dose = 1 + 2*runiform()
gen double w    = 0.5 + runiform()
gen double kappatrue = 0
gen double betatrue  = 0
* pre-treatment (inside the esprelength(2) window = relyr -1 and -2 only; relyr <= -3
* is the POOLED REFERENCE and must stay 0 for the fit to be exact)
qui replace kappatrue =  0.20 if relyr == -1
qui replace betatrue  =  0.50 if relyr == -1
qui replace kappatrue = -0.15 if relyr == -2
qui replace betatrue  = -0.40 if relyr == -2
* post-treatment
qui replace kappatrue = 0.3 + 0.11*i + 0.07*relyr if relyr >= 0 & relyr < .
qui replace betatrue  = 0.8 + 0.25*i - 0.13*relyr if relyr >= 0 & relyr < .
gen double y = 2 + 0.6*i + 0.35*t + 0.9*dose + kappatrue + betatrue*dose
