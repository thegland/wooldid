* shared generator: tiny staggered binary panel with a continuous dose
* args: none. leaves data in memory.
clear
set obs 60
gen long unit = _n
gen double ttre = .
replace ttre = 5 if unit <= 15
replace ttre = 7 if unit > 15 & unit <= 30
replace ttre = 9 if unit > 30 & unit <= 45
* units 46-60 never treated
gen double ufe = rnormal()
gen double dose = 0.5 + runiform()
expand 12
bysort unit: gen double t = _n
gen double tfe = 0.1*t
gen byte post = (t >= ttre) & ttre < .
gen double y = ufe + tfe + 2.0*post*dose + rnormal()*0.5
gen double i = unit
