* genB.do -- Design B: ADMIN SHAPE. i = state (14), t = 1..7, several individual rows
* per (i,t) cell, sample weights varying WITHIN cell and across cells.
* This is the shape in which W_it / W_i / T_i normalizers actually do work and in which
* att != att_it and att_i != att_itime.
clear
set obs 14
gen int i = _n
gen double ui = rnormal()*2
gen int ttre = .
replace ttre = 3 if i <= 4
replace ttre = 5 if i > 4 & i <= 8
* i in 9..14 never treated
expand 7
bysort i: gen int t = _n
* ragged number of individual records per (i,t) cell: 2 to 5
gen int nrec = 2 + int(runiform()*4)
expand nrec
bysort i t: gen int r = _n
gen double relyr = t - ttre
gen double tau = 0
replace tau = 0.8 + 0.5*relyr if relyr >= 0 & relyr < .
gen double y = ui + 0.4*t + tau + rnormal()*0.6
* sample weights varying within cell
gen double w = 0.3 + runiform()*3
gen double cwm = 0.2 + runiform()*2
keep i t ttre y w cwm
order i t ttre y w cwm
sort i t
