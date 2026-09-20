* What exactly does wooldid's replay gate reject?
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fuzz/data"
use "`D'/ds_canon.dta", clear
qui wooldid y unit t ttre, att att_i cluster(region)
dis "SETUP_RC=" _rc
dis "=== bare ==="
capture noisily qui wooldid
dis "R_BARE=" _rc
dis "=== trailing comma ==="
capture noisily qui wooldid ,
dis "R_COMMA=" _rc
dis "=== level(95) ==="
capture noisily qui wooldid , level(95)
dis "R_LEVEL=" _rc
dis "=== single space ==="
capture noisily qui wooldid  
dis "R_SPACE=" _rc
dis "=== estimates store + replay ==="
qui estimates store mm
capture noisily qui estimates replay mm
dis "R_REPLAY=" _rc
dis "=== estimates table (a common consumer) ==="
capture noisily estimates table mm
dis "R_ETAB=" _rc
dis "=== estimates stats ==="
capture noisily estimates stats mm
dis "R_ESTATS=" _rc
dis "=== estout/esttab if installed ==="
capture noisily esttab mm
dis "R_ESTTAB=" _rc
capture noisily coefplot mm
dis "R_COEFPLOT=" _rc
