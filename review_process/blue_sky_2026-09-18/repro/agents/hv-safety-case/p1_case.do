clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid
dis "BUILDCHECK_DONE"

do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-safety-case/mkdata.do"

* ---------- leg 1: lowercase baseline, must succeed ----------
dis "=== LEG1 safety(off) ==="
cap noi wooldid y i t ttre, att safety(off)
dis "LEG1_RC=" _rc
cap dis "LEG1_B=" el(e(b),1,1)

dis "=== LEG2 safety(on) ==="
cap noi wooldid y i t ttre, att safety(on)
dis "LEG2_RC=" _rc

dis "=== LEG3 no safety option ==="
cap noi wooldid y i t ttre, att
dis "LEG3_RC=" _rc

* ---------- case variants ----------
foreach v in OFF Off oFF ON On {
  dis "=== LEGCASE safety(`v') ==="
  cap noi wooldid y i t ttre, att safety(`v')
  dis "LEGCASE_`v'_RC=" _rc
}

* ---------- whitespace ----------
dis "=== LEGWS safety( off ) ==="
cap noi wooldid y i t ttre, att safety( off )
dis "LEGWS_RC=" _rc

* ---------- the ordering bite: a safety-keyed ERROR fires before validation ----------
dis "=== LEGORD safety(OFF) unconditionalse att_i ==="
cap noi wooldid y i t ttre, att_i unconditionalse safety(OFF)
dis "LEGORD_RC=" _rc

dis "=== LEGORD2 safety(off) unconditionalse att_i (control: should run) ==="
cap noi wooldid y i t ttre, att_i unconditionalse safety(off)
dis "LEGORD2_RC=" _rc

* ---------- does a CAUTION print before the rejection? ----------
dis "=== LEGCAU safety(OFF) unconditionalse aweight ==="
gen double wt = 1 + runiform()
cap noi wooldid y i t ttre [aw=wt], att unconditionalse safety(OFF)
dis "LEGCAU_RC=" _rc

dis "ALLDONE_P1"
