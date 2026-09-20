* SELF-CONTAINED. Mechanism test for FW-01(b): is the weight EXPRESSION TEXT pasted
* verbatim into the estimand-weight arithmetic?
* Sharp asymmetric prediction if it is:
*   [aw = wf + zero]  ->  w_att_it = wf + zero/w_per_it = wf  ->  att_it == att EXACTLY
*   [aw = zero + wm]  ->  w_att_it = zero + wm/w_per_it      ->  IDENTICAL to [aw = wm]
* If instead the expression were materialised first, BOTH would behave like a plain
* weight variable and neither identity would hold.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

set obs 40
gen long i = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 4 if grp==1
replace ttre = 6 if grp==2
replace ttre = 8 if grp==3
gen double ai = rnormal()
expand 10
bys i: gen int t = _n
gen double y = ai + 0.3*t + rnormal() + 1.5*(t>=ttre & ttre<.)
gen double wf = 20 + 10*runiform()
gen double wm = 15 + 12*runiform()
gen double wtot = wf + wm
gen double zero = 0
gen double cwm = 0.5 + runiform()
assert wtot == wf + wm

prog def grab
  args tag
  tempname B
  matrix `B' = e(b)
  local names : colnames `B'
  dis "COLS_`tag': `names'"
  foreach c in att att_it att_i att_itime {
    local cn = colnumb(`B', "main_`c'")
    if `cn' < . {
      local v = `B'[1, `cn']
      dis "VAL_`tag'_`c' = " %21x (`v') "   dec " %20.15f (`v')
    }
  }
  tempname V
  matrix `V' = e(V)
  foreach c in att att_it att_i att_itime {
    local cn = colnumb(`V', "main_`c'")
    if `cn' < . {
      local s = sqrt(`V'[`cn', `cn'])
      dis "SE_`tag'_`c' = " %21x (`s') "   dec " %20.15f (`s')
    }
  }
end

dis _n "=== A: truth, [aw=wtot], customvce(robust) ==="
capture noisily wooldid y i t ttre [aw=wtot], att att_it att_i att_itime customvce(robust)
dis "A_RC=" _rc
grab A

dis _n "=== B: suspect, [aw=wf+wm], customvce(robust)  (wf+wm==wtot on every row) ==="
capture noisily wooldid y i t ttre [aw=wf+wm], att att_it att_i att_itime customvce(robust)
dis "B_RC=" _rc
grab B

dis _n "=== C: [aw=wf+zero] -- predict att_it == att bitwise ==="
capture noisily wooldid y i t ttre [aw=wf+zero], att att_it att_i att_itime customvce(robust)
dis "C_RC=" _rc
grab C

dis _n "=== D: [aw=wf] -- the honest comparison for C ==="
capture noisily wooldid y i t ttre [aw=wf], att att_it att_i att_itime customvce(robust)
dis "D_RC=" _rc
grab D

dis _n "=== E: [aw=zero+wm] -- predict IDENTICAL to F ==="
capture noisily wooldid y i t ttre [aw=zero+wm], att att_it att_i att_itime customvce(robust)
dis "E_RC=" _rc
grab E

dis _n "=== F: [aw=wm] ==="
capture noisily wooldid y i t ttre [aw=wm], att att_it att_i att_itime customvce(robust)
dis "F_RC=" _rc
grab F

dis _n "=== G: customweightmultiplier(cwm) + [aw=wf+wm] -- predict even att breaks ==="
capture noisily wooldid y i t ttre [aw=wf+wm], att customweightmultiplier(cwm) customvce(robust)
dis "G_RC=" _rc
grab G

dis _n "=== H: customweightmultiplier(cwm) + [aw=wtot] ==="
capture noisily wooldid y i t ttre [aw=wtot], att customweightmultiplier(cwm) customvce(robust)
dis "H_RC=" _rc
grab H

dis _n "=== I: unconditionalse leg, safety(off), [aw=wf+wm] vs [aw=wtot] ==="
capture noisily wooldid y i t ttre [aw=wf+wm], att att_it robust unconditionalse safety(off)
dis "I_RC=" _rc
grab I
capture noisily wooldid y i t ttre [aw=wtot], att att_it robust unconditionalse safety(off)
dis "J_RC=" _rc
grab J

dis "P2_DONE"
