* SKEPTIC probe of EC-3 / EC-5 / EC-8: failure sites, surviving e(), replay refusals
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31337
set sortseed 31337
which wooldid

clear
set obs 96
gen long id = _n
gen byte clu = ceil(id/4)
gen byte sgs = mod(id,2)                 /* short subgroup labels */
gen double sg8 = 10000000 + mod(id,2)    /* 8-digit subgroup labels */
gen double ui = rnormal()
gen double dose = 0.5 + runiform()
gen int ttre = .
replace ttre = 2004 if mod(id,4)==0
replace ttre = 2007 if mod(id,4)==1
replace ttre = 2010 if mod(id,4)==2
expand 13
bys id: gen int t = 2000 + _n - 1
gen double y = 1 + ui + 0.05*(t-2000) + rnormal()
replace y = y + 1.5 + 0.4*(t-ttre) if t>=ttre & ttre<.
save base_p4, replace

*========== EC-3: is the e() after the 8-digit r(198) IDENTICAL to the clean run? =====
* clean run with SHORT labels
cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sgs) jointtests adversarial cluster(clu)
dis "SHORT_RC=" _rc
mat bS = e(b)
mat VS = e(V)
local kS = colsof(bS)
dis "SHORT K=`kS' names: `: colnames bS'"
* run with 8-DIGIT labels (same partition of units, only relabelled)
cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg8) jointtests adversarial cluster(clu)
local rc8 = _rc
dis "LONG_RC=`rc8'"
dis "LONG e(cmd)=|" e(cmd) "|"
cap mat bL = e(b)
cap mat VL = e(V)
cap dis "LONG K=" colsof(bL) " names: `: colnames bL'"
* compare the SHARED coefficients bit-for-bit
cap mata:
  bs = st_matrix("bS"); bl = st_matrix("bL")
  ns = st_matrixcolstripe("bS"); nl = st_matrixcolstripe("bL")
  printf("EC3 cols short=%f long=%f\n", cols(bs), cols(bl))
  for (a=1; a<=cols(bs); a++) {
    nm = ns[a,1]+":"+ns[a,2]
    hit = 0
    for (bq=1; bq<=cols(bl); bq++) {
      nm2 = nl[bq,1]+":"+nl[bq,2]
      if (nm2 == nm) { hit=1
        if (bs[1,a] != bl[1,bq]) printf("EC3 DIFF %s  %21x vs %21x\n", nm, bs[1,a], bl[1,bq])
      }
    }
    if (!hit) printf("EC3 ONLY-IN-SHORT %s\n", nm)
  }
end
* is the surviving e() usable?
cap lincom main_att
dis "EC3 lincom rc=" _rc " b=" r(estimate) " se=" r(se)
cap estimates store AFTERLONG
dis "EC3 est store rc=" _rc
dis "EC3 LASTRC=|${WOOLDID_LASTRC}|"
cap noisily wooldid
dis "EC3 bare-replay rc=" _rc

*========== EC-8: estimates restore then bare replay ==========
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
use base_p4, clear
cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2)
dis "EC8 good rc=" _rc
cap estimates store GOODM
dis "EC8 store rc=" _rc
cap noisily wooldid y id t ttre, att subgroups(nosuchvarhere)
dis "EC8 fail rc=" _rc "  LASTRC=|${WOOLDID_LASTRC}|"
cap noisily estimates restore GOODM
dis "EC8 restore rc=" _rc " ecmd=|" e(cmd) "| K=" colsof(e(b)) " LASTRC=|${WOOLDID_LASTRC}|"
cap noisily wooldid
dis "EC8 replay-after-restore rc=" _rc
* control: does any other Stata command behave this way?
cap noisily reg y i.t
cap estimates store R1
cap noisily reg y i.nosuchvar
dis "EC8 reg fail rc=" _rc
cap noisily estimates restore R1
cap noisily reg
dis "EC8 reg replay rc=" _rc
* does a fresh successful wooldid clear it?
cap noisily wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2)
dis "EC8 refit rc=" _rc " LASTRC=|${WOOLDID_LASTRC}|"
cap noisily wooldid
dis "EC8 replay-after-refit rc=" _rc

*========== EC-5: which failure sites destroy e()? ==========
clear all
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
use base_p4, clear
cap prog drop setbase
prog def setbase
  cap noisily wooldid y id t ttre, att att_i esfixedbaseperiod esprelength(2) espostlength(2) cluster(clu)
  cap estimates store USERM
  matrix usermat = J(2,2,7)
  scalar userscal = 42
  dis "BASE rc=" _rc " ecmd=|" e(cmd) "| K=" colsof(e(b)) " N=" e(N)
end
cap prog drop chkafter
prog def chkafter
  args tag
  dis "`tag' AFTER: ecmd=|" e(cmd) "| K=" colsof(e(b)) " N=" e(N) " Nclust=" e(N_clust)
  cap estimates describe USERM
  dis "`tag'   USERM rc=" _rc
  cap confirm matrix usermat
  dis "`tag'   usermat rc=" _rc
  cap confirm scalar userscal
  dis "`tag'   userscal rc=" _rc
  cap count if e(sample)
  dis "`tag'   e(sample) N=" r(N)
  cap noisily wooldid
  dis "`tag'   bare-replay rc=" _rc
  macro drop WOOLDID_LASTRC
end

setbase
cap noisily wooldid y id t ttre, att contreat(dose) esfixedbaseperiod esprelength(2) espostlength(2) cluster(clu)
dis "F_contreat rc=" _rc
chkafter F_contreat

setbase
cap noisily wooldid y id t ttre, att ccc_absorb(id) esfixedbaseperiod esprelength(2) espostlength(2) cluster(clu)
dis "F_cccabs rc=" _rc
chkafter F_cccabs

setbase
cap noisily wooldid y id t ttre, att ccc_absorb(id) cluster(clu)
dis "F_cccabs_noes rc=" _rc
chkafter F_cccabs_noes

setbase
cap noisily wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) subgroups(sg8) cluster(clu)
dis "F_sg8 rc=" _rc
chkafter F_sg8

setbase
cap noisily wooldid y id t ttre, att nosuchoption
dis "F_badopt rc=" _rc
chkafter F_badopt

setbase
cap noisily wooldid y id t nosuchvar, att
dis "F_missvar rc=" _rc
chkafter F_missvar

dis "ALLDONE_P4"
