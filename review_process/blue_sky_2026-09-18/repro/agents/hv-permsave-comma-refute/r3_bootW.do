* r3_bootW.do - SKEPTIC's independent repro of hv-permsave-comma-02
* Which of boot_jointW / boot_jointWobs / boot_sgW / boot_sgWobs / boot_advW / boot_advWobs
* are ereturned, with and without boottstat, and what the permute() twin does.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 171717
set sortseed 171717
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w3"
which wooldid

clear
set obs 120
gen long unit = _n
gen byte grp = mod(_n, 5)
gen ttre = .
replace ttre = 2005 if grp == 1
replace ttre = 2008 if grp == 2
replace ttre = 2011 if grp == 3
gen byte sg = mod(_n, 2)
gen double ue = rnormal()*0.8
expand 14
bys unit: gen int year = 1999 + _n
gen double y = 1.2 + 0.05*(year-1999) + ue + rnormal()*0.5
replace y = y + (0.4 + 0.2*sg)*(year - ttre + 1)/3 if year >= ttre & ttre < .

prog def whichW
  syntax , tag(string)
  foreach m in boot_jointW boot_jointWobs boot_sgW boot_sgWobs boot_advW boot_advWobs ///
               perm_jointW perm_jointWobs perm_sgW perm_sgWobs perm_advW perm_advWobs {
    capture confirm matrix e(`m')
    dis "`tag'|`m'|" cond(_rc==0,"PRESENT","absent")
  }
  dis "`tag'|e(boot_jt_family)=" "`e(boot_jt_family)'"
  dis "`tag'|e(joint_horiz)=" e(joint_horiz)
  dis "`tag'|e(joint_att_allsgssame)=" e(joint_att_allsgssame)
  dis "`tag'|e(adversarial_pval_att)=" e(adversarial_pval_att)
end

* ---- RUN 1: boottstat (studentized)
capture noisily wooldid y unit year ttre, att cluster(unit) subgroups(sg) jointtests ///
   esprelength(3) espostlength(3) esfixedbaseperiod adversarial ///
   bootstrap(standard) reps(30) boottstat dots(0)
dis "R1_RC=" _rc
if _rc == 0 whichW, tag(R1_BOOTTSTAT)

* ---- RUN 2: same, no boottstat (fixedmetric)
capture noisily wooldid y unit year ttre, att cluster(unit) subgroups(sg) jointtests ///
   esprelength(3) espostlength(3) esfixedbaseperiod adversarial ///
   bootstrap(standard) reps(30) dots(0)
dis "R2_RC=" _rc
if _rc == 0 whichW, tag(R2_FIXEDMETRIC)

* ---- RUN 3: permute() twin
capture noisily wooldid y unit year ttre, att cluster(unit) subgroups(sg) jointtests ///
   esprelength(3) espostlength(3) esfixedbaseperiod adversarial ///
   permute(shuffle) permreps(30) dots(0)
dis "R3_RC=" _rc
if _rc == 0 whichW, tag(R3_PERMUTE)

dis "ALLDONE"
