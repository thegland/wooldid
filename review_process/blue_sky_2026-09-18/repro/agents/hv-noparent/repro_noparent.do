*=====================================================================
* repro_noparent.do -- SELF-CONTAINED.
* Claim: several wooldid options are accepted, do nothing at all, and
* say nothing at all, when their parent option is absent.  This
* contradicts wooldid.sthlp line 1181 ("No accepted option is a silent
* no-op").  No repo data, no test suite, no helper files.
*
* Run:  stata-mp -e do repro_noparent.do   then read repro_noparent.log
* Verdict lines start with  ===VERDICT
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set maxvar 20000
which wooldid          // must print the pinned build, version 1.4

* ---- dump the FULL e() in exact hex, minus the command-line echo ----
cap prog drop edump
prog def edump
  args fn
  tempname fh
  file open `fh' using "`fn'", write replace text
  local ss : e(scalars)
  foreach s of local ss {
    file write `fh' "S `s' = `= string(e(`s'), "%21x")'" _n
  }
  local ms : e(macros)
  foreach m of local ms {
    if "`m'" != "cmdline" file write `fh' `"M `m' = `e(`m')'"' _n
  }
  local xs : e(matrices)
  foreach x of local xs {
    tempname M
    matrix `M' = e(`x')
    forval r = 1/`= rowsof(`M')' {
      local line
      forval c = 1/`= colsof(`M')' {
        local line `line' `= string(`M'[`r',`c'], "%21x")'
      }
      file write `fh' "X `x' `r':`line'" _n
    }
  }
  file close `fh'
end

* ---- run once: named log (so wooldid's own output is isolated) + e() dump
cap prog drop R
prog def R
  syntax , nm(string) opts(string asis) [yv(string) wt(string asis)]
  if "`yv'" == "" local yv y
  cap log close cx
  qui log using "R_`nm'.log", text replace name(cx)
  cap noi wooldid `yv' unit t ttre `wt', `opts'
  local rc = _rc
  qui log close cx
  cap noi edump "R_`nm'.txt"
  dis "  rc(`nm') = `rc'"
end

* ---- compare a pair and print a verdict ----------------------------
cap prog drop V
prog def V
  args tag a b
  tempname h
  local same 1
  local msg 1
  qui checksum "R_`a'.txt"
  local ca = r(checksum)
  local na = r(filelen)
  qui checksum "R_`b'.txt"
  if r(checksum) != `ca' | r(filelen) != `na' local same 0
  qui checksum "R_`a'.log"
  local la = r(checksum)
  local ma = r(filelen)
  qui checksum "R_`b'.log"
  if r(checksum) != `la' | r(filelen) != `ma' local msg 0
  dis "===VERDICT `tag': results-identical=`same'  output-identical=`msg'"
end

*=====================================================================
* DATA 1 -- plain staggered binary panel (no dose)
*=====================================================================
clear
set obs 60
gen long unit = _n
expand 10
bysort unit: gen int t = 1999 + _n
gen double ttre = .
qui replace ttre = 2003 if unit <= 18
qui replace ttre = 2005 if unit > 18 & unit <= 36
set seed 111
gen byte sg2 = mod(unit,2)
gen double wgt = 0.5 + runiform()
gen double y = 0.07*unit + 0.03*(t-1999) + 0.8*(ttre < . & t >= ttre) + rnormal()*0.3
qui save "rp_bin.dta", replace

local B att robust
foreach z in pfx1 pfx2 {
  cap erase "`z'_ES.dta"
}

dis "--- 1. saveplots() with no plotting option at all ---"
use "rp_bin.dta", clear
R , nm(n1_base) opts(`B') wt([aw=wgt])
use "rp_bin.dta", clear
R , nm(n1_orph) opts(`B' saveplots(pfx1)) wt([aw=wgt])
V saveplots_no_plot_option n1_base n1_orph

dis "--- 2. makeplots + saveplots() but no event-study horizons ---"
use "rp_bin.dta", clear
R , nm(n2_base) opts(`B' makeplots) wt([aw=wgt])
use "rp_bin.dta", clear
R , nm(n2_orph) opts(`B' makeplots saveplots(pfx2)) wt([aw=wgt])
V saveplots_makeplots_no_horizons n2_base n2_orph

dis "--- 3. ccnoabsorb with no contreat() ---"
use "rp_bin.dta", clear
R , nm(n3_orph) opts(`B' ccnoabsorb) wt([aw=wgt])
V ccnoabsorb_no_contreat n1_base n3_orph

dis "--- 4. contreatpoly(3) with no contreat() ---"
use "rp_bin.dta", clear
R , nm(n4_orph) opts(`B' contreatpoly(3)) wt([aw=wgt])
V contreatpoly_no_contreat n1_base n4_orph

dis "--- 5. suppresssgprimaryeffects with no subgroups() ---"
use "rp_bin.dta", clear
R , nm(n5_orph) opts(`B' suppresssgprimaryeffects) wt([aw=wgt])
V suppresssg_no_subgroups n1_base n5_orph
dis "    (expect results-identical=0 only because e(wdid_suppresssgpe) echoes the option;"
dis "     every estimate below is identical -- see the diff of R_n1_base.txt / R_n5_orph.txt)"

*=====================================================================
* DATA 2 -- repeated cross-section with a continuous dose, so the
*           whole continuous-treatment machinery estimates
*=====================================================================
clear
set obs 60
gen long unit = _n
expand 10
bysort unit: gen int t = 1999 + _n
expand 4
gen double ttre = .
qui replace ttre = 2003 if unit <= 18
qui replace ttre = 2005 if unit > 18 & unit <= 36
set seed 2468
gen byte post = (ttre < . & t >= ttre)
gen double d = runiform()*10
gen double dose = cond(post, d, 0)
gen double wgt = cond(d > 6, 40, 1)
gen byte sg2 = mod(unit,2)
gen double y = 0.07*unit + 0.03*(t-1999) + (0.5 + 0.15*d)*post + rnormal()*0.3
gen byte evertreated = (ttre != .)
qui gquantiles latvar = dose, xtile nquantiles(4) by(evertreated)
qui replace latvar = 0 if latvar >= .
qui save "rp_cts.dta", replace

local C att robust contreat(dose) contreatcontrols(untreatedZero) lattice(byTreated 4)

use "rp_cts.dta", clear
R , nm(n6_base) opts(`C') wt([aw=wgt])

dis "--- 6. cfxplottypes() with no makecfxplots ---"
use "rp_cts.dta", clear
R , nm(n6_orph) opts(`C' cfxplottypes(post 2)) wt([aw=wgt])
V cfxplottypes_no_makecfxplots n6_base n6_orph

dis "--- 7. cfxplottypes(GARBAGE) with no makecfxplots: the option's own"
dis "       validator is inside the makecfxplots block, so nonsense passes ---"
use "rp_cts.dta", clear
R , nm(n7_orph) opts(`C' cfxplottypes(psot 1.5 nosuchvariable)) wt([aw=wgt])
V cfxplottypes_GARBAGE_no_makecfxplots n6_base n7_orph
dis "    control: the SAME garbage WITH makecfxplots must be refused r(198):"
use "rp_cts.dta", clear
cap noi wooldid y unit t ttre [aw=wgt], `C' makecfxplots cfxplottypes(psot 1.5 nosuchvariable)
dis "===VERDICT cfxplottypes_garbage_WITH_parent: rc=" _rc " (expected 198)"

dis "--- 8. cfxplotnose with no makecfxplots ---"
use "rp_cts.dta", clear
R , nm(n8_orph) opts(`C' cfxplotnose) wt([aw=wgt])
V cfxplotnose_no_makecfxplots n6_base n8_orph

dis "--- 9. makecfxplotsbysg with subgroups() but no makecfxplots ---"
use "rp_cts.dta", clear
R , nm(n9_base) opts(`C' subgroups(sg2)) wt([aw=wgt])
use "rp_cts.dta", clear
R , nm(n9_orph) opts(`C' subgroups(sg2) makecfxplotsbysg) wt([aw=wgt])
V makecfxplotsbysg_no_makecfxplots n9_base n9_orph

dis "--- 10. latticeignoreweights on a PRE-MADE lattice variable (the"
dis "        1-token lattice() form, where no quantile is ever computed) ---"
use "rp_cts.dta", clear
R , nm(n10_base) opts(att robust contreat(dose) contreatcontrols(untreatedZero) lattice(latvar)) wt([aw=wgt])
use "rp_cts.dta", clear
R , nm(n10_orph) opts(att robust contreat(dose) contreatcontrols(untreatedZero) lattice(latvar) latticeignoreweights) wt([aw=wgt])
V latticeignoreweights_premade_lattice n10_base n10_orph

dis "--- 11. latticeignoreweights with no lattice() at all ---"
use "rp_cts.dta", clear
R , nm(n11_base) opts(att robust contreat(dose) contreatcontrols(untreatedZero)) wt([aw=wgt])
use "rp_cts.dta", clear
R , nm(n11_orph) opts(att robust contreat(dose) contreatcontrols(untreatedZero) latticeignoreweights) wt([aw=wgt])
V latticeignoreweights_no_lattice n11_base n11_orph

dis "--- CONTROL: latticeignoreweights DOES work when the lattice is"
dis "    actually computed from quantiles.  Must NOT be identical. ---"
use "rp_cts.dta", clear
R , nm(n12_orph) opts(`C' latticeignoreweights) wt([aw=wgt])
V CONTROL_latticeignoreweights_generated_lattice n6_base n12_orph

*=====================================================================
* FILE EVIDENCE: saveplots() without its parent writes nothing,
*                saveplots() with its parent writes files.
*=====================================================================
foreach z in pfx3 pfx4 {
  foreach s in _ES.dta _ES_att.gph {
    cap erase "`z'`s'"
  }
}
use "rp_bin.dta", clear
cap noi wooldid y unit t ttre [aw=wgt], att robust esprelength(2) espostlength(2) makeplots saveplots(pfx3)
dis "  with-parent rc=" _rc
use "rp_bin.dta", clear
cap noi wooldid y unit t ttre [aw=wgt], att robust esprelength(2) espostlength(2) saveplots(pfx4)
dis "  no-parent   rc=" _rc
foreach z in pfx1 pfx2 pfx3 pfx4 {
  local got
  foreach s in _ES.dta _ES_att.gph {
    cap confirm file "`z'`s'"
    if _rc == 0 local got `got' `z'`s'
  }
  dis "===VERDICT files written under saveplots(`z'): `got'"
}
dis "REPRO_NOPARENT_COMPLETE"
