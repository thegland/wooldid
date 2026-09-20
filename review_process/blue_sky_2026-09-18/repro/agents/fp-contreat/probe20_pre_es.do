clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
local G "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/gen_es.do"

cap prog drop truths
prog def truths, rclass
  syntax varname , IFC(string) [WVAR(string)]
  tempvar q ww wperit wperi tperi tag
  if "`wvar'" == "" qui gen double `ww' = 1
  else qui gen double `ww' = `wvar'
  qui sum `varlist' [aw=`ww'] if `ifc'
  return scalar att = r(mean)
  qui gegen double `wperit' = total(`ww') if `ifc', by(i t)
  qui gen double `q' = `ww'/`wperit' if `ifc'
  qui sum `varlist' [aw=`q'] if `ifc'
  return scalar att_it = r(mean)
  qui gegen double `wperi' = total(`ww') if `ifc', by(i)
  qui replace `q' = `ww'/`wperi' if `ifc'
  qui sum `varlist' [aw=`q'] if `ifc'
  return scalar att_i = r(mean)
  qui gegen byte `tag' = tag(i t) if `ifc'
  qui gegen double `tperi' = total(`tag') if `ifc', by(i)
  qui replace `q' = `ww'/(`wperit'*`tperi') if `ifc'
  qui sum `varlist' [aw=`q'] if `ifc'
  return scalar att_itime = r(mean)
end

foreach wt in "none" "aw" {
 foreach est in att att_it att_i att_itime {
   qui do "`G'"
   local wgt
   local wv
   if "`wt'" == "aw" {
     local wgt [aw=w]
     local wv w
   }
   truths betatrue , ifc(ttre<. & inlist(relyr,-1,-2)) wvar(`wv')
   local T_pre = r(`est')
   cap qui wooldid y i t ttre `wgt' , `est' contreat(dose) contreatcontrols(basic) robust ///
        esprelength(2) espostlength(3) regtol(14)
   if _rc == 0 {
     mat M = e(wdidmainresults)
     local nr = rowsof(M)
     local rn : rownames M
     local re : roweq M
     if "`est'" == "att_it" & "`wt'" == "none" {
       dis "LAYOUT_RN=`rn'"
       dis "LAYOUT_RE=`re'"
     }
     forval r = 1/`nr' {
       local nm : word `r' of `rn'
       local eq : word `r' of `re'
       if "`eq'" == "pretreatment" & substr("`nm'",1,3) == "ame" {
         dis "PRE|`wt'|`est'|est=" %16.10f M[`r',1] "|truth=" %16.10f `T_pre' "|reld=" %10.3e abs(M[`r',1]-`T_pre')/abs(`T_pre')
       }
       if "`eq'" == "ES_ame" {
         local h = M[`r',7]
         local hest = "`est'"
         if "`est'" == "att_i" | "`est'" == "att_itime" local hest att_it
         truths betatrue , ifc(ttre<. & relyr==`h') wvar(`wv')
         local Th = r(`hest')
         dis "ESH2|`wt'|`est'|relyr=`h'|est=" %16.10f M[`r',1] "|truth=" %16.10f `Th' "|reld=" %10.3e abs(M[`r',1]-`Th')/max(abs(`Th'),1e-12)
       }
     }
   }
   else dis "PRE|`wt'|`est'|RC=" _rc
 }
}
dis "DONE20"
