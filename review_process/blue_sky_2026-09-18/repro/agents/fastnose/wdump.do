* wdump.do -- dump wooldid's e() surface, one line per item, exact hex for numbers.
* NOTE: e(scalars)/e(macros)/e(matrices) are EMPTY in this Stata for every estimation
* command (regress included), and Mata st_dir("e()",...) errors r(3300) here, so the
* member list is an explicit candidate list grepped out of the pinned ado.  MACRO
* candidates are written as strings only -- never evaluated with `local x = e(name)',
* which aborts the do-file on e(cmdline) ("wooldid not found", r(111)).
cap program drop wdump
program define wdump
  syntax , tag(string) file(string)

  local FAM att att_it att_i att_itime ame ame_it
  local SGV 0 1 2 3 4 5

  * ---------------- matrices ----------------
  local MCAND b V V_analytic V_boot b_boot
  local MCAND `MCAND' wdidmainresults pwsgs gteffects escomp histogramestimates summarystats
  local MCAND `MCAND' perm_p perm_ci perm_sd perm_n perm_kernel perm_se_shift perm_comparator_b
  local MCAND `MCAND' p ci_l ci_u
  local MCAND `MCAND' cfxatt_0p5 cfxatt_1 cfxatt_2 cfxame_0p5 cfxame_1 cfxame_2
  local MCAND `MCAND' cfxatt_post cfxatt_pre cfxame_post cfxame_pre
  foreach f of local FAM {
    local MCAND `MCAND' b_es_`f' V_es_`f' adversarial_line_`f' esband_`f'
  }
  foreach s of local SGV {
    local MCAND `MCAND' wdidsgresults_`s' escomp_s`s'
    foreach f of local FAM {
      local MCAND `MCAND' b_es_`f'_s`s' V_es_`f'_s`s' adversarial_line_`f'_s`s'
    }
  }
  local MCAND : list uniq MCAND

  * ---------------- numeric scalars ----------------
  local SCAND N N_clust Nclus N_dropped df_r df_r_analytic level rank
  local SCAND `SCAND' r2 r2adj r2_within r2_withinadj r2pseud
  local SCAND `SCAND' escomp_refN escomp_refW escomp_refncoh escomp_refnunits escomp_trtN escomp_trtW
  local SCAND `SCAND' joint_horiz joint_horiz_i joint_Choriz joint_Choriz_i
  local SCAND `SCAND' imdf_min imdf_max perm_analyticse perm_completed perm_reps
  local SCAND `SCAND' contreat_mean contreat_med contreat_sd ypre_mean ypre_sd ypre_med ypre_iqr
  local SCAND `SCAND' wdid_fcncell wdid_nptb wdid_nprat wdid_wprat
  foreach f of local FAM {
    local SCAND `SCAND' numpre_`f' es_ndegen_`f' adversarial_pval_`f' joint_pre_`f'
  }
  foreach s of local SGV {
    foreach f of local FAM {
      local SCAND `SCAND' numpre_`f'_s`s' adversarial_pval_`f'_s`s'
    }
    local SCAND `SCAND' joint_horiz_s`s' joint_horiz_i_s`s' joint_Choriz_s`s' joint_Choriz_i_s`s'
  }
  local SCAND : list uniq SCAND

  * ---------------- macros (string only) ----------------
  local KCAND cmd cmdline estimate_type vce vcetype clustvar properties depvar
  local KCAND `KCAND' esband_families esband_engine esband_level
  local KCAND `KCAND' wdid_vcvstate wdid_vcespec wdid_fastnose wdid_fcol wdid_fcguard
  local KCAND `KCAND' wdid_fcskip wdid_sglist wdid_sgvar wdid_estimands wdid_espre wdid_esrel
  local KCAND `KCAND' wdid_semielas wdid_cetype wdid_timetrends wdid_poisson wdid_contreat
  local KCAND `KCAND' wdid_suppressjt wdid_suppresssgpe wdid_adversarial wdid_sumstats
  local KCAND `KCAND' wdid_sgdegen wdid_sgdegenact wdid_sgdegenr wdid_sgnopre
  local KCAND `KCAND' wdid_vblocks wdid_vblockn wdid_cfxmats wdid_cfxlatvar
  local KCAND `KCAND' perm_method perm_pconv anticipation
  local KCAND : list uniq KCAND

  tempname fh
  file open `fh' using "`file'", write append text
  local nseen 0
  foreach nm of local MCAND {
    cap confirm matrix e(`nm')
    if _rc != 0 continue
    local ++nseen
    tempname MX
    matrix `MX' = e(`nm')
    local nr = rowsof(`MX')
    local nc = colsof(`MX')
    local rn : rowfullnames `MX'
    local cn : colfullnames `MX'
    file write `fh' `"`tag'|MATDIM|`nm'|`nr' `nc'"' _n
    file write `fh' `"`tag'|MATRN|`nm'|`rn'"' _n
    file write `fh' `"`tag'|MATCN|`nm'|`cn'"' _n
    forval r = 1/`nr' {
      forval c = 1/`nc' {
        local v = `MX'[`r',`c']
        if missing(`v') file write `fh' "`tag'|E|`nm'[`r',`c']|MISS" _n
        else            file write `fh' "`tag'|E|`nm'[`r',`c']|" %21x (`v') _n
      }
    }
    matrix drop `MX'
  }
  foreach nm of local SCAND {
    cap confirm matrix e(`nm')
    if _rc == 0 continue
    local sv `"`e(`nm')'"'
    if `"`sv'"' == "" continue
    local ++nseen
    local nv = e(`nm')
    if missing(`nv') file write `fh' "`tag'|NV|`nm'|MISS" _n
    else             file write `fh' "`tag'|NV|`nm'|" %21x (`nv') _n
  }
  foreach nm of local KCAND {
    local sv `"`e(`nm')'"'
    if `"`sv'"' == "" continue
    local ++nseen
    file write `fh' `"`tag'|SV|`nm'|`sv'"' _n
  }
  file write `fh' `"`tag'|COUNT||`nseen'"' _n
  file close `fh'
end
