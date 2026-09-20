* wddump.do -- helper: dump e() results of a wooldid run to a csv keyed name,hex,dec
* usage: local DUMPTAG <tag> ; local DUMPFILE <path> ; do wddump.do
cap prog drop wdd_mat
prog def wdd_mat
  syntax , handle(string) mname(string) label(string)
  cap confirm matrix `mname'
  if _rc != 0 {
    file write `handle' "MISSING,`label',,,,," _n
    exit
  }
  tempname M
  matrix `M' = `mname'
  local rn : rownames `M'
  local cn : colnames `M'
  local nr = rowsof(`M')
  local nc = colsof(`M')
  forval r = 1/`nr' {
    local rlab : word `r' of `rn'
    forval c = 1/`nc' {
      local clab : word `c' of `cn'
      file write `handle' "`label',`rlab',`clab',`r',`c'," (strofreal(`M'[`r',`c'],"%21x")) "," (strofreal(`M'[`r',`c'],"%22.17e")) _n
    }
  }
end

cap prog drop wdd_run
prog def wdd_run
  syntax , file(string)
  file open wdh using "`file'", write replace text
  file write wdh "block,row,col,ri,ci,hex,dec" _n
  tempname B
  cap matrix `B' = e(b)
  if _rc == 0 wdd_mat , handle(wdh) mname(`B') label(b)
  foreach m in b_es_att b_es_att_it b_es_att_i b_es_att_itime gteffects wdidmainresults escomp {
    tempname M`m'
    cap matrix `M`m'' = e(`m')
    if _rc == 0 wdd_mat , handle(wdh) mname(`M`m'') label(`m')
  }
  foreach s in N N_clust df_r {
    file write wdh "scalar,`s',,,," (strofreal(e(`s'),"%21x")) "," (strofreal(e(`s'),"%22.17e")) _n
  }
  file close wdh
end
