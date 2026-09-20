* probe13 -- MEMBERSHIP under permexact, read off the dumps probe03 wrote.
* The whole plain-share convention (and the "cannot print 0" floor) turns on the observed
* assignment's own enumerated replication carrying the observed statistic BIT-EXACTLY.
clear all
set more off
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics"
mata:
void membchk(string scalar bc, string scalar tag)
{
  real matrix B
  real colvector r0, io
  real rowvector b0
  real matrix br
  real scalar i, j
  B  = st_data(., bc)
  r0 = st_data(., "_rep")
  io = st_data(., "isobs")
  b0 = select(B, r0 :== 0)
  br = select(B, (io :== 1) :& (r0 :> 0))
  if (rows(br) == 0) {
    printf("  %s: NO replication carries the observed assignment\n", tag)
    return
  }
  for (i=1;i<=rows(br);i++) {
    printf("  %s: identical-assignment rep %g bitwise-equal to rep 0 = %g ; maxreldif=%g\n",
           tag, i, (br[i,.] == b0), max(reldif(br[i,.], b0)))
    for (j=1;j<=cols(b0);j++) {
      if (br[i,j] != b0[1,j]) printf("    col %g rep=%21x  rep0=%21x\n", j, br[i,j], b0[1,j])
    }
  }
}
end
foreach f in x_def x_abs x_ase x_ts {
  dis "================ `f' ================"
  use `f'.dta, clear
  qui ds onset_*
  local ons `r(varlist)'
  cap qui ds onset_* _rep se_*, not
  if _rc != 0 qui ds onset_* _rep, not
  local bc `r(varlist)'
  dis "rows=" _N "  onsetcols=`: word count `ons''  coefcols=`: word count `bc''"
  qui tostring `ons', replace force
  qui egen str100 akey = concat(`ons'), punct("|")
  qui levelsof akey if _rep == 0, local(obskey) clean
  qui count if _rep > 0 & akey == "`obskey'"
  dis "identical-assignment replications = " r(N)
  qui duplicates report akey if _rep > 0
  dis "distinct assignments among replications = " r(unique_value)
  gen byte isobs = (akey == "`obskey'")
  mata: membchk("`bc'", "`f'")
}
dis "PROBE13_DONE"
