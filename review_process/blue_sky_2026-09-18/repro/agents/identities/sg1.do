clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
use "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/identities/base.dta", clear
cap gen byte one = 1
dis "SENT_A"
qui wooldid y i t ttre , att att_it att_i att_itime esprelength(2) espostlength(6) cluster(i)
mat A = e(b)
dis "SENT_B"
qui wooldid y i t ttre , att att_it att_i att_itime esprelength(2) espostlength(6) cluster(i) subgroups(one)
mat B = e(b)
mat list A, format(%21x)
mat list B, format(%21x)
local an : colnames A
local bn : colnames B
dis "A names: `an'"
dis "B names: `bn'"
mata:
a = st_matrix("A"); b = st_matrix("B")
na = tokens(st_local("an")); nb = tokens(st_local("bn"))
for (j=1;j<=cols(na);j++) {
  for (k=1;k<=cols(nb);k++) {
    if (nb[k] == na[j]+"_SG1_" | nb[k] == subinstr(na[j],"main_","main_SG1_") ) {
      printf("%s vs %s  reldiff=%g\n", na[j], nb[k], abs(a[1,j]-b[1,k])/max((abs(a[1,j]),1e-30)))
    }
  }
}
for (j=1;j<=cols(nb);j++) printf("Bname %g = %s = %21x\n", j, nb[j], b[1,j])
end
dis "SENT_END"
