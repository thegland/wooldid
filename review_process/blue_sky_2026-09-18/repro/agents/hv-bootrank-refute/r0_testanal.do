clear all
set more off
set processors 2
set seed 4242
include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-bootrank-refute/anal.mata"
mata:
  B = 80; K = 5
  A = rnormal(B, K, 0, 1)
  A[1..30, 4] = J(30,1,.)
  A[1..20, 5] = J(20,1,.)
  st_matrix("AA", A)
  ok = rowmissing(A) :== 0
  CC = select(A, ok)
  d = CC :- mean(CC)
  st_matrix("VV", (d'*d)/(rows(CC)-1))
  nav = colsum(A :< .)
  st_matrix("NAV", nav)
end
mat colnames AA = main_att c2 c3 c4 c5
mat colnames VV = main_att c2 c3 c4 c5
mat rownames VV = main_att c2 c3 c4 c5
mat colnames NAV = main_att c2 c3 c4 c5
mata: wd_analyze("AA", "VV", "NAV", 20, 100, "T")
dis "R0_DONE"
