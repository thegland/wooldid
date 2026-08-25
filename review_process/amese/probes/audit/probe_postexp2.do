* probe_postexp2 — POSTEXP addendum: symmetric-but-INDEFINITE V cases.
clear all
set more off
set linesize 255
cap program drop postexp
program define postexp, eclass
  args case
  tempname b V
  matrix `b' = (1, 2, 3)
  * S6: symmetric, indefinite via off-diagonals at tiny scale (eigs 3e-22, -1e-22, 1)
  if "`case'" == "S6" matrix `V' = (1e-22, 2e-22, 0 \ 2e-22, 1e-22, 0 \ 0, 0, 1)
  * S7: symmetric, indefinite at O(1) scale (eigs 3, -1, 1)
  if "`case'" == "S7" matrix `V' = (1, 2, 0 \ 2, 1, 0 \ 0, 0, 1)
  * S8: symmetric PSD with mixed scales incl an exact-zero block plus healthy block
  if "`case'" == "S8" matrix `V' = (0, 0, 0 \ 0, 0, 0 \ 0, 0, 1)
  matrix colnames `b' = x1 x2 x3
  matrix colnames `V' = x1 x2 x3
  matrix rownames `V' = x1 x2 x3
  ereturn post `b' `V'
end
foreach c in S6 S7 S8 {
  cap noi postexp `c'
  dis "POSTEXP2|`c'|rc=" _rc
  cap noi mata: d = diagonal(st_matrix("e(V)")); printf("POSTEXP2|`c'|eV diag: %12.0g %12.0g %12.0g|allzero=%g\n", d[1], d[2], d[3], all(st_matrix("e(V)") :== 0))
}
dis "PROBE_POSTEXP2_DONE"
