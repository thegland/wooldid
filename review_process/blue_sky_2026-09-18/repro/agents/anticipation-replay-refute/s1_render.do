clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---- A: how does syntax , Level(real 0) render its local? ----
cap prog drop s1p
prog def s1p
  syntax [, Level(real 0)]
  dis "SYNTAXLOCAL=[`level']"
end
foreach L in 99.9 99.99 99.5 92.5 98.9 79.9 99.91 88.8 44.4 12.3 95 {
  dis "--- L=`L'"
  s1p , level(`L')
}

* ---- B: how does `global g = c(level)' render? ----
foreach L in 99.9 99.99 99.5 92.5 98.9 79.9 99.91 88.8 44.4 12.3 95 {
  cap noi set level `L'
  dis "SETLEVEL `L' rc=" _rc " c(level)=" c(level)
  global G0 = c(level)
  dis "GLOBAL_RENDER [`L'] = [$G0]"
  local LL = c(level)
  dis "LOCAL_RENDER  [`L'] = [`LL']"
}
set level 95

* ---- C: does set level accept the 16-digit rendering? ----
foreach S in 99.90000000000001 99.99000000000001 99.5 92.5 {
  cap noi set level `S'
  dis "SETSTR [`S'] rc=" _rc
}
set level 95
dis "S1_DONE"
