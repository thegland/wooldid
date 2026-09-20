clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

* ---- mechanism probe: what does syntax's real do to level()? ----
cap prog drop lvlecho
prog def lvlecho
  syntax , [Level(real 0)]
  dis "ECHO_LEVEL=[`level']"
  dis "ECHO_LEN=" strlen("`level'")
end
foreach L in 95 99.5 97.5 90.25 99.9 99.99 99.95 68.27 {
  dis "---- input `L' ----"
  lvlecho , level(`L')
}

* ---- plain Stata set level ----
foreach L in 95 99.5 97.5 90.25 99.9 99.99 {
  cap noi set level `L'
  dis "SETLEVEL_`L'_RC=" _rc " c(level)=" c(level)
}
set level 95
