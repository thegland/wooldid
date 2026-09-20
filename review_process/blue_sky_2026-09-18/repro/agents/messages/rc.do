clear all
set more off
set processors 2
set seed 20260916
set sortseed 20260916
foreach rc in 1 9 100 111 119 198 301 322 411 459 498 503 602 907 {
  di "RCPROBE_BEGIN `rc'"
  cap noisily error `rc'
  di "RCPROBE_RC=" _rc
  di "RCPROBE_END `rc'"
}
di "SENTINEL_DONE"
