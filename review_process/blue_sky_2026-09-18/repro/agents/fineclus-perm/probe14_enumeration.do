* probe14: is the enumerated randomization space EXACTLY the set of distinct multiset
* permutations of the observed cluster-level onsets? (hand enumeration in Mata)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"

* --- design A: 5 clusters, onsets {4,4,6,.,.} -> 5!/(2!1!2!) = 30
* --- design B: 6 clusters, onsets {4,4,6,6,.,.} -> 6!/(2!2!2!) = 90
foreach d in A B {
  if "`d'" == "A" {
    local file permA
    local nc 5
    local hand 30
  }
  else {
    local file B_NONE
    local nc 6
    local hand 90
  }
  use "`file'.dta", clear
  drop if _rep == 0
  gen str key = ""
  forval s = 1/`nc' {
    qui replace key = key + cond(missing(onset_`s'), "NT", string(onset_`s',"%9.0g")) + "|"
  }
  qui duplicates report key
  dis "ENUM `d': rows=" _N " distinct keys=" r(unique_value) " hand count=" `hand'
  * every key must be a permutation of the observed multiset (same counts of each onset)
  gen byte okmult = 1
  forval s = 1/`nc' {
    * counts
  }
  gen int n4 = 0
  gen int n6 = 0
  gen int nNT = 0
  forval s = 1/`nc' {
    qui replace n4 = n4 + (onset_`s' == 4)
    qui replace n6 = n6 + (onset_`s' == 6)
    qui replace nNT = nNT + missing(onset_`s')
  }
  qui count if n4 != cond("`d'"=="A",2,2) | n6 != cond("`d'"=="A",1,2) | nNT != cond("`d'"=="A",2,2)
  dis "ENUM `d': multiset violations=" r(N)
  * observed assignment present exactly once?
  if "`d'" == "A" local obskey "4|4|6|NT|NT|"
  else local obskey "4|4|6|6|NT|NT|"
  qui count if key == "`obskey'"
  dis "ENUM `d': observed assignment appears " r(N) " time(s)"
}
