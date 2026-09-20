*==============================================================================
* repro_AR2.do -- SELF-CONTAINED.
* Claim: under permute() (and bootstrap()) the replay path prints a blank line
* before the engine block that the live run does not, so replay is NOT
* byte-for-byte identical as sthlp:1458 states.
* Mechanism: ado:48-51 in the replay branch
*     if "`e(perm_method)'" != "" dis as result " "
*     wooldid_display , part(perm)
*     if "`e(boot_method)'" != "" dis as result " "
*     wooldid_display , part(boot)
* The live shell (ado:297, ado:507) calls the same part() with NO preceding blank.
*
* Run: stata-mp -e do repro_AR2.do   -- the do-file prints the offending hunk.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local HERE "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid

quietly {
    set obs 120
    gen long i = _n
    gen byte cohort = mod(_n,4)
    gen double ttre = .
    replace ttre = 2004 if cohort==1
    replace ttre = 2006 if cohort==2
    replace ttre = 2008 if cohort==3
    gen double a_i = rnormal()
    expand 12
    bysort i: gen int t = 1999 + _n
    gen double y = a_i + 0.1*(t-2000) + 0.8*(ttre<. & t>=ttre) + 0.3*rnormal()
}

cap prog drop LIVE
prog def LIVE
  wooldid y i t ttre, att cluster(i) espre(3) espost(3) permute(shuffle) permreps(25) rseed(77) dots(0) jointtests
end
cap prog drop REPLAY
prog def REPLAY
  wooldid
end

log using "`HERE'/repro_AR2_live.txt", replace text name(L1)
LIVE
log close L1
log using "`HERE'/repro_AR2_replay.txt", replace text name(L2)
REPLAY
log close L2

dis as result "==== the 3 lines before 'Randomization Inference' in each file ===="
cap prog drop SHOWTAIL
prog def SHOWTAIL
  args fn tag
  tempname fh
  local buf3 ""
  local buf2 ""
  local buf1 ""
  file open `fh' using "`fn'", read text
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', "Randomization Inference") > 0 {
      dis as result "`tag' -3: [" `"`macval(buf3)'"' "]"
      dis as result "`tag' -2: [" `"`macval(buf2)'"' "]"
      dis as result "`tag' -1: [" `"`macval(buf1)'"' "]"
      dis as result "`tag'  0: [" `"`macval(line)'"' "]"
      file close `fh'
      exit
    }
    local buf3 `"`macval(buf2)'"'
    local buf2 `"`macval(buf1)'"'
    local buf1 `"`macval(line)'"'
    file read `fh' line
  }
  file close `fh'
  dis as error "`tag': marker not found"
end
SHOWTAIL `HERE'/repro_AR2_live.txt LIVE
SHOWTAIL `HERE'/repro_AR2_replay.txt REPLAY
dis "REPRO_AR2_DONE"
