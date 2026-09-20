clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"
cd "`WD'"
which wooldid
* Under bootdisplay(wald) the printed pointwise interval is b +/- invt(G-1)*sd
* (bootdf(clusters) is the DEFAULT with cluster()).  The simultaneous band is
* b +/- c*sd with c a raw bootstrap quantile that never sees bootdf().
* So c < invt(G-1,.975) makes the "simultaneous" band NARROWER than the
* SYMMETRIC pointwise interval it accompanies -- a case the help's non-nesting
* caveat explicitly does not cover (it is scoped to rank/asymmetric displays).
foreach GG in 8 12 20 {
  foreach MM in 1 3 {
    clear
    local gG `GG'
    local gUnits 3
    local gT 8
    local gnever = ceil(`GG'/4)
    local gonsets "2003 2004 2005"
    include "`WD'/gen_small.do"
    qui levelsof state, local(ss)
    cap noi wooldid y state t ttre , cluster(state) bootstrap(standard) reps(999) rseed(77) bootdisplay(wald) att esprelength(`MM') espostlength(`MM') simultaneousbands dots(0)
    local rc = _rc
    if `rc' != 0 {
      dis "W_G`GG'_M`MM'_RC=`rc'"
      continue
    }
    local cc = e(esband_crit_att)
    local df = e(boot_df)
    local zc = invt(`df', 0.975)
    dis "W_G`GG'_M`MM' : bootdf=`df'  zc=" %9.6f `zc' "  bandcrit=" %9.6f `cc' "  band_narrower=" (`cc' < `zc')
    local BFAM att
    local BTAG W_G`GG'_M`MM'
    include "`WD'/bandcmp.do"
  }
}
dis "P11_DONE"
