* ---- shared checker: include this AFTER adopath/use ----
mata:
void wdpcmp()
{
    real matrix B, C, D, Dm
    string matrix bs, cs
    real scalar k, j, md, n, d, nmiss
    string scalar mn

    B  = st_matrix("e(b)")
    bs = st_matrixcolstripe("e(b)")
    C  = st_matrix("e(perm_comparator_b)")
    cs = st_matrixcolstripe("e(perm_comparator_b)")
    md = 0
    mn = "none"
    n  = 0
    for (k=1; k<=cols(C); k++) {
        for (j=1; j<=cols(B); j++) {
            if (bs[j,1]==cs[k,1] & bs[j,2]==cs[k,2]) {
                if (C[1,k]<. & B[1,j]<.) {
                    d = abs(C[1,k]-B[1,j])
                    if (abs(B[1,j])>1e-12) d = d/abs(B[1,j])
                    n = n + 1
                    if (d>md) {
                        md = d
                        mn = cs[k,2]
                    }
                }
                j = cols(B) + 1
            }
        }
    }
    st_numscalar("wdp_ncmp", n)
    st_numscalar("wdp_maxrd", md)
    st_local("wdp_at", mn)
    D = st_matrix("e(permdist_b)")
    if (rows(D)>0) {
        nmiss = sum(D :>= .)
        Dm = editmissing(D, 0)
        st_numscalar("wdp_dr", rows(D))
        st_numscalar("wdp_dc", cols(D))
        st_numscalar("wdp_dm", nmiss)
        st_numscalar("wdp_s1", sum(Dm))
        st_numscalar("wdp_s2", sum(Dm:*Dm))
    }
    else {
        st_numscalar("wdp_dr", 0)
        st_numscalar("wdp_dc", 0)
        st_numscalar("wdp_dm", .)
        st_numscalar("wdp_s1", .)
        st_numscalar("wdp_s2", .)
    }
}
end

cap prog drop pchk
prog def pchk
  syntax anything(name=tag) [, OPTS(string asis) WT(string asis) IFQ(string asis) YV(string) REPS(integer 30) IV(string) CL(string)]
  if "`yv'" == "" local yv y
  if "`iv'" == "" local iv unit
  if "`cl'" == "" local cl unit
  timer clear 9
  timer on 9
  cap noi wooldid `yv' `iv' t ttre `wt' `ifq' , att cluster(`cl') permute(shuffle) permreps(`reps') rseed(777) `opts'
  local rc = _rc
  timer off 9
  qui timer list 9
  dis "PCHK|`tag'|RC=" `rc' "|SEC=" %6.1f r(t9)
  if `rc' != 0 exit
  cap confirm matrix e(perm_comparator_b)
  if _rc != 0 {
    dis "PCHK|`tag'|NOCOMPARATOR"
    exit
  }
  mata: wdpcmp()
  dis "PCHK|`tag'|NCMP=" scalar(wdp_ncmp) "|MAXRELDIFF=" %21x scalar(wdp_maxrd) "|AT=`wdp_at'"
  dis "PCHK|`tag'|PDIST=" scalar(wdp_dr) "x" scalar(wdp_dc) "|NMISS=" scalar(wdp_dm) "|S1=" %21x scalar(wdp_s1) "|S2=" %21x scalar(wdp_s2)
  dis "PCHK|`tag'|P_ATT=" %21x e(perm_p_att) "|B_ATT=" %21x _b[main_att] "|SD=" %21x e(perm_sd_att)
end
