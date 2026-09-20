* richer comparator dump: prints the 4 worst columns
mata:
void wdpcmp2(string scalar tag)
{
    real matrix B, C
    string matrix bs, cs
    real scalar k, j, n
    real colvector d
    string colvector nm
    B  = st_matrix("e(b)")
    bs = st_matrixcolstripe("e(b)")
    C  = st_matrix("e(perm_comparator_b)")
    cs = st_matrixcolstripe("e(perm_comparator_b)")
    n = 0
    d = J(0,1,.)
    nm = J(0,1,"")
    for (k=1; k<=cols(C); k++) {
        for (j=1; j<=cols(B); j++) {
            if (bs[j,1]==cs[k,1] & bs[j,2]==cs[k,2]) {
                if (C[1,k]<. & B[1,j]<.) {
                    real scalar dd
                    dd = abs(C[1,k]-B[1,j])
                    if (abs(B[1,j])>1e-12) dd = dd/abs(B[1,j])
                    n = n + 1
                    d  = d  \ dd
                    nm = nm \ (cs[k,1]+":"+cs[k,2])
                }
                j = cols(B) + 1
            }
        }
    }
    printf("CMP2|%s|NB=%f|NC=%f|NMATCH=%f\n", tag, cols(B), cols(C), n)
    if (n>0) {
        real colvector o
        real scalar q
        o = order(-d, 1)
        for (q=1; q<=min((4,n)); q++) {
            printf("CMP2|%s|W%f|%s|reldiff=%21x\n", tag, q, nm[o[q]], d[o[q]])
        }
    }
}
end
cap prog drop pchk2
prog def pchk2
  syntax anything(name=tag) [, OPTS(string asis) WT(string asis) IFQ(string asis) YV(string) REPS(integer 30) IV(string) CL(string)]
  if "`yv'" == "" local yv y
  if "`iv'" == "" local iv unit
  if "`cl'" == "" local cl unit
  cap noi wooldid `yv' `iv' t ttre `wt' `ifq' , att cluster(`cl') permute(shuffle) permreps(`reps') rseed(777) `opts'
  local rc = _rc
  dis "CMP2|`tag'|RC=" `rc'
  if `rc' != 0 exit
  cap confirm matrix e(perm_comparator_b)
  if _rc != 0 {
    dis "CMP2|`tag'|NOCOMPARATOR"
    exit
  }
  mata: wdpcmp2("`tag'")
end
