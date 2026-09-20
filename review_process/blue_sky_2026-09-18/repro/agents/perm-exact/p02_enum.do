* p02_enum.do -- EXACT ORACLE on wooldid's assignment enumerator itself.
* Calls mata wooldid_permx() directly over many designs and reconciles, for each:
*    M   = the log-factorial possible-assignment count wooldid reports (wdp_npd)
*    C_A = cols() of the enumerated matrix A
*    D_A = number of DISTINCT columns of A
*    S_A = 1 iff every column preserves the per-stratum multiset of the observed onsets
*    O_A = number of columns equal to the observed assignment
* An exact enumeration requires M == C_A == D_A, S_A == 1 and O_A == 1.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"

which wooldid

* wooldid's mata functions are NOT visible to the session after an ordinary call, so the
* three enumeration functions are extracted VERBATIM from the pinned build
*   sed -n '13624,13785p;13790,13840p' build/wooldid.ado > permx_src.mata.do  (md5 e5cdfc392613bb2c1b11a61a077cf7c2)
* and run here unmodified.  Any disagreement found below is a property of the shipped code.
mata: mata set matastrict off
do permx_src.mata.do
mata: mata describe wooldid_permx()

mata:
// independent structural verifier
real rowvector wdchk(string scalar onsstr, string scalar strstr, real scalar exact,
                     real scalar reps)
{
    real colvector ons, str, obs
    real matrix A
    real scalar M, CA, DA, SA, OA, i, j, dup, si, ok
    real rowvector out
    ons = strtoreal(tokens(onsstr))'
    str = strtoreal(tokens(strstr))'
    if (st_frameexists("wdfr")) st_framedrop("wdfr")
    st_framecreate("wdfr")
    st_framecurrent("wdfr")
    st_addobs(rows(ons))
    st_addvar("double", "o")
    st_addvar("double", "s")
    st_store(., "o", ons)
    st_store(., "s", str)
    st_framecurrent("default")
    st_local("wdp_npd", "")
    st_local("wdp_err", "")
    // pname/lname are tempnames in the real caller (matrices that may not exist),
    // never "": mata `&' evaluates both operands, so st_matrix("") would error
    wooldid_permx("wdfr", "o", "s", reps, exact, "WDA", "WDPNOP", "WDPNOL")
    M = strtoreal(st_local("wdp_npd"))
    if (st_local("wdp_err") != "") {
        out = (M, -1, -1, -1, -1)
        return(out)
    }
    A = st_matrix("WDA")
    CA = cols(A)
    // distinct columns
    DA = 0
    for (i = 1; i <= CA; i++) {
        dup = 0
        for (j = 1; j < i; j++) {
            if (all((A[., i] :== A[., j]) :| (A[., i] :>= . :& A[., j] :>= .))) dup = 1
        }
        if (!dup) DA = DA + 1
    }
    // per-stratum multiset preservation, and count of observed columns
    SA = 1
    OA = 0
    for (i = 1; i <= CA; i++) {
        for (si = 1; si <= rows(uniqrows(str)); si++) {
            ok = 1
            if (sort(select(A[., i], str :== uniqrows(str)[si]), 1) !=
                sort(select(ons,     str :== uniqrows(str)[si]), 1)) ok = 0
            if (!ok) SA = 0
        }
        if (all((A[., i] :== ons) :| (A[., i] :>= . :& ons :>= .))) OA = OA + 1
    }
    out = (M, CA, DA, SA, OA)
    return(out)
}
end

capture program drop runcase
program define runcase
  args tag onsstr strstr
  tempname r
  mata: st_matrix("`r'", wdchk("`onsstr'", "`strstr'", 1, 0))
  local M   = `r'[1,1]
  local CA  = `r'[1,2]
  local DA  = `r'[1,3]
  local SA  = `r'[1,4]
  local OA  = `r'[1,5]
  dis "CASE|`tag'|ons=`onsstr'|str=`strstr'|M=`M'|colsA=`CA'|distinct=`DA'|multisetok=`SA'|nobs=`OA'|err=`wdp_err'"
end

* ---------------- cases -------------------------------------------------
runcase "c01_basic"        "4 7 . . . ."            "1 1 1 1 1 1"
runcase "c02_ties"         "4 4 7 7 . ."            "1 1 1 1 1 1"
runcase "c03_singletonstr" "4 7 . . . ."            "1 1 1 1 1 2"
runcase "c04_allmissstr"   "4 7 . . . ."            "1 1 1 2 2 2"
runcase "c05_uneven"       "4 . 7 . . ."            "1 1 2 2 2 2"
runcase "c06_allsingle"    "1 2 3 4 5 ."            "1 1 1 1 1 1"
runcase "c07_onetreated"   "4 . . . . ."            "1 1 1 1 1 1"
runcase "c08_nonevertr"    "4 4 7 7 9 9"            "1 1 1 1 1 1"
runcase "c09_extmiss"      "4 7 .a . . ."           "1 1 1 1 1 1"
runcase "c10_extmiss2"     "4 .a .b . . ."          "1 1 1 1 1 1"
runcase "c11_strmissing"   "4 7 . . . ."            "1 1 . . 2 2"
runcase "c12_strcross"     "4 . 7 . . ."            "2 2 1 1 1 1"
runcase "c13_9clus7560"    "1 1 2 2 3 3 . . ."      "1 1 1 1 1 1 1 1 1"
runcase "c14_10clusOVER"   "1 1 2 2 3 3 . . . ."    "1 1 1 1 1 1 1 1 1 1"
runcase "c15_bigcancel"    "5 5 5 5 5 5 5 5 5 5 5 ." "1 1 1 1 1 1 1 1 1 1 1 1"
runcase "c16_negons"       "-3 -1 . . . ."          "1 1 1 1 1 1"
runcase "c17_fracons"      "2000.5 2001.5 . . . ."  "1 1 1 1 1 1"
runcase "c18_3strata"      "4 . 7 . 9 ."            "1 1 2 2 3 3"
runcase "c19_strall1clus"  "4 7 9 . . ."            "1 2 3 4 5 6"
runcase "c20_manylevels"   "1 2 3 4 5 6 7 ."        "1 1 1 1 1 1 1 1"

dis "ALLDONE"
