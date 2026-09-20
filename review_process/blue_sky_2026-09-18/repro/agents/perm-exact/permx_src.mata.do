cap mata: mata drop wooldid_permx()
cap mata: mata drop wooldid_permx_enum()
cap mata: mata drop wooldid_permx_msp()
mata:
void wooldid_permx(string scalar fr, string scalar onsv, string scalar strv,
                   real scalar reps, real scalar exact, string scalar assignname,
                   string scalar pname, string scalar lname)
{
    real matrix A
    real colvector ons, str, idx, levels, v
    real scalar C, si, lnM, M, R, r, k, capexact, capwarn, oc
    string scalar nps
    st_local("wdp_err", "")
    st_local("wdp_caut1", "")
    st_local("wdp_caut2", "")
    // BS-23 (author ruling 2026-09-07): the accumulator ROW that holds the observed
    // assignment. 0 = "not an enumeration" - sampled draws, and the permute(probabilities)
    // generator, which returns before the exact branch below.
    st_local("wdp_obsrow", "0")
    string scalar curfr
    curfr = st_framecurrent()
    st_framecurrent(fr)
    ons = st_data(., onsv)
    str = st_data(., strv)
    st_framecurrent(curfr)
    C = rows(ons)
    st_local("wdp_caut0", "")
    if (pname != "" & rows(st_matrix(pname)) > 0) {
        wooldid_permxp(ons, st_matrix(pname), st_matrix(lname), reps, assignname)
        return
    }
    // possible-assignment count: product over strata of multiset-permutation counts
    lnM = 0
    for (si = 1; si <= rows(uniqrows(str)); si++) {
        v = select(ons, str :== uniqrows(str)[si])
        lnM = lnM + lnfactorial(rows(v))
        levels = uniqrows(v)
        for (k = 1; k <= rows(levels); k++) {
            lnM = lnM - lnfactorial(sum(v :== levels[k] :| (v :>= . :& levels[k] :>= .)))
        }
    }
    M = (lnM < 700 ? round(exp(lnM)) : .)
    if (lnM < 34) nps = strofreal(round(exp(lnM)), "%20.0f")
    else nps = strofreal(exp(lnM - floor(lnM/ln(10))*ln(10)), "%9.4f") + "e+" + strofreal(floor(lnM/ln(10)), "%3.0f")
    st_local("wdp_npd", strofreal(M, "%20.0f"))
    st_local("wdp_nps", strtrim(nps))
    // TWO OF THESE FIVE NUMBERS ARE ARITHMETIC, THREE ARE POLICY. Arithmetic: 700 is the
    // double-overflow guard (ln of the largest double is about 709.78), which is what makes
    // the count deliberately MISSING rather than +inf; 34 is the exact-integer boundary
    // (exp(34) is about 5.8e14, under 2^53), so below it round() and the %20.0f print are
    // lossless and above it the code must switch to the mantissa/exponent string. Policy,
    // and tunable: capexact is the permexact ceiling - it bounds the C x M enumeration
    // matrix and the length of the rep loop, not just the user's patience - capwarn the
    // slow-run warning, and the 0.1*M below the sampling-overlap caution
    capexact = 10000
    capwarn = 1000
    if (exact) {
        if (M >= . | M > capexact) {
            st_local("wdp_err", "Error: permexact requires enumerating every assignment, but there are " + strtrim(nps) + " (ceiling: " + strofreal(capexact) + "). Use permreps().")
            return
        }
        if (M > capwarn) st_local("wdp_caut1", "permexact is enumerating " + strtrim(nps) + " assignments (over " + strofreal(capwarn) + "); this may be slow.")
        R = M
    }
    else {
        R = reps
        if (M < . & reps > 0.1*M) st_local("wdp_caut2", "permreps(" + strofreal(reps) + ") exceeds 10% of the " + strtrim(nps) + " possible assignments; consider permexact.")
    }
    st_local("wdp_R", strofreal(R))
    // NOTE the possible-assignment count and the enumeration are produced by two
    // INDEPENDENT algorithms: R comes from the log-factorial formula above, the columns of
    // A from the recursive enumerator below (this J() allocation is thrown away in the
    // exact branch). Only their AGREEMENT makes the rep loop correct, and nothing
    // reconciles them. The loop is driven by wdp_R and reads the matrix with el(), which
    // returns MISSING out of range rather than erroring (measured, Stata 19), and a missing
    // onset is wooldid's never-treated code - so R > cols(A) yields replications in which
    // every cluster is silently a control, and R < cols(A) yields a non-exhaustive
    // enumeration whose completeness flag still claims exact shares. The two routines agree
    // because both bucket every missing onset as ONE never-treated level (the
    // `:| (v :>= . :& levels[k] :>= .)' clause in the counter, the same test in
    // wooldid_permx_msp): exact for `.', wrong for extended missings (.a, .b), which
    // uniqrows() keeps as distinct levels
    A = J(C, R, .)
    if (exact) {
        // enumerate per-stratum multiset permutations, then take the cartesian product
        A = wooldid_permx_enum(ons, str)
        // BS-23 (author ruling 2026-09-07): locate the observed assignment among the
        // enumerated columns. It is always one of them (the enumeration is every distinct
        // multiset permutation), but NOT necessarily the first: each stratum is enumerated
        // from its SORTED onset vector, so the identity column sits wherever the observed
        // order falls. Missing onsets (never-treated) are one level here, exactly as in
        // wooldid_permx_msp's own multiset test.
        for (oc = 1; oc <= cols(A); oc++) {
            if (all((A[., oc] :== ons) :| (A[., oc] :>= . :& ons :>= .))) {
                st_local("wdp_obsrow", strofreal(oc))
                oc = cols(A)
            }
        }
    }
    else {
        for (r = 1; r <= R; r++) {
            v = ons
            for (si = 1; si <= rows(uniqrows(str)); si++) {
                idx = selectindex(str :== uniqrows(str)[si])
                v[idx] = jumble(ons[idx])
            }
            A[., r] = v
        }
    }
    st_matrix(assignname, A)
}

real matrix wooldid_permx_enum(real colvector ons, real colvector str)
{
    real matrix A, B, S1
    real colvector idx
    real scalar si, n1, n2, a, b
    // per-stratum enumerations, combined by cartesian product
    A = J(rows(ons), 0, .)
    for (si = 1; si <= rows(uniqrows(str)); si++) {
        idx = selectindex(str :== uniqrows(str)[si])
        S1 = wooldid_permx_msp(sort(ons[idx], 1))
        if (cols(A) == 0) {
            A = J(rows(ons), cols(S1), .)
            for (b = 1; b <= cols(S1); b++) A[idx, b] = S1[., b]
        }
        else {
            n1 = cols(A)
            n2 = cols(S1)
            B = J(rows(ons), n1*n2, .)
            for (a = 1; a <= n1; a++) {
                for (b = 1; b <= n2; b++) {
                    B[., (a-1)*n2 + b] = A[., a]
                    B[idx, (a-1)*n2 + b] = S1[., b]
                }
            }
            A = B
        }
    }
    return(A)
}

real matrix wooldid_permx_msp(real colvector v)
{
    real matrix R, S
    real colvector rest, lv
    real scalar k, j
    // all distinct permutations of multiset v (v sorted), columns = permutations
    if (rows(v) == 1) return(v)
    lv = uniqrows(v)
    R = J(rows(v), 0, .)
    for (k = 1; k <= rows(lv); k++) {
        j = selectindex(v :== lv[k] :| (v :>= . :& lv[k] :>= .))[1]
        if (j == rows(v)) rest = v[|1 \ j-1|]
        else if (j == 1) rest = v[|2 \ .|]
        else rest = v[|1 \ j-1|] \ v[|j+1 \ .|]
        S = wooldid_permx_msp(rest)
        R = R, (J(1, cols(S), lv[k]) \ S)
    }
    return(R)
}
end
cap mata: mata drop wooldid_permxp()
mata:
void wooldid_permxp(real colvector ons, real matrix P, real rowvector lvl,
                    real scalar reps, string scalar assignname)
{
    real matrix A
    real colvector cum
    real scalar C, K, c, r, k, u, lnS, zp, obsk
    string scalar nps
    C = rows(ons)
    K = cols(lvl)
    lnS = 0
    zp = 0
    for (c = 1; c <= C; c++) {
        lnS = lnS + ln(sum(P[c, .] :> 0))
        obsk = 0
        for (k = 1; k <= K; k++) {
            if (ons[c] == lvl[k]) obsk = k
        }
        if (obsk == 0 | (obsk > 0 ? P[c, obsk] <= 0 : 0)) zp = 1
    }
    if (lnS < 34) nps = strofreal(round(exp(lnS)), "%20.0f")
    else nps = strofreal(exp(lnS - floor(lnS/ln(10))*ln(10)), "%9.4f") + "e+" + strofreal(floor(lnS/ln(10)), "%3.0f")
    st_local("wdp_npd", strofreal((lnS < 700 ? round(exp(lnS)) : .), "%20.0f"))
    st_local("wdp_nps", strtrim(nps))
    st_local("wdp_R", strofreal(reps))
    if (exp(lnS) < . & reps > 0.1*exp(lnS)) st_local("wdp_caut2", "permreps(" + strofreal(reps) + ") exceeds 10% of the " + strtrim(nps) + " possible assignments in the design's support.")
    if (zp) st_local("wdp_caut0", "the OBSERVED assignment has zero probability under the specified design, so the reported p-value is a tail probability under that design distribution, not a randomization test that includes the observed assignment.")
    A = J(C, reps, .)
    for (r = 1; r <= reps; r++) {
        for (c = 1; c <= C; c++) {
            // this loop can leave A[c, r] MISSING if u exceeds the last cumulative bin,
            // and a missing onset means never-treated downstream - so the cluster silently
            // becomes a control in that replication, rc 0, no message. It cannot happen
            // because wooldid_permp (NFPERM-06) rescales every row so its cumulative sum
            // ends at 1 to within about K*eps, which the 1e-12 cushion below covers; a P
            // matrix arriving by any other route must be normalized first. `k = K + 1' is
            // the deliberate loop break, written as a counter assignment
            cum = runningsum(P[c, .])'
            u = runiform(1, 1)
            for (k = 1; k <= K; k++) {
                if (u <= cum[k] + 1e-12) {
                    A[c, r] = lvl[k]
                    k = K + 1
                }
            }
        }
    }
    st_matrix(assignname, A)
}
end
