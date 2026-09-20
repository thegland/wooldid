mata:
void gtoracle(string scalar tag, string scalar gtname, real scalar lvl)
{
    real matrix GT
    real matrix B
    real matrix V
    real matrix cells
    real matrix srt
    string matrix bs
    string matrix cs
    real scalar j, nb, n, i, k, df, g, est, se, pv, crit, lo, hi, al
    real scalar pcoh, pry, pest, pse, plb, pub, ppv
    real scalar mb, mse, mp, mlb, mub, dmax
    real scalar rb, rse, rp, rlb, rub
    B  = st_matrix("e(b)")
    V  = st_matrix("e(V)")
    df = st_numscalar("e(df_r)")
    bs = st_matrixcolstripe("e(b)")[.,2]
    nb = rows(bs)
    cells = J(0,3,.)
    for (j=1; j<=nb; j++) {
        if (regexm(bs[j], "^([0-9]+)\.([^#]+)#([0-9]+)\.([^#]+)#1\.([^#]+)$")) {
            cells = cells \ (strtoreal(regexs(1)), strtoreal(regexs(3)), j)
        }
    }
    n = rows(cells)
    GT = st_matrix(gtname)
    cs = st_matrixcolstripe(gtname)[.,2]
    pcoh = 0
    pry  = 0
    pest = 0
    pse  = 0
    plb  = 0
    pub  = 0
    ppv  = 0
    for (k=1; k<=rows(cs); k++) {
        if (cs[k]=="i")     pcoh = k
        if (cs[k]=="relyr") pry  = k
        if (cs[k]=="estimate") pest = k
        if (cs[k]=="se")    pse  = k
        if (cs[k]=="p")     ppv  = k
        if (substr(cs[k],1,3)=="lb_") plb = k
        if (substr(cs[k],1,3)=="ub_") pub = k
    }
    printf("OR|%s|NCOEF|%f|NROW|%f|df_r|%21x\n", tag, n, rows(GT), df)
    if (n != rows(GT)) printf("OR|%s|COUNTMISMATCH|coef=%f table=%f\n", tag, n, rows(GT))
    srt = sort(cells, (1,2))
    mb  = 0
    mse = 0
    mp  = 0
    mlb = 0
    mub = 0
    dmax = min((n, rows(GT)))
    al = (100-lvl)/200
    for (i=1; i<=dmax; i++) {
        g  = srt[i,1]
        j  = srt[i,3]
        est = B[1,j]
        se = .
        if (V[j,j] < . & V[j,j] > 0) se = sqrt(V[j,j])
        pv = .
        lo = .
        hi = .
        if (se < .) {
            if (df < .) {
                pv = 2*ttail(df, abs(est/se))
                crit = invttail(df, al)
            }
            else {
                pv = 2*normal(-abs(est/se))
                crit = invnormal(1-al)
            }
            lo = est - crit*se
            hi = est + crit*se
        }
        if (GT[i,pcoh] != g) printf("OR|%s|KEYMISMATCH|row=%f table_i=%f coef_i=%f\n", tag, i, GT[i,pcoh], g)
        rb = reldif(GT[i,pest], est)
        rse = 0
        if (!(GT[i,pse] >= . & se >= .)) rse = reldif(GT[i,pse], se)
        rp = 0
        if (!(GT[i,ppv] >= . & pv >= .)) rp = reldif(GT[i,ppv], pv)
        rlb = 0
        if (!(GT[i,plb] >= . & lo >= .)) rlb = reldif(GT[i,plb], lo)
        rub = 0
        if (!(GT[i,pub] >= . & hi >= .)) rub = reldif(GT[i,pub], hi)
        if (rb  > mb)  mb  = rb
        if (rse > mse) mse = rse
        if (rp  > mp)  mp  = rp
        if (rlb > mlb) mlb = rlb
        if (rub > mub) mub = rub
    }
    printf("OR|%s|MAXRELDIF|b=%10.3e se=%10.3e p=%10.3e lb=%10.3e ub=%10.3e\n", tag, mb, mse, mp, mlb, mub)
}
mata mosave gtoracle(), dir("/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects") replace
end
