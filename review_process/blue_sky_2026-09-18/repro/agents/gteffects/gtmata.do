mata:
mata clear
void gtident(string scalar tag, real scalar rtol)
{
    real matrix GT, B
    string matrix cs, bs
    real scalar nr, i, k, pon, pct, pry, pest, pw, bcol
    real scalar num, den, nmissw, npost, rec, tgt, rd, njm, nbad
    string vector ests
    string scalar v

    GT = st_matrix("e(gteffects)")
    B  = st_matrix("e(b)")
    cs = st_matrixcolstripe("e(gteffects)")[.,2]
    bs = st_matrixcolstripe("e(b)")[.,2]
    nr = rows(GT)
    printf("GT|%s|ROWS|%f\n", tag, nr)
    pon=0; pct=0; pry=0; pest=0
    for (k=1;k<=rows(cs);k++) {
        if (cs[k]=="onset") pon=k
        if (cs[k]=="caltime") pct=k
        if (cs[k]=="relyr") pry=k
        if (cs[k]=="estimate") pest=k
    }
    njm=0; nbad=0
    for (i=1;i<=nr;i++) {
        if (GT[i,pon]>=.) njm++
        else if (abs((GT[i,pct]-GT[i,pon]) - GT[i,pry])>1e-9) nbad++
    }
    printf("GT|%s|JOINMISS|%f\n", tag, njm)
    printf("GT|%s|CALTIME_MINUS_ONSET_NE_RELYR|%f\n", tag, nbad)
    ests = ("att","att_it","att_i","att_itime")
    for (k=1;k<=4;k++) {
        v = "w_" + ests[k]
        pw = 0
        for (i=1;i<=rows(cs);i++) if (cs[i]==v) pw=i
        bcol = 0
        for (i=1;i<=cols(B);i++) if (bs[i]=="main_"+ests[k]) bcol=i
        if (bcol==0) {
            printf("GT|%s|%s|SKIP|not-in-e(b)\n", tag, ests[k])
            continue
        }
        if (pw==0) {
            printf("GT|%s|%s|SKIP|no-weight-col\n", tag, ests[k])
            continue
        }
        num=0; den=0; nmissw=0; npost=0
        for (i=1;i<=nr;i++) {
            if (GT[i,pry]>=0 & GT[i,pry]<.) {
                npost++
                if (GT[i,pw]>=. | GT[i,pest]>=.) {
                    nmissw++
                    continue
                }
                num = num + GT[i,pw]*GT[i,pest]
                den = den + GT[i,pw]
            }
        }
        if (den==0) {
            printf("GT|%s|%s|ZERODEN|npost=%f nmissw=%f\n", tag, ests[k], npost, nmissw)
            continue
        }
        rec = num/den
        tgt = B[1,bcol]
        rd  = abs(rec-tgt)/max((abs(tgt),1e-12))
        printf("GT|%s|%s|%s|reldiff=%10.3e rebuilt=%21x reported=%21x npost=%f nmissw=%f sumw=%12.6g\n",
               tag, ests[k], (rd>rtol ? "FAIL" : "ok"), rd, rec, tgt, npost, nmissw, den)
    }
}

void gtidentpre(string scalar tag, real scalar rtol)
{
    real matrix GT, B
    string matrix cs, bs
    real scalar nr,i,k,pry,pest,pw,bcol,num,den,rec,tgt,rd,npre
    string vector ests
    GT = st_matrix("e(gteffects)")
    B  = st_matrix("e(b)")
    cs = st_matrixcolstripe("e(gteffects)")[.,2]
    bs = st_matrixcolstripe("e(b)")[.,2]
    nr = rows(GT)
    pry=0; pest=0
    for (k=1;k<=rows(cs);k++) {
        if (cs[k]=="relyr") pry=k
        if (cs[k]=="estimate") pest=k
    }
    ests = ("att","att_it","att_i","att_itime")
    for (k=1;k<=4;k++) {
        pw=0
        for (i=1;i<=rows(cs);i++) if (cs[i]=="w_"+ests[k]) pw=i
        bcol=0
        for (i=1;i<=cols(B);i++) if (bs[i]=="pretreatment_"+ests[k]) bcol=i
        if (bcol==0 | pw==0) continue
        num=0; den=0; npre=0
        for (i=1;i<=nr;i++) {
            if (GT[i,pry]<0) {
                npre++
                if (GT[i,pw]>=. | GT[i,pest]>=.) continue
                num = num + GT[i,pw]*GT[i,pest]
                den = den + GT[i,pw]
            }
        }
        if (den==0) {
            printf("GT|%s|PRE_%s|ZERODEN|npre=%f\n", tag, ests[k], npre)
            continue
        }
        rec=num/den
        tgt=B[1,bcol]
        rd=abs(rec-tgt)/max((abs(tgt),1e-12))
        printf("GT|%s|PRE_%s|%s|reldiff=%10.3e rebuilt=%21x reported=%21x npre=%f\n",
               tag, ests[k], (rd>rtol ? "FAIL":"ok"), rd, rec, tgt, npre)
    }
}

void gtidentes(string scalar tag, real scalar rtol)
{
    real matrix GT, B
    string matrix cs, bs
    real scalar nr,i,k,pry,pest,pw,num,den,rec,tgt,rd,q,cand,ry,bcol
    real colvector rys
    string scalar nm
    string vector ests
    GT = st_matrix("e(gteffects)")
    B  = st_matrix("e(b)")
    if (cols(B)==0) return
    cs = st_matrixcolstripe("e(gteffects)")[.,2]
    bs = st_matrixcolstripe("e(b)")[.,2]
    nr = rows(GT)
    pry=0; pest=0
    for (k=1;k<=rows(cs);k++) {
        if (cs[k]=="relyr") pry=k
        if (cs[k]=="estimate") pest=k
    }
    ests = ("att","att_it","att_i","att_itime")
    for (k=1;k<=4;k++) {
        pw=0
        for (i=1;i<=rows(cs);i++) if (cs[i]=="w_"+ests[k]) pw=i
        if (pw==0) continue
        rys = uniqrows(GT[.,pry])
        for (q=1; q<=rows(rys); q++) {
            ry = rys[q]
            for (cand=1; cand<=2; cand++) {
                nm = (cand==1 ? "ES_"+ests[k]+"_post"+strofreal(ry) : "ES_"+ests[k]+"_pre"+strofreal(-ry))
                bcol=0
                for (i=1;i<=cols(B);i++) if (bs[i]==nm) bcol=i
                if (bcol==0) continue
                num=0; den=0
                for (i=1;i<=nr;i++) {
                    if (GT[i,pry]==ry & GT[i,pw]<. & GT[i,pest]<.) {
                        num = num + GT[i,pw]*GT[i,pest]
                        den = den + GT[i,pw]
                    }
                }
                if (den==0) {
                    printf("GT|%s|ES_%s|ZERODEN\n", tag, nm)
                    continue
                }
                rec=num/den
                tgt=B[1,bcol]
                rd=abs(rec-tgt)/max((abs(tgt),1e-12))
                printf("GT|%s|ES_%s|%s|reldiff=%10.3e rebuilt=%21x reported=%21x\n",
                       tag, nm, (rd>rtol?"FAIL":"ok"), rd, rec, tgt)
            }
        }
    }
}
mata mosave gtident(), dir("/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects") replace
mata mosave gtidentpre(), dir("/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects") replace
mata mosave gtidentes(), dir("/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects") replace
end
