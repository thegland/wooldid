* bandcmp.do -- compare e(esband_<fam>) with the pointwise interval in
* e(wdidmainresults). Requires locals: BFAM (family, e.g. att), BTAG.
mat __MRB = e(wdidmainresults)
mat __SBB = e(esband_`BFAM')
mata:
  MR = st_matrix("__MRB")
  rs = st_matrixrowstripe("__MRB")
  SB = st_matrix("__SBB")
  cs = st_matrixcolstripe("__SBB")
  fam = st_local("BFAM")
  tag = st_local("BTAG")
  nviol = 0
  for (j = 1; j <= cols(SB); j++) {
    want = "ES_" + fam
    hit = 0
    for (r = 1; r <= rows(rs); r++) {
      if (rs[r,1] == want & rs[r,2] == cs[j,2]) hit = r
    }
    if (hit == 0) {
      printf("%s BANDROW_NOTFOUND %s %s\n", tag, want, cs[j,2])
      continue
    }
    plb = MR[hit,5]; pub = MR[hit,6]
    blb = SB[1,j];   bub = SB[2,j]
    vlo = (blb > plb + 1e-12)
    vhi = (bub < pub - 1e-12)
    if (vlo | vhi) nviol++
    printf("%s h=%s  point=[%12.8f,%12.8f]  band=[%12.8f,%12.8f]  narrow_lo=%f narrow_hi=%f  halfratio=%9.5f\n",
           tag, cs[j,2], plb, pub, blb, bub, vlo, vhi, (bub-blb)/(pub-plb))
  }
  printf("%s BANDVIOLATIONS=%f\n", tag, nviol)
end
