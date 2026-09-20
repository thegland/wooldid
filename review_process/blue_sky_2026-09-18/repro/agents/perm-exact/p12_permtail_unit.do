* p12_permtail_unit.do -- SELF-CONTAINED Mata unit test of wooldid_permtail's share branch.
*
* THE CONTRADICTION.  wooldid_permtail (ado:15992-16000) takes the PLAIN share whenever
*     if (exact | memb) { plo = m0/n ; phi = m1/n }   else  add-one
* where `exact' is COMPLETENESS of the enumeration and `memb' is MEMBERSHIP (the observed
* assignment's own replication carries the observed statistic).  But wooldid_permstat
* (ado:16147) sets `wdp_addone' on  !memb  ALONE:
*     if (M < . & !memb & sobs < . & sum(s :< .) >= 2) nonm = nonm + " " + <element name>
* and the caution the caller then prints (ado:8297) says:
*     "... their reference set does not contain the observed statistic and their p-values
*      use the add-one form (1 + count)/(n + 1) rather than an exact share"
* For an element with exact = 1 and memb = 0 that sentence is false: the code uses the
* plain share over a reference set that does NOT contain the observed statistic, and the
* reported p is then 0 whenever the observed statistic sits outside the draws.
* The add-one form exists precisely to floor that at 1/(n+1).
*
* This file extracts the three kernels VERBATIM from the pinned build
*     { echo "mata:"; sed -n '15908,15915p;15937,15943p;15945,16062p' build/wooldid.ado;
*       echo "end"; } > permtail_src.mata.do        (md5 7b878e5b1e4484a0527f855ec6a65154)
* and calls wooldid_permtail directly, so nothing here depends on my notes.
*
* NOT DEMONSTRATED HERE: that wooldid can be driven into exact=1 with memb=0 from the
* command line.  I could not: membership held on every element of every design I built,
* and the two estimation paths whose fits are not bit-reproducible run-to-run
* (unconditionalse, customvce(cluster ...)) are refused by permute().  This is a
* contradiction between a branch and the message that describes it, of unknown reachability.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-exact"
which wooldid
mata: mata set matastrict off
do permtail_src.mata.do

mata:
  // a complete enumeration of 30 assignments, drawn statistics 1.00 .. 1.29
  draws = 1 :+ (0::29)/100
  al = 0.025
  zc = invnormal(1 - al)
  // the observed statistic sits ONE ULP BELOW the smallest draw -- i.e. the identity
  // replication's own value differs from the comparator's in the last bit, which is
  // exactly what memb = 0 records
  sobs = draws[1] - 2*epsilon(draws[1])
  printf("\nobserved statistic below every draw; n = %f\n", rows(draws))
  r1 = wooldid_permtail(draws, sobs, 1, 1, al, zc, 0)   // exact=1 memb=1 (the normal case)
  r0 = wooldid_permtail(draws, sobs, 1, 0, al, zc, 0)   // exact=1 memb=0 (the flagged case)
  ra = wooldid_permtail(draws, sobs, 0, 0, al, zc, 0)   // exact=0 memb=0 (true add-one)
  printf("PERMTAIL exact=1 memb=1 : p = %21x  (%f)\n", r1[3], r1[3])
  printf("PERMTAIL exact=1 memb=0 : p = %21x  (%f)   <-- the element wdp_addone names\n", r0[3], r0[3])
  printf("PERMTAIL exact=0 memb=0 : p = %21x  (%f)   <-- what the caution PROMISES, 2/(n+1)\n", ra[3], ra[3])
  printf("p_lo/p_hi at exact=1 memb=0 : %f  %f\n", r0[14], r0[15])
  printf("p_lo/p_hi at exact=0 memb=0 : %f  %f\n", ra[14], ra[15])
  // same thing under the one-sided and absolute-value conventions
  printf("permside(left)  exact=1 memb=0 : %f   (add-one branch would give %f)\n",
         wooldid_permtail(draws, sobs, 1, 0, al, zc, 3)[3],
         wooldid_permtail(draws, sobs, 0, 0, al, zc, 3)[3])
  printf("permabsp        exact=1 memb=0 : %f   (add-one branch would give %f)\n",
         wooldid_permtail(draws, sobs, 1, 0, al, zc, 1)[3],
         wooldid_permtail(draws, sobs, 0, 0, al, zc, 1)[3])
  printf("Monte-Carlo se of the share at exact=1 memb=0 : %f (missing = treated as exact)\n", r0[16])
end
dis "ALLDONE"
