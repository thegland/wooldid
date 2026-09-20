clear all
set more off
set processors 1
do cmplib.do
cap noi cmpf , f1(s1) f2(s2) tag("p7a1 COUNT sorted vs shuffled")
dis "RC1=" _rc
cap noi cmpf , f1(s3) f2(s4) tag("p7a2 REID sorted vs shuffled")
dis "RC2=" _rc
cap noi cmpf , f1(s1) f2(s3) tag("p7a3 sorted COUNT vs REID")
dis "RC3=" _rc
cap noi cmpf , f1(s2) f2(s4) tag("p7a4 shuffled COUNT vs REID")
dis "RC4=" _rc
cap noi cmpf , f1(n20) f2(n40) tag("p7b reps20 vs first20 of reps40") maxrep(20)
dis "RC5=" _rc
cap noi cmpf , f1(k1) f2(k2) tag("p7c keep-i COUNT vs REID")
dis "RC6=" _rc
cap noi cmpf , f1(bX) f2(bY) tag("p5 fweight-count vs physical(aw=1)")
dis "RC7=" _rc
cap noi cmpf , f1(bW) f2(bZ) tag("p5 count no-wt vs count aw=1")
dis "RC8=" _rc
cap noi cmpf , f1(bX) f2(bW) tag("p5 boottstat-count vs plain-count")
dis "RC9=" _rc
cap noi cmpf , f1(c1) f2(c2) tag("p12 collapse vs nofastcollapse")
dis "RC10=" _rc
dis "CMP3_DONE"
