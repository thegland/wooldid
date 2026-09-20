clear all
set more off
set processors 1
do cmplib.do
cap noi cmpf , f1(bX) f2(bY) tag("P5 X fweight-count leg vs Y physical leg (aw=1)")
dis "RC1=" _rc
cap noi cmpf , f1(bW) f2(bZ) tag("P5 W count no-weight vs Z count aw=1")
dis "RC2=" _rc
cap noi cmpf , f1(bX) f2(bW) tag("P5 X boottstat-count vs W plain-count")
dis "RC3=" _rc
dis "CMP2_DONE"
