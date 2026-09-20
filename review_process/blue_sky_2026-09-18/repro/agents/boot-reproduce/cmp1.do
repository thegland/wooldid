clear all
set more off
set processors 1
do cmplib.do
cap noi cmpf , f1(bA) f2(bB) tag("A count clu=unit vs B reid clu=unit")
dis "RC1=" _rc
cap noi cmpf , f1(bA) f2(bC) tag("A count vs C reid-via-boottstat")
dis "RC2=" _rc
cap noi cmpf , f1(bD) f2(bE) tag("D count clu=region vs E reid clu=region")
dis "RC3=" _rc
cap noi cmpf , f1(bG) f2(bH) tag("G count nocluster vs H reid nocluster")
dis "RC4=" _rc
dis "--- column inventories ---"
foreach f in bA bB bC bD bE bG bH {
  qui use "`f'.dta", clear
  qui ds
  dis "`f': N=" _N " vars=" wordcount("`r(varlist)'")
  dis "`f' vars: " "`r(varlist)'"
}
dis "CMP1_DONE"
