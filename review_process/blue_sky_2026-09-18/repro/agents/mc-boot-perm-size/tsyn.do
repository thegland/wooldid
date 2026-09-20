clear all
set more off
cap prog drop p1
prog def p1
  syntax , G(integer 10) T(integer 8) SEED(integer 1) [ HETVAR(real 1) ]
  dis "OK1 g=`g' t=`t' seed=`seed' hetvar=`hetvar'"
end
cap noi p1 , g(10) t(8) seed(5)
dis "RC1=" _rc
cap prog drop p2
prog def p2
  syntax , NGRP(integer 10) NPER(integer 8) SEED(integer 1) [ HETVAR(real 1) ]
  dis "OK2 ngrp=`ngrp' nper=`nper' seed=`seed' hetvar=`hetvar'"
end
cap noi p2 , ngrp(10) nper(8) seed(5)
dis "RC2=" _rc
