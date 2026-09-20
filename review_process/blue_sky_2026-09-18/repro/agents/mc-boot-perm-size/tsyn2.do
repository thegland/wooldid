clear all
set more off
cap prog drop a1
prog def a1
  syntax , NGRP(integer 10)
  dis "A1 ok `ngrp'"
end
cap noi a1 , ngrp(7)
dis "RCA1=" _rc
cap prog drop a2
prog def a2
  syntax [, NGRP(integer 10) HETV(real 1)]
  dis "A2 ok `ngrp' `hetv'"
end
cap noi a2 , ngrp(7)
dis "RCA2=" _rc
cap prog drop a3
prog def a3
  syntax [, NGRP(integer 10) NPER(integer 8) SEEDX(integer 1) HETV(real 1)]
  dis "A3 ok `ngrp' `nper' `seedx' `hetv'"
end
cap noi a3 , ngrp(7) nper(9) seedx(3)
dis "RCA3=" _rc
