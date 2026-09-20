clear all
local s1 = 0.1
local s2 = 0.2
local s3 = 0.3
local s4 = 0.4
tempname V RR BB
matrix `V' = J(4,4,.)
forval a = 1/4 {
  matrix `V'[`a',`a'] = (`s`a'')^2
}
mat list `V'
dis "OK1"
matrix `RR' = (1,-1,0,0 \ 1,0,-1,0 \ 1,0,0,-1)
matrix `BB' = (1 \ 2 \ 3 \ 4)
dis "OK2"
