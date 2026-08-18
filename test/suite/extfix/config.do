* Copied from the colleagues' validation bundle
* /Users/hegland/Documents/Research/wooldid_dev/wooldid_external/wooldid_test_comparison_code/tests/config.do
* on 2026-08-18; pure Stata, deterministic, no RNG (config.do sets seeds but the
* fixture formulas draw nothing). Content unmodified except this header.
version 17.0

* Deterministic suite-wide Stata settings. Case files may override a setting only
* when the override is part of the behavior being tested and is recorded.
set more off
set rmsg off
set trace off
set varabbrev off
set linesize 255
set level 95
set type double
set dp period
set graphics on
set scheme s2color
set emptycells keep
set matsize 800
set maxiter 16000
set seed 8675309
set sortseed 314159265
