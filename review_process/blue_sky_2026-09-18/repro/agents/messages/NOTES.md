# messages recon — NOTES

## State
- extract scripts: extract.py, extract2.py, extract3.py, seq.py (in this dir). msgs3.json = full site table.
- 545 message sites. 218 `dis as error`, 130 `Caution:`, 8 `Note:`, 190 error/exit NNN.
- No `assert` at top level, no mata errprintf, no `///` continuations -> line-based extraction is complete.
- syntax lines 555 (worker) and 4572 (validate) verified BYTE-IDENTICAL (md5 f5b09b83055d7a48bcea4695c623b43d).

## Stata canned text appended by `error NNN` (measured, rc.log)
1 -> (nothing)   9 -> (nothing)   100/111/198 -> "invalid syntax"   119 -> "statement out of context"
301 -> "last estimates not found"   322 -> "something that should be true of your estimation results is not"
411 -> "nonpositive values encountered"   459 -> "something that should be true of your data is not"
498 -> (nothing)   602 -> "file already exists"   907 -> long maxvar advice block

## Confirmed leads so far
1. **contreatwithin+contreat() silently sets safety=off GLOBALLY** (wooldid.ado:4829
   `if "contreatwithin"=="contreatwithin" & "contreat"!="" local safety off`), and `safety`
   is exported to the worker by the c_local blanket at :5508. So the i-positive-integer
   check (:4957), the float>2^24 caution (:4974), negative-weights (:5245), lattice
   nonneg (:5290), cfxlattice nonneg (:5319), the t-ttre grid check (:5358), and every
   worker screen keyed on `safety` are ALL disabled. No message says so.
2. **safety() value is validated at :5473, ~500 lines AFTER 20+ checks consult it.**
   Check is `!inlist("safetyinput","on","off")` -> case-sensitive. `safety(Off)`/`safety(OFF)`
   behaves as safety ON until then; if any earlier safety-gated error fires first the user
   is told "Set option safety(off) to disregard this error" when they believe they did.
3. Bare error/exit sites (no own message): 1388 (r459 from amescreen bycoh), 11947 (r322,
   but messages printed inside the loop), 6226 exit 0, 13018 exit 0, and 7 `error 1`
   (Break passthrough, prints nothing - intentional).
4. Duplicate texts across mutually-exclusive branches (NOT double-fire): 839 vs 12151
   (esfixedbaseperiod vs reltime), 12095 vs 12117 (timetrends vs not). Verified aliases in
   reltime :12012 so safety() DOES reach it.

## TODO
- permctl/bootctl message sequence
- option names in messages the user cannot type
- dead/unreachable messages (qui-wrapped)
- empirical caution counts on 3-5 realistic configs

## FINAL STATE (deliverable written)
Output: /private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/messages.md (883 lines, incl. Appendix A = inventory.txt)

Confirmed leads, in the file as L0-L8:
 L0 maxvar preflight (:1371 si8demand) omits lattice() -> r(907) with NO wooldid message,
    only StataCorp canned advice. ux2.log config F. REPRODUCIBLE.
 L1 contreatwithin+contreat sets working safety=off (:4829); :4957 (i positive-integer) and
    :4974 (float caution) still gate on the WORKING value -> negative i dies r(452)
    "__000005: factor variables may not contain negative values". REGRESSION vs the P6
    ruling at :1346-1350 which fixed the identical class at :1355. mech2.log M2d/M2e.
 L2 safety() validated at :5473, consulted from :4618 -> safety(Off) yields
    "Set option safety(off) to disregard this error". mech.log M1a/M1c.
 L3 all validate/permctl/bootctl cautions precede all worker errors; permctl orders its
    low-reps caution BEFORE ~20 refusals while bootctl orders its cautions last.
 L4 :1371 (96w) + :1443 (90w) fire on the standard i=unit/cluster(i)/jointtests run; :1443
    is gated by NEITHER safety setting.
 L5 :839/:12151 delete rows and name nothing.
 L6 error 1 prints nothing (measured).
 L7 ":4778 latticenum" names an internal local.
 L8 syntax accepts fw, worker refuses at :629.

Empirical counts (ux2.log, 150 units x 12 periods): A=1, B=2, C=4+1note, D=3, E=3, F=3.
