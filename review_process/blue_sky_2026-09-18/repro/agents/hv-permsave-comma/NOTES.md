# hv-permsave-comma — running notes

## Assignment
Verify from scratch: "permsave() lacks the comma check bootsave() got; a Stata-idiomatic
call writes a bizarre filename." Then hunt other perm*/boot* twin asymmetries.

## Build assertion
`which wooldid` in p1_comma.log:59-62 →
/private/.../scratchpad/bluesky/build/wooldid.ado, `*! version 1.4 2026-08-19`. Pinned build, NOT v1.3.

## STATUS

### F1 — CLAIM CONFIRMED, mechanism exactly as stated. (p1_comma.do / p1_comma.log)
- T1  `permsave(draws, replace)`   → **rc 0**, prints
      `permsave(): the permutation draws were written to draws, replace.dta (20 replication rows...)`
      and a file literally named `draws, replace.dta` (56213 bytes, 21 obs, real draws) lands in cwd.
- T1b same command a 2nd time → **rc 602**:
      `Error: permsave(draws, replace) - draws, replace.dta already exists in the working
       directory; add the replace option to overwrite it.`
      → the user is told to add the option they already typed. This is the UX punchline.
- T2  `bootsave(bd2, replace)`     → **rc 198** with the helpful RF-6 message
      `Error: bootsave() takes a plain file name and no comma; replace is a separate option
       of wooldid, so write bootsave(bd2  replace) replace.`   No file written.
- T6  `use draws` → **r(601) file draws.dta not found** (the obvious next step fails).
- replace semantics are otherwise IDENTICAL on both twins:
  T3a perm first=0, T3b perm second no-replace=602, T3c perm+replace=0;
  T4a boot first=0, T4b boot second no-replace=602.
- ado:5866 comment claims "permsave(): the name rules are bootsave()'s" — now FALSE.
  ado:6284 dates the bootsave comma check to RF-6, 2026-09-13. RF-6 shipped on ONE twin.
- sthlp: `grep -i comma wooldid.sthlp` → NO hit for either option. The help documents
  neither refusal, so this is UX/consistency, not doc-vs-code. sthlp:495 does state the
  two options' name rules in parallel ("bootsave() strips a trailing .dta in the same way").

## NEXT / IN FLIGHT
- p2: other twin asymmetries.
  * dots(-5) under permute() vs bootstrap()  [ado:~6270 comment admits it]
  * permsave file layout vs bootsave file layout (permtstat vs boottstat, subgroups())
  * **LEAD (self-found, from reading)**: D-10 cautions exist ONLY on permsave.
    wooldid_permsave (ado:8877-8880) prints two cautions when the onset_ block is dropped
    or short. wooldid_bootsave (ado:6523-6571) prints NONE, though sthlp:495 says bootsave
    "skips its own w_<clustervalue> columns on the same kind of name collision".
    If bootsave drops w_ columns silently that is doc-adjacent + the same one-twin pattern.
  * permsave + permonsetcohorts (§7, one day old).

## FINAL STATE (all probes complete)
- F1 CONFIRMED (seeded) — p1_comma.do, p4_dropcol.do §E
- F2 FOUND (self) — f2_jointwobs.do : boot_sgW/boot_jointWobs/boot_sgWobs absent under boottstat
- F3 FOUND (self) — f3_dots.do : dots(-5) silent under permute(), r(198) under bootstrap()
- 5 negative results recorded in FINDINGS.md (frame-local e(); rngstate replay both engines;
  permonsetcohorts+permsave; layout parity; forced w_ collision did not fire)
- Files: p1_comma.do/.log p2_twins.do/.log p3_replay.do/.log p4_dropcol.do/.log
  f2_jointwobs.do/.log f3_dots.do/.log FINDINGS.md
