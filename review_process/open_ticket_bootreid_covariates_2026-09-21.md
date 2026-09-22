# TICKET — substitute the re-identified i into covariate payloads under the bootstrap re-id leg

**Opened 2026-09-21 on the author's instruction,** left as a reviewer note on card #7 of the message
review page:

> *"if the i variable visibly appears in the covariate list, it should have the new re-id'd i
> substituted for it in the covariate list for that given bootstrap replication. Please also update
> the help file to note that covariates that depend on i (but where the i variable is replaced with
> some transformation of it) may not work with bootstrap in the cases where physical re-id is used
> (specifying those cases). Please then delete this caution."*

**Status: OPEN. The caution is already deleted (that part of the instruction is done). The fix and
the help note are not built.** Anchor on literal text; the file moves.

## What the caution said, and why deleting it was right

> Caution: a covariate specification names the i variable; under resampling, duplicated clusters'
> copies share their original i codes inside fe()/controls() terms, and only the cell and
> aggregation identities are re-identified.

It warned about a real defect instead of fixing it. Worse, **its detector was close to useless.**
The scan looped over the tokens of `fe()`, `interactivecontrols()`, `controls()`,
`coarsecohortcontrols()` and `ccc_absorb()` and asked whether a token *string-equals* the i
variable's name. But `fe()` and `controls()` are declared `varlist fv`, so the realistic i-keyed
specifications, `fe(i.unit)` or `controls(i.coh##c.yr)`, tokenize to something that never equals
the bare name. `coarsecohortcontrols()` is a plain string that Stata has not unabbreviated, so an
abbreviation slips past too. The caution therefore stayed silent in most of the cases it existed for.

## The defect itself

The physical re-identification leg mints a new identifier, `group()` of the draw index and the
original i, and passes it as the **positional** i of each replication's nested `wooldid` call. The
option string is replayed almost verbatim: exactly one rewrite happens to it today, a `subinstr` on
`cluster()`. So a duplicated cluster's two copies share their original i codes wherever i appears
inside a covariate option, while the cell and aggregation identities are re-identified. The
replication is then fitting a model whose fixed effects pool copies that the resampling meant to
keep apart.

## The fix, in four parts

**Part 1: substitute in the shell, not in the leg.** Where the five option locals still exist
separately, emit a *second* complete option string alongside the one already built, in which each of
the five payloads has had the i token replaced by the re-id'd variable, and pass it to the
replication runner as a new argument. In the leg, choose between the two strings on the existing
nesting flag and leave the `cluster()` rewrite exactly as it is, operating on whichever string was
chosen.

This is preferred over mirroring the `cluster()` surgery inside the leg. An in-code comment there
already warns that string surgery on the assembled option string is fragile, and the alternative
needs the five raw locals plumbed down anyway, so it is strictly more work.

**Part 2: match tokens properly.** Do not copy the old equality test. Use boundary-anchored
replacement on each payload so that `i.unit`, `unit#c.t`, `i.unit##c.t` and a bare `unit` all
rewrite, while `unitsize` and `myunit` do not. `coarsecohortcontrols()` needs care because it is
declared string and arrives un-unabbreviated.

**Part 3: the help-file note.** It must name exactly when the physical re-id leg is taken, since
only those runs are affected. The leg is selected under `bootstrap(standard)` when the run requests
`att_it`, `att_i` or `att_itime`, or `boottstat` together with `cluster()`. Confirm that list against
the code before writing, and say what remains unfixed (Part 4).

**Part 4: state what substitution cannot reach.** A covariate the user built from i *before* calling
wooldid, whose values encode i but whose name is not i, is invisible to any scan. So is a variable
merged in on i. Those keep their original codes across duplicated copies and no rewrite can find
them. This is the residual the help note has to be honest about, and it is the reason the note is
part of the fix rather than an optional extra.

## Verification plan

1. A test that builds a run with `fe(i.unit)` under the re-id leg, and asserts the replication's
   fixed effects separate duplicated copies. The cleanest oracle is that two draws of the same
   cluster must not be collinear.
2. A test that the boundary matcher rewrites `i.unit` and leaves `unitsize` alone.
3. The `cluster()` rewrite must be unchanged: pin it.
4. Full battery, expecting movement **only** on configs that take the re-id leg with an i-keyed
   covariate. If any other config moves, the substitution is too greedy. Bless deliberately.

## Note for the author

This changes bootstrap standard errors on affected runs, so it is a number-mover and wants its own
verified commit and a re-bless. It is also the kind of change where the battery's job is to show
exactly which configs move, not to come back identical.
