# hv-noparent — verification of "Eight options are silent no-ops without their parent"

## STATE: code read complete, Stata probes pending

## Code-read map (pinned build, wooldid.ado 17021 lines / wooldid.sthlp 1664 lines)

### Options that DO have a requires-parent screen (claim does NOT apply)
| option | parent | screen | ado line |
|---|---|---|---|
| lattice() | contreat() | ERROR r198 | 4818 |
| makecfxplots | contreat() | ERROR r198 | 5268 |
| cfxlattice() | makecfxplots | CAUTION | 5307 |
| esrelativeto() | esfixedbaseperiod | CAUTION | 5371 |
| poisexpresults | poisson | ERROR r198 | 5479 |
| contreatwithin | contreat() | ERROR | 5456 |
| altmethod_adversarial | adversarial | ERROR | 5492 |
| pretrendpower() | esfixedbaseperiod / esprelength | ERROR | 10318/10323 |
| reps/ties/bootalpha/bootdf/bootdisplay/boottstat/bootsave | bootstrap() | ERROR r198 | 6209-6228 |
| permreps/permexact/permtstat/permtci/permstrata/permprobs/permsave/permanalyticse/permabsp/permstepdown | permute() | ERROR r198 | 4882-4894 |

### Candidate ORPHAN-SILENT options (no screen found by reading)
| option | parent | only use site(s) | gated by |
|---|---|---|---|
| saveplots() | makeplots / makecfxplots / histogramcohorteffects | 2494,12390,13120,13325,13440 | all inside plot programs |
| cfxplottypes() | makecfxplots | 5192-5197 | `if "makecfxplots" != ""` |
| cfxplotnose | makecfxplots | 2424,2437 | cfx block |
| makecfxplotsbysg | makecfxplots + subgroups() | 2550 ; 995/12294 (length caution only) | cfx block |
| latticeignoreweights | lattice() 2-token form | 1037 | `contreat!=""` & `lattice!=""` & wordcount==2 |
| ccnoabsorb | contreat()/contreatcontrols() | 1319 | locals only set under contreat |
| contreatpoly() | contreat() | 1108-1109, 5074 | `if "contreat" != ""` |
| contreatelasticitytype() | contreat() | 5034-5038 | code COMMENT at ado:5029 admits "no requires-parent screen and is inert" |
| suppresssgprimaryeffects | subgroups() | 2796,3236,3604,4012,4329 | all nested under `if subgroups!=""` |
| replace | saveplots/permsave/bootsave | n/a | inert with no file writer |

### THE KEY DOC CONTRADICTION (found while reading, not seeded)
sthlp:1181, in the "refuse / caution / skip policy" section, states as policy:
  "No accepted option is a silent no-op, and no check downgrades itself"
but sthlp:976 (makecfxplotsbysg) and sthlp:978 (cfxplotnose) BOTH say the option is
"silently inert" / "inert: nothing reads it, and no message is printed".
=> the help contradicts itself, and the remaining candidates above are silent no-ops
   that are not documented ANYWHERE.

## TODO
- [ ] Stata probe 1: orphan battery (rc + e() identity + log message absence + file absence)
- [ ] Stata probe 2: cfxplottypes() validation bypass (garbage accepted w/o makecfxplots)
- [ ] Stata probe 3 (independent): does saveplots() write what the help promises, and do
      the saved _ES.dta numbers match e()?

---
## STATE 2: ALL PROBES RUN. FINDINGS.md is complete and final.

Three findings written up in FINDINGS.md, each with a self-contained repro that a skeptic can
run on its own:
  F1  repro_noparent.do          nine silent orphan no-ops vs sthlp:1181's "No accepted option
                                 is a silent no-op"                            medium / certain
  F2  repro_cet_mislabel.do      contreatelasticitytype() without contreat() relabels the
                                 results table as elasticities; numbers unchanged.
                                 ado:5029's comment says it "is inert" - false.  medium-high / certain
  F3  repro_semielas_domain.do   semielasticity / eyex-eydx has no domain screen; on an outcome
                                 crossing zero it returns a full inference row with no caution
                                 (att 1.328), or r(459), or an opaque Mata r(3598)
                                 "variable __00001O not found"                  medium / probable

REFUTED from the seeded list: contreatcontrols() IS screened. lattice()/makecfxplots without
contreat() ARE screened. makecfxplotsbysg and cfxplotnose ARE documented as silently inert
(sthlp:976/978), so only the behaviour half of the recon claim holds for those two.

FALSE ALARM I caught myself: latticeignoreweights looked inert even WITH its parent on two
panels; p7_latticewt_controls.do proves it is exactly correct on a shape where the lattice is
not collinear with the unit-by-period cell dummies. Not filed.

Working probe files, in order: p1_orphans.do (battery), p2 (contreat family + label),
p3 (semielasticity crash isolation), p4 (contreat orphans + minimal semielas repro),
p5 (cfx orphans + lattice weight), p6 (lattice weight mechanism), p7 (lattice controls).
Every .log and every E*/L* dump is kept beside them.
