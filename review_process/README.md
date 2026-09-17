# review_process

Working records for the wooldid review and verification campaigns. What lives here is the
material that answers **what the software does and why**: the decision record, the live defect
ledger, the operating brief, current tickets, and the standing analyses.

What does **not** live here any more is the week-by-week process material - sweep reports, chunk
reports, implementation plans and ship records for work that has shipped. Those were removed on
2026-09-16 because they had stopped being read: an audit of fourteen documents that day turned up
**191 stale claims**, several of which had already caused real errors. A record nobody maintains
is worse than no record, because its stale status lines get mistaken for current ones - on
2026-09-16 a stale line in one of these files caused thirteen already-ruled decisions to be
reported to the author as still open.

## The files removed on 2026-09-16, and how to read them

Nothing was lost. Git keeps every one of them in full. To read one:

```sh
git log --diff-filter=D --oneline -- review_process/        # find the removing commit
git show <commit>^:review_process/<file>                    # print the file as it last stood
git show <commit>^:review_process/<file> > /tmp/<file>      # or write it out
```

Removed, with what replaced them as the live record:

| removed | where the live answer is now |
|---|---|
| `bug_sweep_2026-09-06.md` | decisions in `rulings_2026-09-07.md`; anything unfixed in `open_defects.md` |
| `consistency_sweep_2026-09-06.md` | same |
| `practice_sweep_2026-09-06.md` | same, plus the help file itself |
| `qr_ship_2026-09-06.md` | `rulings_2026-09-07.md` |
| `quality_review_2026-09-05.md` | `rulings_2026-09-07.md`, `open_defects.md` |
| `bootstrap_plan.md`, `bootstrap_report_2026-08-30.md` | `bootstrap_coverage_2026-09-02.md` (kept) and the `bootstrap()` help entry |
| `f6_permutation_plan.md` | `permute_coverage_2026-09-02.md` (kept) and the `permute()` help entry |
| `si20_ship_2026-09-05.md` | the SI-20 section of `wooldid_solved_issues.md` |
| `mini_assessment_2026-08-23.md` | `wooldid_solved_issues.md` |
| `documentation_review_2026-09-10.md` | the rulings table at the head of `open_defects.md` |
| `efficiency_review_2026-08-31.md`, `mlib_decision_2026-08-31.md` | `guard_cost_2026-09-16.md` for runtime; the decision itself is in history |
| `inv1_v_es_ame_2026-08-24.md` | the code comments at the `e(V_es_ame)` construction |

Several surviving documents still cite the removed files by name. Those citations are not broken -
they name a document that git still holds - and they were deliberately left rather than rewritten,
since rewriting forty citations by hand invites exactly the sort of unverified edit this cleanup
was meant to stop.

## A standing rule this directory learned the hard way

**A status line is only true on the day it is written.** Three separate documents here carried
lines - "nothing in it has been decided", "OPEN, ruled, NOT started", "Nothing below is implemented
yet" - that were accurate when typed and badly wrong within days. If you write a status claim into
a document that will outlive the week, date it and say what would make it false.
