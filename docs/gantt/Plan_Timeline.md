# Plan timeline

A Gantt chart of the build and test steps of the final-draft work: where
each step stands, what it took, and an estimate of what remains. Open
`plan_timeline.html` by double-click; it is one self-contained file and
needs no network or server.

## What it is for

In Tom's words (2026-11-12): "I don't want any calendar target dates. I
just want estimates to aid in planning. I am not shooting for a specific
done datetime. I just want visibility of the timeline of the current plan.
Agents which are building or testing are not beholden to the gantt chart.
It is just to help the user estimate how the plan is shaping up."

- **It is for visibility only.** It shows the current plan and an estimate
  of how long its steps take. It sets no finish date for anything.
- **To every agent building or testing: nobody is bound by it.** A step
  takes as long as doing it properly takes. Never shorten, skip, reorder or
  hurry work to match a bar. A step that runs longer than its estimate is
  not a problem to fix, only an estimate to correct at the next update.
  Quality comes first and the plan is the guide (`CLAUDE.md`, "Current
  work").
- **The estimates change.** They are corrected after each merge, and
  blocks B, C and D will move most as they are planned.
- **No calendar dates for future work.** Past work is shown at its real
  times from the git log, because that is a record. Everything to the
  right of the "now" line is labelled in elapsed working time from now
  (+1:00, +2:00, ...), as if the work ran without pause: it excludes
  nights, idle time and waiting for Tom. No finish date or time is printed
  for a block or the plan; the page gives remaining estimated working
  hours per block instead.

## Files

| File | What it is |
|---|---|
| `timeline.toml` | The data: one `[[task]]` per bar. The only file edited by hand. Its header lists every field. |
| `build_timeline.py` | Reads the data, validates it, schedules the future, writes the page. Standard library only. |
| `plan_timeline.html` | The generated page. Committed, so it can be opened without running anything. |
| `pyproject.toml`, `uv.lock`, `.python-version` | The uv project on Python 3.13, no third-party dependencies. `.venv` is gitignored. |

## How to read the chart

- **Task code** `<stage>_<step>_<kind>[n]`, e.g. `A_0c_T`:
  - *stage*: the block letter, A, B, C or D (blocks of
    `../StochasticPlasticDynamicalSystemPaper/docs/final_draft_to_do.md`).
  - *step*: the step as the plan numbers it. Block A uses Plan A's numbers
    (`docs/Plan_A_chained_vessicles.md` section 7: 0a, 0b, 0c, 0d, 1, 2, 3,
    4, 5, 7, 8; then the lead's 9 (the spread), 10, 11, 12 for the finish).
    Step 0 of block A also holds the work of 2026-11-10 that came before
    Plan A (test audit, regression tiers). Block C uses Plan C's steps
    (`docs/Plan_C_memory_capacity_gain_bias.md` section 6: C0 the three
    fixes and their re-captures, C1 the stage, C2 the pilot, C3 the main
    runs, C4 figures and report; then 5 the close-out and Tom's review, 6
    the medium run, 7 its figures). `C_3_U` is Tom's decision of
    2026-11-13 morning that made block C a quick test (FD-010, decision
    21); the main run's lost-login interruption, 03:48 to 08:11, is
    left as a gap between bars, not drawn as working time. Block B's step 0 is its pilots and
    their report; its later steps, and all of block D, are not planned in
    detail and are numbered 1, 2, ... in order.
  - *the medium run*: there is one, for blocks A, B and C together, drawn
    in block C (`C_6_U` Tom decides to start it, `C_6_B` the run): Tom,
    2026-11-12, "we will not be doing a medium run until we get done with
    stage C", so it follows block C's close-out and Tom's review, and
    block B's fast run.
  - *kind*: P = Plan, B = Build, T = Test, U = User: nothing proceeds
    until Tom has decided, reviewed, or talked it through (a hold point,
    Tom fixing the tolerances). Writing a plan or a design stays P, a run
    that produces results stays B even when only Tom starts it.
  - *n*: an optional number when one step has several tasks of the same
    kind (`A_4_B1` to `A_4_B5` are FD-006's five display packages).
  The code is the task's id in the data; dependencies refer to it.
- **Status** is shown by a word and a symbol in the status column and by
  the bar's pattern, not by colour alone: solid grey = done (actual
  times); amber stripes = running; outlined = planned; purple dots =
  waiting on Tom (a decision, a review, or a run only Tom starts); dashed
  outline and "(guess)" = rough, not yet planned in detail. The kind U and
  the status "waiting on Tom" are different things: a task's kind says
  what the step is and never changes, its status says where it stands and
  changes as the work moves (a U task is done once Tom has spoken, and a
  B task such as Tom's medium run waits on Tom without being a U).
- **Beside each bar**: done tasks show how long they took ("took 1 h 06
  min", with ≈ when the start was inferred), and "est. X, took Y" when an
  estimate was recorded before they ran. Future tasks show their estimate.
- **Idle stretches** in the past (3 h or more with no recorded work, such
  as nights) are drawn as a narrow grey "idle" column.
- **Hover** over a bar or a row (or tab to it) for what the step does and
  why, who does it, the issue, status, times, the prior estimate, what it
  waits for, how it is checked, and the merge commit.
- **Estimates against actuals**, near the top of the page: how many
  finished tasks had an estimate recorded before they ran, and the median
  and range of actual / estimate. It is information for planning, not a
  score.

## Where the numbers come from

**Actual times of finished tasks** come from git (`git log --merges
--date=iso main`, the branch commits, and the lead's status commits in
`docs/ISSUES.md` / `docs/WORKLOG.md`). `start` and `end` keep the exact
times; `took_min` is the measured minutes (the script checks it against
them); bars are drawn on whole quarter hours (start rounded down, end up).
Each done task says in `rule` how it was measured:

- `brief to merge`: start = the lead's commit that set the issue in
  progress or recorded the brief; end = the merge commit. Used for
  worker builds. (First commit on the branch to merge would understate
  them: a worker usually commits once, at the end, so FD-017 would read
  2 min.)
- `previous merge to merge`: the start was not recorded; it is taken as
  the merge or lead commit that made the task possible. Approximate.
- `merge to gates commit`: from the merge to the lead's next commit
  recording that the gates passed. Used for the lead's tests.
- `commits on main`: work committed straight to main (before the worker
  arrangement, and the lead's own steps); first to last commit, or a
  status commit.
- `worklog`: times stated in `docs/WORKLOG.md`.

`approx_start = true` marks an inferred start; the chart then writes ≈.

**Estimates.** `estimate` is the working estimate in minutes, a multiple
of 15, for tasks not finished. For block A they are the lead's; where a
duration is a guess (`guess = true`) the hover box and the bar say so.
Block C's are the lead's (from Plan C and its measured
pilot) up to Tom's review; blocks B and D, and the medium run, are
guesses.

**Estimate of record.** When a task starts, its estimate at that moment
is copied into `estimate_of_record` (with `estimate_of_record_source`)
and frozen. When it finishes, the actual is recorded beside it, and
neither is edited afterwards. Estimates are never filled in after the
fact: a finished task with no estimate recorded before it ran shows its
actual time alone and "no prior estimate". This is what keeps the
estimates-against-actuals summary honest.

## How the future is scheduled

The script places it; nothing about the future is typed in. Done tasks
sit at their actual times. A running task starts at its actual start and
lasts its estimate; if it has already run longer, its bar extends to a
quarter hour past now and a dotted tick marks where the estimate ended.
Every other task starts at the later of now (rounded up to the quarter
hour) and the end of the last task it depends on, and lasts its estimate.
Resources are not modelled: tasks whose dependencies allow it are drawn in
parallel, so the lead's serial work must be expressed as dependencies.

## Update procedure (after every merge into main)

Run from `docs/gantt/`. Do not start MATLAB; this needs only git and uv.

1. **Read what happened.** `git log --date=iso --format='%h %ad %s' -30
   main` and `git log --merges --date=iso main`; the top entry of
   `docs/WORKLOG.md`; the issue's status in `docs/ISSUES.md`.
2. **Mark what started.** For each task now under way: `status =
   "running"`, `start =` the time it began (the commit or worklog time that
   shows it; `approx_start = true` if inferred), and copy its current
   `estimate` into `estimate_of_record` with `estimate_of_record_source`
   saying where that estimate was recorded. If no estimate existed before
   it began, leave `estimate_of_record` out. A task that is under way but
   waiting on Tom (a hold point, his own run) keeps `status = "tom"` and
   gets a `start` and a frozen estimate the same way; it is then drawn
   like a running task.
3. **Mark what finished.** `status = "done"`, `end =` the merge (or the
   lead's gates commit for a test), `took_min =` end - start in minutes,
   `rule =` which of the rules above gave the times, `merge`, `commits`.
   Never change `estimate_of_record` once set, and never add one to a
   finished task.
4. **Correct the plan.** Adjust `estimate` of tasks not yet started; add,
   split or remove future tasks if the plan changed (a new code follows
   the pattern; update every `depends` that pointed at a removed code).
   For blocks B, C, D, replace rough bars with real steps as
   `ISSUES.md` gains their issues, and drop `guess` where the lead gives an
   estimate. Write plain sentences for `what` and `gate`: the reader is a
   neuroscientist following the build, not the code.
5. **Rebuild and check.** `uv run python build_timeline.py` (it uses the
   clock for "now"; `--now 2026-11-12T10:45-06:00` reproduces a given
   page exactly). It refuses data with an unknown dependency, a cycle, a
   duplicate or malformed code, an estimate that is not a multiple of 15,
   a done task without times, `took_min` or `rule`, a frozen estimate on a
   task that has not started, or a name too long for its column. Open the
   page, then commit `timeline.toml` and `plan_timeline.html` together
   with a message naming the merge it records.

The words "deadline", "target", "due", "on track", "late" and "behind" do
not belong in the data or on the page; check with
`grep -n -i -E "\b(due|deadline|target|on track|late|behind)\b" timeline.toml plan_timeline.html`.
