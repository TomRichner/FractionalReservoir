# User Notes

Notes TR asks the agent to write down so they are not forgotten, and so he can
come back to them later if he decides to. **Nothing here is a commitment.** An
entry is a reminder, not a task: it is not a backlog item, it is not approval to
act, and it does not become work until TR says so.

This file is for TR to read. The agent writes an entry only when asked to, and
otherwise leaves it alone. Commit messages remain the primary record of what was
actually done.

Each entry is stamped with the date, branch @ commit, hostname, and the agent
and session that wrote it. Newest first.

---

## Update CLAUDE.md's equations, and write an SFA version of the parameter table

| | |
|---|---|
| Noted | 2026-08-20 11:05 · `dev` @ `fb3639e` · R5469844 · Claude Code (Opus 5), session 5f04ea68 |

Three related pieces of documentation drift, raised by TR. **Not started.**

1. **The equation block in `CLAUDE.md` ("What this project is") is out of date.**
   It is written for `SRNNModel2`: a scalar `b_i` depressing every outgoing synapse,
   no per-route synapses, no STF, and no noise term (`dx_i` is shown as a plain ODE
   even though the model now supports additive Wiener noise on `x`).

2. **`SRNNCellTypePairs` is the primary class and CLAUDE.md should say so.** It is
   currently presented alongside `SRNNModel2` as a duck-typed sibling, with
   "default to `SRNNModel2` + `ParamSpaceAnalysis2`" as the closing advice. That no
   longer reflects where the work is.

3. **Write an SFA version of `docs/EquationsParametersDocs/cell_type_pair_parameter_table.md`.**
   That document covers STD and STF per route, and the additive Wiener noise on `x`,
   but has **no mention of SFA at all** — no `a_{i,k}` states, no `c`, no `tau_a`.
   Wanted: the same treatment for adaptation, in the per-cell-type notation.

---

## Separate the network seed from the noise seed for reps

| | |
|---|---|
| Noted | 2026-08-20 10:32 · `dev` @ `13d66d6` · R5469844 · Claude Code (Opus 5), session 5f04ea68 |
| Origin | FR-005, raised 2026-08-13 · `dev` @ `ac59b42` · R5456622 |

`noise_seed` derives from `rng_seeds(1)`, so a rep varies the network **and**
the noise realisation together. That is the right default — it is what makes
the noise path shared across the four adaptation conditions at a grid point —
but it means the within-level spread in a reps histogram mixes two sources and
neither can be attributed.

Wanted: a way to hold one fixed while varying the other, so "how much of this
spread is the network and how much is the noise" becomes answerable.
