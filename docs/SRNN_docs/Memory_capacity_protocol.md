Here’s a protocol that’s essentially their memory-capacity method, adapted to your continuous-time reservoir with SFA + STD.

I’ll first restate the measurement idea, then give you a concrete MATLAB-style protocol you can drop on top of your existing class.

---

## 1. Conceptual summary (what they do)

In the paper, the memory task is:

* Drive the reservoir with a scalar random input (u(t)) uniformly in ([0,1]).

* For each delay (d) (measured in time steps), train a **separate linear readout** to reconstruct (u(t-d)) from the reservoir state at time (t).

* For each (d), compute a coefficient of determination (R^2_d) between the true delayed input and the readout output:

  [
  R^2_d
  = \frac{\mathrm{cov}\big(u(t-d),,y_d(t)\big)^2}
  {\mathrm{var}\big(u(t)\big),\mathrm{var}\big(y_d(t)\big)}.
  ]

* The **memory capacity** is the sum of these (R^2_d) values up to some maximum delay (d_{\max}) (they cap at (d_{\max}=70)).

So your job is to:

1. Run your continuous-time reservoir with input (u(t)),
2. Sample the state,
3. Fit linear readouts for each delay,
4. Compute (R^2_d) and sum.

---

## 2. Your reservoir equations

You’ve got (continuous time):

[
\begin{aligned}
\dot{x}*i &= \frac{-x_i + u_i(t) + \sum*{j=1}^{J} w_{ij}, b_j r_j}{\tau_d}[4pt]
r_i &= \phi!\left(
x_i - a_{0_i}
- c \sum_{k=1}^{K} a_{ik}
\right)[4pt]
\dot{a}*{ik} &= \frac{-a*{ik} + r_i}{\tau_k}[4pt]
\dot{b}*i &= \frac{1-b_i}{\tau*{rec}}
- \frac{b_i, r_i}{\tau_{rel}}.
\end{aligned}
]

You can treat (r(t)) (possibly plus other internal variables if you want) as the reservoir state for readout.

---

## 3. Protocol to measure memory capacity

### 3.1. Fix basic simulation settings

Choose:

* Time step for sampling / input updates: (\Delta t)

  * E.g. (\Delta t) smaller than or comparable to (\tau_d); if you already have a fixed ODE step, just sample every fixed number of integrator steps.
* Total number of samples:

  * Washout: (T_{\text{wash}}) (e.g. 1000–2000 steps)
  * Training: (T_{\text{train}}) (e.g. 5000–10000 steps)
  * Test: (T_{\text{test}}) (e.g. 5000 steps)
* Maximum delay: (d_{\max}=70) (to match the paper).

Define discrete sampling times (t_n = n\Delta t), with (n = 1,\dots,T_{\text{total}}) where (T_{\text{total}} = T_{\text{wash}} + T_{\text{train}} + T_{\text{test}}).

### 3.2. Input signal and input weights

1. Generate a scalar input sequence:
   [
   u_n \sim \mathcal{U}(0,1),\quad n=1,\dots,T_{\text{total}},
   ]
   as in the paper.

2. Define a fixed input weight vector (W_\text{in} \in \mathbb{R}^{N}):

   * Pick some fraction (f_{\text{in}}) of neurons that receive input (e.g. 30%, matching their (f_\text{in}=0.3)).
   * For those neurons, draw weights i.i.d. from
     [
     W_{\text{in},i} \sim \mathcal{U}\left[-\frac{\sigma_{\text{in}}}{2},\frac{\sigma_{\text{in}}}{2}\right],
     ]
     with (\sigma_{\text{in}}) your input scaling parameter (they tune (\sigma_{\text{in}}) per task).

3. Map scalar input into your equations via
   [
   u_i(t_n) = W_{\text{in},i}, u_n.
   ]
   In MATLAB terms, at each sample time:

   ```matlab
   u_vec = W_in * u(n);   % N x 1
   ```

   and feed `u_vec` into your class as the `u_i(t)` term.

Use the **same** (W) and (W_\text{in}) in all conditions (baseline vs SFA/STD) so differences in memory capacity are attributable to adaptation, not different random initializations.

### 3.3. Run the reservoir and record states

For a given parameter condition (e.g. “no adaptation”, “SFA only”, “SFA+STD”):

1. Reset all states: (x_i=0), (a_{ik}=0), (b_i=1) (or whatever baseline you prefer).

2. For each (n = 1,\dots,T_{\text{total}}):

   * Set (u_i(t)) according to (u_n) as above (piecewise-constant over ([t_n,t_{n+1}))).
   * Integrate your ODE from (t_n) to (t_{n+1}) with your class (Euler/RK etc).
   * Compute (r_i(t_{n+1})) via your nonlinearity.
   * Store:

     * either just (r(t_{n+1})) (standard RC practice),
     * or an augmented state such as ([r; x; a; b]) if you want a richer feature set (but keep this fixed across conditions if you want clean comparisons).

3. Discard the first (T_{\text{wash}}) samples to remove transients.

Collect:

* (R_{\text{all}} \in \mathbb{R}^{N_\text{feat} \times T_{\text{eff}}}): reservoir features over the remaining (T_{\text{eff}} = T_{\text{train}} + T_{\text{test}}) steps.
* (u_{\text{eff}} \in \mathbb{R}^{T_{\text{eff}}}): corresponding inputs.

Split into training and test segments:

* Training: indices (1,\dots,T_{\text{train}})
* Test: indices (T_{\text{train}}+1,\dots,T_{\text{eff}})

### 3.4. Build delayed targets

For each delay (d=1,\dots,d_{\max}):

1. Define training targets
   [
   u^{(d)}*{\text{train}}(n)
   = u*{\text{eff}}(n-d),
   ]
   for (n = d+1,\dots,T_{\text{train}}).
   (You’ll lose the first (d) samples for that delay.)

2. Similarly define test targets:
   [
   u^{(d)}*{\text{test}}(n)
   = u*{\text{eff}}(T_{\text{train}} + n - d),
   ]
   for (n = d+1,\dots,T_{\text{test}}).

In code, you can implement this with index shifts rather than writing explicit formulas each time.

### 3.5. Train linear readouts (ridge regression)

For each delay (d):

1. Extract the corresponding slices of the reservoir features, aligned with the target:

   * Training design matrix
     (X^{(d)}*{\text{train}} \in \mathbb{R}^{T*{\text{train},d} \times N_\text{feat}})
     (rows are time points, columns are features).
   * Training target vector
     (y^{(d)}*{\text{train}} \in \mathbb{R}^{T*{\text{train},d}}).

2. Fit a linear readout with ridge regression, as in the paper:
   [
   w^{(d)} = \underset{w}{\arg\min},
   |X^{(d)}*{\text{train}} w - y^{(d)}*{\text{train}}|^2

   * \eta |w|^2,
     ]
     with regularization (\eta) (they use (\eta=10^{-7})).

   Closed-form:
   [
   w^{(d)} = \Big( {X^{(d)}*{\text{train}}}^\top X^{(d)}*{\text{train}}
   + \eta I \Big)^{-1}
   {X^{(d)}*{\text{train}}}^\top
   y^{(d)}*{\text{train}}.
   ]

3. On the test set, compute predictions
   [
   \hat{y}^{(d)}*{\text{test}} = X^{(d)}*{\text{test}}, w^{(d)}.
   ]

### 3.6. Compute (R^2_d) and memory capacity

For each delay (d), compute the coefficient of determination on the **test** segment as in Eq. (6):

[
R^2_d
= \frac{\mathrm{cov}\big(u^{(d)}*{\text{test}},,\hat{y}^{(d)}*{\text{test}}\big)^2}
{\mathrm{var}\big(u^{(d)}*{\text{test}}\big),
\mathrm{var}\big(\hat{y}^{(d)}*{\text{test}}\big)}.
]

Then define total memory capacity (for that parameter condition):

[
\text{MC} = \sum_{d=1}^{d_{\max}} R^2_d.
]

Use the same (u_n), (W), (W_\text{in}), and readout procedure across conditions.

---

## 4. Comparing “with” vs “without” SFA / STD

You want to compare at least two conditions:

1. **Baseline reservoir (no SFA, no STD)**

   * Set (c = 0) and do not evolve (a_{ik}) (or equivalently set all (\tau_k\to\infty) and initial (a_{ik}=0)).
   * Fix (b_i(t) \equiv 1): either do not integrate (\dot{b}*i), or set parameters so that the ODE collapses to (b_i=1) (e.g. infinite (\tau*{\text{rec}}), zero depression term).

2. **SFA only**

   * Include the (a_{ik}) dynamics with finite (\tau_k) and (c>0).
   * Keep (b_i(t)\equiv 1) as above.

3. **SFA + STD**

   * Full system as written, with finite (\tau_k), (\tau_{\text{rec}}), (\tau_{\text{rel}}).

For each condition:

* Re-run the full protocol (Sections 3.3–3.6),
* Average MC over multiple random seeds (different (W) and/or (W_\text{in})) to get stable estimates.

You can also explore how MC depends on adaptation time constants ((\tau_k, \tau_{\text{rec}}, \tau_{\text{rel}})) and on the gain (c), holding the rest fixed.

---

## 5. Practical tuning notes

* **Dynamical regime:** As in the paper, MC will be maximal when the reservoir avoids silent, saturated, or globally synchronized regimes and stays in an “intermediate” firing-rate regime.

  * You can monitor mean firing rate and entropy (optional) to ensure you’re not stuck in trivial dynamics.
* **Input scaling (\sigma_{\text{in}}):** Memory task in their discrete RC prefers relatively small input scaling (more linear regime); nonlinear tasks require larger (\sigma_{\text{in}}).

  * For memory capacity, start with small (\sigma_{\text{in}}) and adjust to maximize MC in the baseline condition, then keep it fixed while you vary SFA/STD.
* **Sampling vs. internal time constants:** Ideally (\Delta t) is smaller than the smallest of (\tau_d, \tau_k, \tau_{\text{rec}}, \tau_{\text{rel}}) so you resolve the dynamics; but MC is defined in units of the sampling step, so once you choose (\Delta t), keep it fixed across all experiments.

---

If you like, next step I can help you turn this into concrete MATLAB scaffolding (function signatures, matrix shapes, ridge solver) that plugs directly into your existing class.
