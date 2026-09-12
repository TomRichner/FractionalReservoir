# Reading list for the top-K Lyapunov spectrum implementation

Papers to fetch and convert to markdown (Mistral OCR) before designing the
top-K routine described in `Lyapunov_estimation_methods.md`, Section 7.
Drop the converted files in `docs/notes/refs/`. 2026-09-11.

## Essential

1. **Engelken, Wolf & Abbott (2023).** Lyapunov spectra of chaotic recurrent
   neural networks. *Phys. Rev. Research* 5, 043044.
   Open access: https://arxiv.org/abs/2006.02427
   Why: same class of rate network at our scale; full spectra; their
   convergence controls, re-orthonormalisation interval and warm-up; the
   extensivity and symmetry results to reproduce as a check.

2. **Skokos (2010).** The Lyapunov characteristic exponents and their
   computation. *Lect. Notes Phys.* 790, 63-135.
   Open access: https://arxiv.org/abs/0811.0882
   Why: the standard review; K-vector Benettin / Shimada-Nagashima algorithm
   in pseudo-code; error and convergence discussion.

3. **Dieci & Van Vleck (1995).** Computation of a few Lyapunov exponents for
   continuous and discrete dynamical systems. *Appl. Numer. Math.* 17,
   275-291. https://doi.org/10.1016/0168-9274(95)00033-Q
   Why: the partial-spectrum case specifically, including the caveats on
   continuous orthonormalisation when K < N. Paywalled.

## Optional

4. **Carbonell, Biscay & Jimenez (2010).** QR-based methods for computing
   Lyapunov exponents of stochastic differential equations. *Int. J. Numer.
   Anal. Model. Ser. B* 1, 147-171.
   Why: justification for QR on a noisy fiducial trajectory; pitfalls of
   interpolating a stochastic trajectory.

5. **Geist, Parlitz & Lauterborn (1990).** Comparison of different methods for
   computing Lyapunov exponents. *Prog. Theor. Phys.* 83, 875-893.
   Open access: https://doi.org/10.1143/PTP.83.875
   (also https://carretero.sdsu.edu/teaching/M-638/lectures/ProgTheorPhys_1990-Geist-875-93.pdf)
   Why: practical comparison of discrete vs continuous methods; the thin-QR
   (rectangular Q) formulation for K < N.

## Not needed

Benettin, Galgani, Giorgilli & Strelcyn (1980) and Shimada & Nagashima
(1979): the review in item 2 covers their content.
