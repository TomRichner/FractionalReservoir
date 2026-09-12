# Reading list: covariant Lyapunov vectors, spectrum-derived measures, transient divergence

For the next additions to the stability analysis (KS entropy rate, attractor
dimension, transient divergence, leading-vector localisation, covariant
Lyapunov vectors). Companion to `topK_reading_list.md` (already fetched:
Dieci 1995, Engelken 2023, Geist 1990, Skokos 2008) and to
`Cov_Lya_vectors.md`.

Please put PDFs in `C:\Users\m218089\Desktop\github_repos\PDF2md\StabilityPaper\`
as before. Items marked **open** are on arXiv and I can fetch them myself;
they are listed so the set is complete in one place.

## Needed (paywalled or not reliably open)

1. **Ginelli, Poggi, Turchi, Chaté, Livi, Politi (2007).** Characterizing
   dynamics with covariant Lyapunov vectors. *Phys. Rev. Lett.* 99, 130601.
   The algorithm (forward QR, backward push of an upper-triangular matrix
   through the stored R's). Primary citation for the CLV computation.
2. **Hennequin, Vogels, Gerstner (2014).** Optimal control of transient
   dynamics in balanced networks supports generation of complex movements.
   *Neuron* 82, 1394-1406. Non-normal transient amplification in a stable
   network; the citation for using the numerical abscissa / transient gain
   as the "transient divergence" measure. (Also their 2012 *Phys. Rev. E* 86,
   011909, "Non-normal amplification in random balanced neuronal networks",
   if easy.)
3. **Wolfe & Samelson (2007).** An efficient method for recovering Lyapunov
   vectors from singular vectors. *Tellus A* 59, 355-366. The alternative
   CLV algorithm; useful as a cross-check of Ginelli's on one network.
4. **Eckmann & Ruelle (1985).** Ergodic theory of chaos and strange
   attractors. *Rev. Mod. Phys.* 57, 617-656. The standard source for
   Pesin's identity (h_KS = sum of positive exponents), the Kaplan-Yorke
   conjecture, and what the Lyapunov dimension does and does not mean.
   Citation for h_KS and D_KY in the Methods.
5. **Frederickson, Kaplan, Yorke, Yorke (1983).** The Liapunov dimension of
   strange attractors. *J. Differential Equations* 49, 185-207. The
   Kaplan-Yorke formula's original statement (Kaplan & Yorke 1979 is a
   conference proceedings and hard to get; this is the citable one).
6. **Pazó, Szendro, López, Rodríguez (2008).** Structure of characteristic
   Lyapunov vectors in spatiotemporal chaos. *Phys. Rev. E* 78, 016209.
   Localisation of CLVs and the participation-ratio measure Engelken uses.
7. **Froyland, Hüls, Morriss, Watson (2013).** Computing covariant Lyapunov
   vectors, Oseledets vectors, and dichotomy projectors: a comparative
   numerical study. *Physica D* 247, 18-39. Compares Ginelli, Wolfe-Samelson
   and two SVD-based approaches; convergence behaviour and pitfalls.
8. **Lajoie, Lin, Shea-Brown (2013).** Chaos and reliability in balanced
   spiking networks with temporal drive. *Phys. Rev. E* 87, 052901.
   Intermittent reliability in chaotic networks, local Lyapunov exponents;
   the "transient divergence amplifies information" citation the
   Introduction already names (Lajoie). Its 2014 *Front. Comput. Neurosci.*
   companion (8:123, "Structured chaos shapes spike-response noise entropy")
   is open access.

## Open on arXiv (I will fetch these; listed for completeness)

9. **Kuptsov & Parlitz (2012).** Theory and computation of covariant
   Lyapunov vectors. *J. Nonlinear Sci.* 22, 727-762. arXiv:1105.5228. The
   review: theory, both algorithms, adjoint CLVs, angles. The one to read
   first.
10. **Ginelli, Chaté, Livi, Politi (2013).** Covariant Lyapunov vectors.
    *J. Phys. A* 46, 254005. arXiv:1212.3961. Tutorial-length version of
    the 2007 method with worked examples and the hyperbolicity (angle) test.
11. **Noethen (2019).** A projector-based convergence proof of the Ginelli
    algorithm for covariant Lyapunov vectors. *Physica D* 396, 18-34.
    arXiv:1802.08461. Why and how fast the backward pass converges.
12. **(2025/26) On the efficient numerical computation of covariant Lyapunov
    vectors.** arXiv:2512.23002. How long the forward and backward
    transients need to be before the vectors are trustworthy -- directly
    what we must choose for the backward warm-up.
13. **Engelken, Wolf, Abbott (2023).** Lyapunov spectra of chaotic recurrent
    neural networks. *Phys. Rev. Research* 5, 043044. arXiv:2006.02427.
    Already fetched; the reference for extensive chaos, D_KY and h_KS versus
    gain, and the delocalised first CLV (participation ratio N/3).

## Optional

14. **Trefethen & Embree (2005).** *Spectra and Pseudospectra*, Princeton.
    Chapters on the numerical abscissa and transient growth bounds, if we
    want the transient-divergence measure stated rigorously. A book; the
    relevant definitions are also in Hennequin 2014's supplement.
15. **Farrell, Recanatesi, Moore, Lajoie, Shea-Brown (2022).** Gradient-based
    learning drives robust representations in recurrent neural networks by
    balancing compression and expansion. *Nat. Mach. Intell.* 4, 564-573.
    Finite-time expansion/compression in trained RNNs; possibly relevant to
    the PyTorch training-rate result.
