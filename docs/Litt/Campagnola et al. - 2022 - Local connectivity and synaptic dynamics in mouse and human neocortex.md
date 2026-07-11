RESEARCH

# RESEARCH ARTICLE SUMMARY

# NEUROSCIENCE

# Local connectivity and synaptic dynamics in mouse and human neocortex

Luke Campagnola†, Stephanie C. Seeman†, Thomas Chartrand, Lisa Kim, Alex Hoggarth, Clare Gamlin, Shinya Ito, Jessica Trinh, Pasha Davoudian, Cristina Radaelli, Mean-Hwan Kim, Travis Hage, Thomas Braun, Lauren Alfiler, Julia Andrade, Phillip Bohn, Rachel Dalley, Alex Henry, Sara Kebede, Alice Mukora, David Sandman, Grace Williams, Rachael Larsen, Corinne Teeter, Tanya L. Daigle, Kyla Berry, Nadia Dotson, Rachel Enstrom, Melissa Gorham, Madie Hupp, Samuel Dingman Lee, Kiet Ngo, Philip R. Nicovich, Lydia Potekhina, Shea Ransford, Amanda Gary, Jeff Goldy, Delissa McMillen, Trangthanh Pham, Michael Tieu, La'Akea Siverts, Miranda Walker, Colin Farrell, Martin Schroedter, Cliff Slaughterbeck, Charles Cobb, Richard Ellenbogen, Ryder P. Gwinn, C. Dirk Keene, Andrew L. Ko, Jeffrey G. Ojemann, Daniel L. Silbergeld, Daniel Carey, Tamara Casper, Kirsten Crichton, Michael Clark, Nick Dee, Lauren Ellingwood, Jessica Gloe, Matthew Kroll, Josef Sulc, Herman Tung, Katherine Wadhwani, Krissy Brouner, Tom Egdorf, Michelle Maxwell, Medea McGraw, Christina Alice Pom, Augustin Ruiz, Jasmine Bomben, David Feng, Nika Hejazinia, Shu Shi, Aaron Szafer, Wayne Wakeman, John Phillips, Amy Bernard, Luke Esposito, Florence D. D'Orazi, Susan Sunkin, Kimberly Smith, Bosiljka Tasic, Anton Arkhipov, Staci Sorensen, Ed Lein, Christof Koch, Gabe Murphy, Hongkui Zeng, Tim Jarsky*

**INTRODUCTION:** The mammalian neocortex is believed to act as the computational substrate for our highest cognitive abilities, particularly the ability to model the world around us and predict the effects of our actions. Many aspects of cortical structure are repeated across brain regions and conserved across species, suggesting a general-purpose approach to cognition. Within this repeating structure, neurons influence the formation of synaptic connections based on their cell type-specific biases. This results in a stereotyped network architecture in which synapse properties and connectivity are strongly influenced by cell type.

Synapses between cell types transmit information in a way that is highly stochastic and depends on the prior history of activity. The dynamic properties of synapses are also strongly dependent on both the pre- and postsynaptic cell types, suggesting an important role in cortical function. This provides a major source of computational diversity that is often absent in neuroscience modeling studies as well as modern machine-learning architectures.

Neurons are broadly grouped into excitatory and inhibitory classes, each of which can be divided into more specific subclasses. Cortical inhibitory neurons, for example, are commonly divided into Pvalb, Sst, and Vip subclasses and are distributed broadly across cortical layers. In contrast, most excitatory cell subclasses occupy narrower regions across cortical layers.

**RATIONALE:** Understanding the connectivity among cell subclasses and the computations performed at their synapses is an essential step to understanding cortical circuit function.

This has led to experiments in different species, brain regions, ages, etc. that focus on one or a few circuit elements. These efforts offer an excellent depth of insight to isolated regions of

**Intralaminar circuit diagram among major excitatory (Pyr) and inhibitory (Pvalb, Sst, and Vip) cell subclasses aggregated from all layers of mouse primary visual cortex.** Line (axon) thickness depicts the relative weight (strength and probability of connection) of connections between subclasses. Black dots indicate connections that are stronger in layer 2/3 compared with layer 5. Axon color shows the spike-to-spike variance in amplitude of synaptic signaling, which is strongly cell subclass dependent. Excitatory synapse variance depends on the postsynaptic subclass. Pvalb cells project low-variance connections, whereas Sst and Vip project high-variance connections. More saturated axon colors indicate higher confidence measurements.

the circuit but offer a fragmented view of the circuit as a whole. Further, the difficulty of accessing these historical data discourages reuse and reanalysis. We saw an opportunity to expand upon this history and conduct a more comprehensive and standardized survey than has been attempted in the past. By publishing the analyses, tools, and data that characterize cortical connection properties, we encourage a more unified approach to describing cortical function.

**RESULTS:** We used microelectrodes to record the activity of 1731 synaptic connections across diverse cell types in living tissue samples from mouse and human neocortex. We characterized these connections with the aid of a synaptic release model and found that excitatory dynamics aligned with postsynaptic cell subclass, whereas inhibitory dynamics aligned with the presynaptic subclass in ways that were subclass specific. Synaptic variability was a primary driver of these cross-subclass differences in mouse cortex. Compared with the mouse, human excitatory connections were tuned toward stability and reliability pointing toward species differences in cortical function. We further introduced a method to estimate the rate of connectivity between cell types that accounts for differences between experimental preparations. With this approach, we compared connection probabilities across layer, cell subclass, and species. For instance, connectivity between excitatory cells and Vip inhibitory cells was present in layer 2/3 and absent in layer 5/6 of mouse cortex. Likewise, connection probability among layer 4 excitatory cells was high in mouse cortex and nearly absent in human cortex. Overall, we found that layer-specific circuit representations are necessary to capture the diversity of intralaminar connectivity among cortical cell subclasses.

**CONCLUSION:** We have generated a comprehensive dataset describing synaptic connections within each layer in the mouse and human cortex. Our deep characterization of synapses points toward important principles of cortical organization that relate to current topics in computational neuroscience and machine learning. The open distribution of our data, analyses, and tools enables greater realism in constraining network and synapse models. ■

The list of author affiliations is available in the full article online.

*Corresponding author. Email: timj@alleninstitute.org

†These authors contributed equally to this work.

Cite this article as L. Campagnola et al., Science 375, eabj5861 (2022). DOI: 10.1126/science.abj5861

READ THE FULL ARTICLE AT
https://doi.org/10.1126/science.abj5861

Campagnola et al., Science 375, 1144 (2022) 11 March 2022

1 of 1

Downloaded from https://www.scienc.org at Mayo Clinic on November 13, 2025

RESEARCH

RESEARCH ARTICLE

NEUROSCIENCE

# Local connectivity and synaptic dynamics in mouse and human neocortex

Luke Campagnola¹†, Stephanie C. Seeman¹†, Thomas Chartrand¹, Lisa Kim¹, Alex Hoggarth¹, Clare Gamlin¹, Shinya Ito¹, Jessica Trinh¹, Pasha Davoudian¹, Cristina Radaelli¹, Mean-Hwan Kim¹, Travis Hage¹, Thomas Braun², Lauren Alfier¹, Julia Andrade¹, Phillip Bohn¹, Rachel Dalley¹, Alex Henry¹, Sara Kebede¹, Alice Mukora¹, David Sandman¹, Grace Williams¹, Rachael Larsen¹, Corinne Teeter¹‡, Tanya L. Daigle¹, Kyla Berry¹§, Nadia Dotson¹, Rachel Enstrom¹, Melissa Gorham¹, Madie Hupp¹, Samuel Dingman Lee¹, Kiel Ngo¹, Philip R. Nicovich¹, Lydia Potekhina¹, Shea Ransford¹, Amanda Gary¹, Jeff Goldy¹, Delissa McMillen¹, Trangthanh Pham¹, Michael Tieu¹, La'Akea Sivert², Miranda Walker¹, Colin Farrell¹, Martin Schroedter¹, Cliff Slaughterbeck¹, Charles Cobb³, Richard Ellenbogen⁴, Ryder P. Gwinn⁵, C. Dirk Keene⁶, Andrew L. Ko⁷, Jeffrey G. Ojemann⁸, Daniel L. Silbergeld⁹, Daniel Carey¹, Tamara Casper¹, Kirsten Crichton¹, Michael Clark¹, Nick Dee¹, Lauren Ellingwood¹, Jessica Gloe¹, Matthew Kroll¹, Josef Sulc¹, Herman Tung¹, Katherine Wadhwani¹, Krissy Brouner¹, Tom Egdorf¹, Michelle Maxwell¹, Medea McGraw¹, Christina Alice Pom¹, Augustin Ruiz¹, Jasmine Bomben¹, David Feng¹, Nika Hejazini¹, Shu Shi¹, Aaron Szafer¹, Wayne Wakeman¹, John Phillips¹, Amy Bernard¹, Luke Esposito¹, Florence D. D'Orazi¹, Susan Sunkin¹, Kimberly Smith¹, Bosiljka Tasic¹, Anton Arkhipov¹, Staci Sorensen¹, Ed Lein¹, Christof Koch¹, Gabe Murphy¹, Hongkui Zeng¹, Tim Jarsky¹*

We present a unique, extensive, and open synaptic physiology analysis platform and dataset. Through its application, we reveal principles that relate cell type to synaptic properties and intralaminar circuit organization in the mouse and human cortex. The dynamics of excitatory synapses align with the postsynaptic cell subclass, whereas inhibitory synapse dynamics partly align with presynaptic cell subclass but with considerable overlap. Synaptic properties are heterogeneous in most subclass-to-subclass connections. The two main axes of heterogeneity are strength and variability. Cell subclasses divide along the variability axis, whereas the strength axis accounts for substantial heterogeneity within the subclass. In the human cortex, excitatory-to-excitatory synaptic dynamics are distinct from those in the mouse cortex and vary with depth across layers 2 and 3.

The study of cortical connectivity in mammalian model systems established foundational circuit diagrams among neurons classified by their morphology and electrophysiology (1–3). Refinement of cell classes by their long-range projections and genetic markers facilitated more detailed microcircuit representations in rodents (4–8). By contrast, studies of human cortical connectivity (9–14) have been more limited but offer the opportunity to identify features that may contribute to our unique cognitive abilities. Although our knowledge of cell types and circuits in rodents and humans has advanced, a complete description of the connectivity and

synaptic properties among cell subclasses in each cortical layer is still lacking (15).

Cortical synapses are dynamic, varying their response strength in ways that are highly stochastic and also modulated by the prior history of activity. The dynamical properties are influenced by both the presynaptic and postsynaptic cell types (16–20) and endow neuronal networks with essential sources of computational diversity (21, 22). As models in neuroscience push toward greater biological realism, there is a need for large-scale datasets with a deep description of synaptic dynamics among a diversity of cell types. Likewise, the incorporation of biologically inspired sources of computational diversity has shown promise in advancing machine-learning techniques (23). To meet the needs of the theoretical and experimental research communities, these large datasets need to be open and accessible.

Results

Synaptic physiology pipeline

We performed 1931 experiments (Fig. 1A) in acute brain slices targeting intralaminar connections in layer 2 (L2) to L6 of adult mouse (mean age 46.0 ± 4.6 days) primary visual

cortex (VISp; 1715 experiments) and human frontotemporal cortex from neurosurgical excised tissue (216 experiments). We used transgenic mice that express unique reporters in two subclasses (24) (table S1). Six excitatory subclasses were layer or projection-class specific (Nr5a1 and Rorb, L4; Sim1 or mscRE4-FlpO AAV, L5 extratelencephalic [ET]; Tlx3, L5 intratelencephalic [IT]; Ntsr1, L6 corticothalamic [CT]), and three inhibitory subclasses (Pvalb, Sst, and Vip) were assessed in all targeted layers (Fig. 1A). We probed 23,620 potential connections (mouse: 20,949; human: 2671), of which 1731 were connected by chemical synapses, giving an overall connectivity rate of 7.3% (mouse: 1526, 7.3%; human: 205, 7.7%).

Details of the experiments are described in Fig. 1, B and C, and in the materials and methods. In each experiment, up to eight neurons were selected for simultaneous whole-cell patch-clamp recording primarily in current-clamp conditions, with a subset of stimuli administered in voltage-clamp conditions (Fig. 1B). Stimuli were elicited in each patched neuron while all others were recorded for evidence of a postsynaptic response. Dynamics (short-term plasticity [STP] and variability) of individual connections were assessed from stimulus trains at frequencies ranging from 10 to 200 Hz (Fig. 1B). Cells were stimulated with long current pulses to characterize their intrinsic physiology and later stained with biocytin to characterize their morphology (Fig. 1C).

Intralaminar connectivity in mouse VISp
Distance dependence of connectivity

The likelihood that two neurons are connected decreases with increasing intersomatic distance (25–27). We probed connectivity among cells up to ~200 µm apart (Fig. 2A) but could not ensure that intersomatic distances were sampled equally across different connection elements (specific pre-post combinations). To make reliable comparisons, we modeled the spatial profile of connectivity versus lateral somatic distance with a Gaussian (25) and estimated the peak (pmax) and lateral spread (sigma: σ) of connection probability using maximum likelihood estimation (see the materials and methods; fig. S1B). This analysis provided an estimate of the spatial profile of connectivity while being relatively robust to differences in the sampled intersomatic distribution. We initially classified cell pairs among the four combinations of excitatory (E) and inhibitory (I) cell classes (E→E, E→I, I→E, and I→I). These data were well fit by the Gaussian model (Fig. 2A, solid red line) and indicated low peak connectivity among excitatory cells (5%), with moderate connectivity rates among inhibitory cells (12%) and across E-I cell classes (E→I 12%, I→E 15%). In a tracing study of primary visual cortex

¹Allen Institute for Brain Science, Seattle, WA, USA. ²Byte Physics, Berlin, Germany. ³The Ben and Catherine Ivy Center for Advanced Brain Tumor Treatment, Swedish Neuroscience Institute, Seattle, WA, USA. ⁴Department of Neurological Surgery, University of Washington, Seattle, WA, USA. ⁵Epilepsy Surgery and Functional Neurosurgery, Swedish Neuroscience Institute, Seattle, WA, USA. ⁶Department of Pathology, University of Washington, Seattle, WA, USA. ⁷Regional Epilepsy Center at Harborview Medical Center, Seattle, WA, USA.

*Corresponding author. Email: tim@alleninstitute.org

†These authors contributed equally to this work.

‡Present address: Sandia National Laboratories, Albuquerque, NM. §Present address: KB Cajal Neuroscience, Seattle, WA.

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

1 of 13

Downloaded from https://www.scienc.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

A Pipeline Throughput

C Cell Characterization

|  Cell | Layer | Dendrite | Axon Length (μm)  |
| --- | --- | --- | --- |
|  1 | 2/3 | spiny | > 200  |
|  2 | 2/3 | aspiny | > 200  |
|  3 | 2/3 | aspiny | > 200  |
|  4 | 2/3 | spiny | > 200  |
|  5 | 2/3 | aspiny | > 200  |
|  6 | 2/3 | aspiny | > 200  |
|  8 | 2/3 | sparsely spiny | > 200  |

B Connection Characterization

ii. Pair Processing

iii. Connection Analysis

Fig. 1. Synaptic physiology pipeline. (A) Pipeline throughput summary. Large circles indicate data collection statistics (left, mouse primary visual cortex; right, human). Age distribution for each species is shown in the center. Distributions of recorded cells (outer circles) and experiment count with two through eight simultaneously recorded cells (inner circles) are shown at the bottom. (B) Connection characterization. (i) Multipatch experiment. Stimuli sets used to probe connectivity and dynamics are shown at the top. Fluorescent image (right) shows recorded cells with connectivity diagram overlaid. Example electrophysiology recordings of cells 1, 3, and 5 during the 50-Hz stimulus (orange box). (ii) Pair processing. For each recorded pair, presynaptic-spike aligned postsynaptic responses (black traces) were overlaid and averaged (colored trace corresponding to the same connection in the multipatch experiment panel) in both voltage-clamp (top row) and current-clamp (bottom row) conditions.

From this a connection could be identified. (iii) Connection analysis. When a connection was identified, the average (colored) was fit (black) and metrics such as amplitude and rise time (orange lines) were extracted from the fit. STP was quantified as a ratio from fits of individual PSPs. (C) (i) Intrinsic Electrophysiology. Long pulse steps were applied to quantify intrinsic cell properties and electrical connections (fig. S4). Hyperpolarizing steps (left) delivered to an example cell probed the subthreshold I–V relationship to quantify metrics such as input resistance and depolarizing suprathreshold sweeps (right) were used to measure spiking and firing rate properties. (ii) Morphology. Cells were filled with biocytin during recording and stained. Images (20×) of the full slice stained with DAPI allowed identification of the cortical layers, and 63× z-stack images were used to assess morphological properties such as dendritic type and axon length. Images and cells are the same as in (B).

L2/3, E→E connections had a wider lateral extent compared with I→E (28). We additionally found that E→E and I→I connections had a similar spatial extent that was larger than both the E→I and I→E connections.

An overabundance of bidirectional connections relative to unidirectional connections can be evidence for connectivity rules that promote the formation of bidirectional connections (29). A simpler explanation, however, is that bidirectionality is not specifically promoted in cortex but rather is an artifact of merging connectivity results across cell types

(30). We quantified the ratio of connected pairs with and without bidirectional connections and observed that reciprocal connections were three to five times more common than expected for a randomly connected network (Fig. 2A, bottom row, red line) among class-level connections (Fig. 2A, bottom row gray line).

Connection probability measurement: Slicing artifacts and detection limits

In the in vitro slice preparation, some connections may be severed, reducing the measured rate of connectivity (25, 31). To mitigate this

bias, we used thick slices (350  \( \mu \) m), targeted cells deep in the slice (median cell depth = 74  \( \mu \) m; fig. S2B), and focused on local (<200  \( \mu \) m apart) intralaminar connections. By modeling the effects of cell depth and axon length on connection probability, we estimated the size of this bias and adjusted our  \( p_{max} \)  measurements accordingly. This yielded an  \( \sim \) 3 to 6% increase in  \( p_{max} \)  when accounting for severed axons (fig. S2, A and D) and an  \( \sim \) 2 to 20% increase in  \( p_{max} \)  when accounting for depth of the targeted cells from the slice surface (fig. S2, B and D).

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

2 of 13

Downloaded from https://www.sciencet.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

**Fig. 2. Mouse connectivity.**

(A) Cells were divided into two main classes, excitatory and inhibitory, and pairs classified into the four combinations of those two classes. Top row: Connection probability as a function of intersomatic distance fit with a Gaussian (red line) and output parameters $p_{\max}$ and sigma ($\sigma$) describe the maximum connection probability and width of the Gaussian. Also shown is the connection probability as a function of intersomatic distance adjusted for presynaptic axon length, depth of the pair from the slice surface, and detection power of connections using a unified model (dashed red line) or by filtering of the data (dotted red line) (see the results section and the materials and methods). Gray line and area are 40-$\mu$m binned average connection probability and 95% CI. Raster below shows distance distribution of connections probed (bottom) and found (top). Bottom row: Normalized rate of reciprocal connections. Probed pairs are unordered and the number of reciprocal connections counted was normalized to the expected value of connection probability squared for a randomly connected network (solid red line). (B) Connection probability matrix for mouse. Connection probability is estimated using a unified model accounting for all corrections as determined from (A) (dashed red line, 'model'). The shading of each element indicates the 95% CI of the data, with higher contrast indicating smaller CI and lower contrast (toward gray) indicating larger CI. The number of

A

ex→in

in→ex

in→in

B

C

D

E

F

connections found out of the number of connections probed are printed in each element. (C to F) Gaussian fit of connection probability versus intersomatic distance (with CI at $p_{\max}$, shaded region) for two contrasting elements with connections found and connections probed raster below. Cross symbol denotes $p_{\max}$ with all adjustments.

Detecting a connection in electrophysiological recordings is influenced by the background electrical noise, which may obscure synaptic responses in the postsynaptic cell. The effect of noise on detection is reduced by averaging the responses to multiple presynaptic action potentials (APs). Over the course of a typical experiment, most cells were stimulated to evoke hundreds of APs; however, this number varied significantly from cell to cell (485 [192, 623] APs, median [interquartile range, IQR]). To account for the opposing

effects of noise and averaging, we quantified the 'detection power' (analogous to the signal-to-noise ratio; see the materials and methods) for each pair of cells that were probed for connectivity. We observed a range of detection powers across pairs, and these data were used to model the effect of detection power on connectivity. A model of the relationship between detection power and connection probability resulted in a 45% increase in estimated $p_{\max}$ for excitatory connections and up to a 90% increase for inhibitory connections (fig.

S2, D and E), suggesting that the observed connectivity rate is affected more by detection power than by slicing artifacts. Detection power is seldom reported and may explain cases in which the observed connectivity in vitro is higher than that in vivo for the same brain area and connections (32, 33).

We extended our model of connection probability on intersomatic distance (Fig. 2A, solid line) to include the effects of slicing and response detection outlined above (see the materials and methods). The model-adjusted

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

3 of 13

Downloaded from https://www.scienc.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

connectivity rate resulted in a one- to twofold increase in estimated $p_{\text{max}}$ at the cell class level (Fig. 2A, dashed line). To confirm these results, we implemented a similar connectivity adjustment by filtering the data for axon lengths, pair depths, and detection power values above their respective median (Fig. 2A, dotted line), where the impact on connection probability is reduced (fig. S2, A to C). This yielded a comparable increase in estimated $p_{\text{max}}$ (Fig. 2A, dotted versus dashed line), so we used the model to adjust peak connectivity estimates at the subclass level, where filtering would induce undersampling. Inhibitory connections received a higher detection power adjustment than excitatory connections. Adjustments for presynaptic axon length and pair depth were relatively uniform across subclasses (fig. S2E).

### Connectivity among cell subclasses

The intralaminar connectivity among mouse VISp subclasses is summarized in Fig. 2B (also see fig. S1A). We found that intralaminar connection probability varied by layer (Fig. 2, B and C), long-range projection target (Fig. 2D), and cell subclass (Fig. 2, E and F). Connectivity was overall higher in superficial layers than in deep layers. For example, the pyramidal to Pvalb connection in L2/3 ($p_{\text{max}} = 0.79$ [0.57, 0.99], 95% confidence interval [CI]) compared with L5 (IT $p_{\text{max}} = 0.22$ [0.12, 0.35]; $\chi^2 P = 7.2 \times 10^{-24}$) (Fig. 2C and fig. S1A). Recurrent connections between excitatory and Vip cells were common in L2/3 (E→Vip 0.38 [0.21, 0.60]; Vip→E 0.11 [0.03, 0.26]) but rare or absent in deeper layers (Fisher's $P = 6 \times 10^{-4}$; fig. S1A).

Within L5, we found several differences between two excitatory projection classes, IT (labeled by Tlx3) and ET (labeled by Sim1 and mscRNA). ET cells overall have more input from local sources relative to IT cells. ET cells have higher recurrent connectivity ($\chi^2 P = 3.17 \times 10^{-6}$; Fig. 2D) and receive unidirectional input from IT cells, consistent with previous results (34). Sst cells also innervate ET cells at a higher rate than IT cells ($\chi^2 P = 3.18 \times 10^{-5}$; Fig. 2, E and F); a similar connectivity pattern was observed in the rat frontal cortex (35).

Sst is thought to avoid connecting with itself (4, 36, 37); however, we observed connections between Sst neurons in every layer (Fig. 2B). The Sst-IRES-Cre driver can sparsely label fast-spiking interneurons that also express Pvalb (38). Uniform manifold approximation and projection (UMAP) of intrinsic properties from inhibitory Cre types showed that slightly less than half of the recurrent Sst connections had both the pre- and postsynaptic cell in a cluster that was spatially distinct from Pvalb cells, suggesting that these connections came from cells that were intrinsically Sst like (fig. S3A). We confirmed that in at least one case, both cells in the pair had axons extending into L1 and sparsely spiny dendrites, consistent with

Martinotti neurons (fig. S3B). We performed experiments using the Patch-seq method (39) and further confirmed that cells that transcriptionally mapped to the Sst subclass do form connections with each other (five connections were found of 66 probed; fig. S3C).

### Electrical connections

In addition to chemical synaptic transmission among cell subclasses, electrical connections formed by gap junctions were also found between inhibitory subclasses (fig. S4A). The likelihood of electrical connectivity as a function of lateral intersomatic distance could be approximated by a Gaussian but with a narrower profile ($\sigma = 74$ $\mu$m; fig. S4B) compared with chemical connections between inhibitory cells ($\sigma = 127$ $\mu$m, Fig. 2A; $P < 0.001$, fig. S4B). This is consistent with previous reports of electrical connections between nearby Pvalb cells showing that the average distance between electrically coupled cells is short, 40 to 80 $\mu$m (40, 41). Pvalb cells showed the highest rate of electrical connections (76/922, -8.2%), whereas those among Vip cells were the most rare (15/908, 1.6%; $\chi^2 P = 1.8 \times 10^{-20}$), in contrast to a previous report of Vip electrical connections, which were more prevalent (7). Most of the electrical connections were found between like subclasses (164/170, 96%) and were bidirectional (148/170, 87%). The distance of the gap junction from the soma coupled with the potential for some rectification of electrical connections (42) could account for the few cases in which reciprocal electrical connections were not observed. The coupling coefficient of electrically coupled pairs was comparable across subclass (fig. S4C) largely due to the lower input resistance of Pvalb cells (fig. S4D). Estimating junctional conductance revealed stronger electrical connections between Pvalb cells (0.38 [0.20, 0.61] nS, median [IQR]) than either Sst (0.23 [0.19, 0.36] nS; Mann-Whitney [MW] $P = 0.02$) or Vip (0.07 [0.03, 0.19] nS; MW $P = 5 \times 10^{-4}$) (fig. S4C).

### Synaptic strength and kinetics

In addition to connectivity, synaptic properties such as strength, latency, and kinetics (rise time and decay tau) determine the impact of a connection on the postsynaptic neuron and, ultimately, on cortical processing. Although these properties have been described previously for some cortical connections (4, 27, 33, 43–45), measuring them systematically across many cell types enabled us to make direct comparisons.

At the class level, inhibitory connections showed short latencies (median = 1.07 ms), slow kinetics, and relatively strong postsynaptic potentials (PSPs). These trends were largely driven by the subclass of the presynaptic cell and Pvalb in particular. Pvalb connections were extremely fast, with submillisecond latencies (Fig. 3A) highlighted

by Pvalb→L6 excitatory connections (0.94 [0.86, 0.99] ms, median [IQR]). Pvalb cells also elicited large resting-state inhibitory postsynaptic potential (IPSP) amplitudes regardless of postsynaptic cell subclass. Presynaptic Sst cells stood out for having some of the slowest kinetics, independent of postsynaptic target (see Sst→Vip rise time: 7.02 [5.65, 9.75] ms, decay: 50.81 [26.97, 142.41] ms; Sst→L5 IT rise time: 7.46 [5.51, 9.92] ms, decay: 58.59 [27.76, 166.83] ms).

In contrast to inhibitory connections, excitatory connections generally had a long latency (median = 1.49 ms; Kolmogorov-Smirnov [KS] compared with inhibitory $P = 1.4 \times 10^{-24}$), fast kinetics, and weak PSPs, all of which related more to the identity of the postsynaptic cell. E→I connections displayed faster rise times than recurrent excitatory connections (E→I 2.75 ms, E→E 3.88 ms, KS $P = 2.97 \times 10^{-12}$; Fig. 3B). Recurrent excitatory connections showed some of the smallest amplitudes in the resting state (e.g., L5 ET→L5 ET, 0.27 [0.13, 0.49] mV), whereas E→I connections were stronger (Fig. 3, D and E) and generated the single biggest PSP (9.51 mV, L2/3 Pyr→L2/3 Sst). E→I connection properties can be further refined by postsynaptic cell subclass. Connections with postsynaptic Pvalb cells had faster kinetics (see L4 Pyr→Pvalb, Fig. 3; rise time: 1.45 [1.28, 1.94] ms, decay tau: 6.92 [5.4, 8.36] ms) than postsynaptic Sst cells (see L2/3 Pyr→Sst, Fig. 3; rise time: 3.91 [2.69, 4.79] ms, decay tau: 19.41 [14.88, 28.48] ms). A dichotomy between Pvalb and Sst was also apparent in the resting-state amplitude with larger excitatory postsynaptic potentials (EPSPs) to Pvalb cells (0.41 [0.22, 0.79] mV) than Sst cells (0.08 [0.05, 0.25] mV); however, the strongest excitatory connections were onto Vip cells found predominantly in superficial layers (L2/3 Pyr→Vip 0.45 [0.26, 0.84] mV). Although resting-state excitation was weakest onto Sst cells, resting-state PSP amplitude is an underestimate of the potential impact on the postsynaptic cell, particularly for connections with strongly facilitating synapses. When we compared the 90th percentile amplitude (Fig. 3E), which measures near-maximal strength, E→Sst connections were comparatively stronger, and even surpassed E→Pvalb amplitudes in some cases (see L2/3 Pyr→Sst versus L4 Pyr→Pvalb). Facilitation onto inhibitory cells further contributes to the longer tails of 90th percentile amplitudes and the rightward shift of E→I amplitudes (0.73 mV) compared with E→E amplitudes (0.34 mV, KS $P = 1.53 \times 10^{-14}$; Fig. 3E, histograms).

A recent survey of cortical connectivity found that connection strength positively correlated with connection probability (6). This result suggests an interesting principle of connectivity but may also result from the reduced detectability of weaker connections. When we adjusted for detection power, connection

Downloaded from https://www.sciencemag.org at Mayo Clinic on November 13, 2025

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

4 of 13

RESEARCH | RESEARCH ARTICLE

Fig. 3. Synaptic strength and kinetics.

Left to right: E-I subclass matrix; excitatory and inhibitory minimum; median, and maximum average traces (light to dark colors); histograms for the major connection classes (E→E, E→I, I→I, I→E); and summary scatter plots for a subset of matrix elements for each metric PSP latency (A), PSP rise time (B), PSP decay tau (C), PSP resting-state amplitude (D), and PSP 90th percentile amplitude (E). In all matrices, inhibitory cells are merged across layers. All matrices are colorized by the median (text in each element) with the saturation scaled by the SE. Two or more pairs were required to fill in an element.

probability was independent of connection strength for both excitatory (weighted Huber regression $r^2 = -0.1$, $P = 0.8$) and inhibitory ($r^2 = 0.2$, $P = 0.2$) connections.

### Synaptic dynamics

The strength and kinetic properties described above characterize the synapse in response to a single presynaptic spike; however, synapses

are highly dynamic. PSP amplitude evolves in predictable ways over the course of milliseconds to seconds because of STP while also being highly stochastic from one response to the next, as quantified by the coefficient of variance. Overall, synaptic dynamics followed a similar pattern to synaptic strength, wherein excitatory connections were most strongly differentiated by the postsynaptic subclass

and inhibitory connections were differentiated by the presynaptic subclass. The STP of a synapse may result in a transient increase (facilitation), decrease (depression), or no change (pseudolinear) in PSP amplitude over the course of a stimulus train, as seen in our data (Fig. 4A) and as described previously (16, 20, 46–48). The time course of recovery from STP is an equally important property of synapses yet one that is

Downloaded from https://www.science.org at Mayo Clinic on November 13, 2025

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

5 of 13

RESEARCH | RESEARCH ARTICLE

**Fig. 4. Synaptic dynamics. (A)** Depressing, facilitating, and pseudolinear excitatory connections (top to bottom) in 50-Hz train; gray/colored dots are individual PSP amplitudes; black traces are average PSP per pulse. Scatter points for pulse 1 (resting-state adjusted coefficient of variation [aCV]) and pulse 8 (STP-induced aCV) are colored according to the color scale in (D). **(B)** STP matrix. **(C)** Recovery (at 250 ms) matrix. **(D)** Resting-state aCV matrix. All matrices are colored by the median (text in each element) with the saturation scaled by the SE. **(E)** Summary plots for paired-pulse ratio. STP induction ratio (average first pulse amplitude/average of the sixth to eighth pulse amplitude) normalized by the 90th percentile, resting-state aCV, and induced state aCV (top to bottom). Each dot corresponds to the average response from one connection. **(F)** Train-induced STP (top) at four different frequencies (10, 20, 50, and 100 Hz) for each of the elements in (E) (colors maintained). Each dot is the grand average of all connections in the element. For L5 ET→L5 ET, the blue shading highlights the 95% CI as an example. Lower plot shows recovery from STP at six different delays (125, 250, 500, 1000, 2000, and 4000 ms) in a similar manner to the plot above.

not well described. The variable delays that we imposed between the induction and recovery pulses (Fig. 1B, multipatch experiment) of our 50-Hz stimulus showed that at our earliest time point (125 ms), connections were still in their STP-induced state, but by 4 s, they were largely recovered.

Excitatory dynamics were strongly aligned with postsynaptic cell class (48) and further

refined by layer in the case of excitatory targets and by subclass of inhibitory targets. Recurrent excitatory connections were largely depressing (Fig. 4, B and E, L5ET→L5ET) and showed increasing depression with stimulus frequency (Fig. 4F). Recurrent excitatory connections occupied a range of recovery and variability profiles that varied with layer. Superficial layers (e.g., L2/3→L2/3) tended to

recover more quickly (Fig. 4C) and showed a higher degree of variability (Fig. 4D). E→I dynamics depend on the subclass of the postsynaptic target. Excitatory connections to Sst cells were strongly facilitating, consistent with previous reports (16, 47). E→Sst connections were also highly variable in the resting state, likely because of a high initial failure rate (Fig. 4A, middle, pulse 1, Fig. 4E); however,

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

6 of 13

Downloaded from https://www.science.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

strongly facilitating synapses often become more reliable in the induced state (Fig. 4A, middle, pulse 8, Fig. 4E). Excitatory connections onto Pvalb (Fig. 4E, L4 pyr→Pvalb) were largely depressing on average, although a subset of connections in L2/3 showed pseudolinear STP that was similar to in vivo measurements in somatosensory cortex (45). Whereas these patterns of excitatory dynamics were apparent on average, multiple measurements of synaptic dynamics showed high heterogeneity from pair to pair within a connection type (Fig. 4E).

Dynamics of inhibitory connections showed patterns more related to the subclass identity of the inhibitory presynaptic cell. Pvalb connections onto other subclasses, both excitatory (Fig. 4E, Pvalb→L6 pyr) and inhibitory, were strongly depressing (Fig. 4B) and were still depressed at our earliest recovery time point (Fig. 4C). Depressing Pvalb connections showed high reliability at the beginning of a stimulus train (resting-state log variability = -1.01; Fig. 4D) and became more variable in their STP-induced state (log variability = -0.27, MW $P = 3.65 \times 10^{-31}$; Fig. 4E). Connections from Sst and Vip cells had a mixture of dynamics from facilitation to moderate depression (Fig. 4E, Sst→L5 IT). These connections were also faster to recover from STP, particularly Sst connections, which tended to over-recover at the shortest interval (Fig. 4C).

Synaptic interactions between Sst and Vip were an exception to the trends highlighted above. Whereas most inhibitory connections were depressing, Sst→Vip (7) showed the highest degree of facilitation in our dataset (0.28 [0.0, 0.44]). The reciprocal Vip→Sst connection was weakly facilitating, as were recurrent Vip connections. These three connection types also over-recovered on short time scales (Fig. 4C) and took many seconds to fully recover (Fig. 4F). Given the facilitating nature of these connections, it is interesting that they had only a moderate degree of variance (resting-state log variability = -0.03; Fig. 4D) compared with other facilitating connections such as E→Sst (resting-state log variability = 1.66; MW $P = -3.0 \times 10^{-8}$).

### Human intralaminar connectivity

As a complement to the mouse visual cortex, our dataset includes synaptic physiology from human temporal cortex. Although our sampling of human connections covered all cortical layers, our analysis focused on the supragranular layers, which are expanded in anthropoid primate cortex (49). Deep L3 cells have distinct electrophysiology, morphology, and gene expression (including genes involved in connectivity and synaptic signaling), and many of these properties vary continuously with depth between L2 and L3b (50). Dense sampling of L2/3 allowed us to

define L2, L3a, and L3b pyramidal subclasses and demonstrate that these principles of cellular diversity have correlates in synaptic physiology. These subclasses showed distinct synaptic properties, including unique polysynaptic connections from L2 cells and STP that closely follow the continuous variability between L2/3 subclasses.

Distance dependence of connections was modeled and adjusted as with mouse connections but without distinguishing connections by cell class (most connections were recurrent excitatory). Connection probability was estimated to fall off with distance at a lateral spread ($\sigma$) of 130 $\mu$m (Gaussian model fit; fig. S5A), slightly larger than the comparable value in mouse (125 $\mu$m for within-class connections), reflecting that whereas cortical expansion is accompanied by the scaling of neuronal morphology, much of this scaling is axial rather than lateral. Examining the connectivity between subclasses (Fig. 5A), we tested for signatures of functional segregation within supragranular layers, finding a strong bias for recurrent connections over cross-connections between L3a and L3b (Fisher's $P = 5.6 \times 10^{-3}$) and a bias for connections from L2 to L3a over L3b (Fisher's $P = 3.4 \times 10^{-3}$).

Recurrent connectivity within human L4 ($p_{\max} = 0.01$ [0.0, 0.04] 95% CI, 1/145) was significantly lower than observed in the mouse ($p_{\max} = 0.22$ [0.16, 0.28], 44/452; Fisher's $P = 1.4 \times 10^{-4}$). This contrast could be related to age (37), species, or brain area. Furthermore, we observed other connections involving L4 pyramidal cells (e.g., I→E and E→I), suggesting that low excitatory recurrence was not a technical limitation of our dataset but rather reveals a unique property of the human L4 circuit.

### Human synaptic properties

The strength, kinetics, and STP properties of the human connections showed moderate differences across layers and large differences by cell class, largely resembling observations in mouse (Fig. 5B). Recurrent excitatory connections in human cortex (E→E) have longer latency than those with a pre- or postsynaptic inhibitory cell (E→E median 1.73 versus I→E 1.04 ms, KS test $P = 2.2 \times 10^{-5}$; E→I 1.34 ms, $P = 2.4 \times 10^{-3}$), and PSP rise times were faster for E→I than for E→E connections (2.42 versus 4.10 ms, $P = 2.0 \times 10^{-8}$) but slower for I→E connections (5.72 ms, $P = 5.3 \times 10^{-3}$, $P = 4.9 \times 10^{-6}$ versus E→E, E→I). We observe some differences between L2 and L3, including presynaptic L3 cells forming more depressing connections than L2 cells (STP ratios -0.39 versus -0.15, $P = 9.2 \times 10^{-3}$). E→I connections were uniformly depressing, consistent with the identification of those inhibitory cells as fast-spiking Pvalb cells.

For certain properties, we did observe contrasts between human and mouse connections.

The amplitudes of L2/3 excitatory PSPs were generally larger in human. For E→E connections, the contrast was stronger for the 90th percentile response (0.51 versus 0.30 mV, $P = 0.019$), whereas for E→I connections, the contrast was stronger for resting-state response (0.91 versus 0.26 mV, $P = 2.2 \times 10^{-3}$) because of significantly more depressing E→I connections in human (-0.51 versus 0.05, $P = 7.1 \times 10^{-6}$). We also observed faster recovery from STP in human than in mouse excitatory connections, with most fully recovered at 500 ms (Fig. 5C). This contrast has been previously noted in recurrent L2/3 excitatory connections (12), and our observations suggest that it also holds for L2/3 E→I connections.

### Human polysynaptic events

We found short-latency (~3 ms) inhibition (Fig. 5D) after stimulation of excitatory cells in human cortex, indicating activation of polysynaptic connections. Plotting latency versus PSP amplitude of human connections (Fig. 5D) revealed a clear boundary where responses with a latency of $\geq 3$ ms, evoked from a confirmed pyramidal cell, were almost exclusively inhibitory (median latency of IPSPs: 4.19 [3.79, 4.81 ms], compared with monosynaptic EPSPs, which had a latency of <3 ms (median EPSP latency: 1.71 [1.44, 2.22] ms). This potential disynaptic inhibition (dIPSPs) originated in L2 and projected to other L2 ($n = 8$) or L3 ($n = 11$) pyramidal cells, with just one ascending polysynaptic response originating in L3. This directionality was consistent with the directionality of monosynaptic excitation across layers 2 and 3 and supported by higher connection rates from L2 pyramidal cells to inhibitory cells.

### Variation with depth in human L2 and L3

Another property of supragranular neurons observed in human, but not in mouse, is strong depth-driven variability of intrinsic electrophysiological properties (50, 51). Visualizing the electrophysiology feature space by a UMAP projection of 27 electrophysiology features (see the materials and methods), we found that L4-type cells (high input resistance) are situated around the perimeter, indicating distinct properties from L2/3-type cells (Fig. 5E). For the L2/3-type cells, projecting a normalized layer depth coordinate (relative to the L2+L3 thickness) onto this space showed a mostly smooth gradient of electrophysiological properties with layer depth, also verified by direct examination of depth correlations with sag and AP upstroke/downstroke ratio ($P = 7.0 \times 10^{-6}$, $3.3 \times 10^{-5}$) (Fig. 5F). This correlation was not found in mouse L2/3 cells ($P > 0.06$ for both).

STP metrics revealed a similar linear variation with layer depth in connections from L2/3 pyramidal cells across both E→E ($n = 86$ human, $n = 15$ mouse) and E→I ($n = 20$ human,

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

7 of 13

Downloaded from https://www.sciencemag.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

Fig. 5. Human data. (A) Connection probability measured from human cortical tissue. Inhibitory cells are identified by morphology as asphyx or sparsely spiny cells, grouped across layer. (B) Kinetics, strength, and dynamics matrices are organized by layer for excitatory cells, with inhibitory cells grouped across layer. Each element is colored by the median (text in each element) with the saturation scaled to the SE. Two or more pairs were required to fill in an element. Latency, rise tau, and resting-state amplitude are quantified from fits of the average PSP response. (C) Train-induced STP (top) across frequencies for a subset of connection types. Each dot is the median of all connections in the element, with shading for the 95% CI (bootstrapped) shown for a single example connection type. Recovery from STP at different delays is shown in the lower plot. (D) Example polysynaptic circuit from one experiment in which cell 1 forms a short-latency (~2 ms) monosynaptic excitatory connection to cell 2 and delayed (~4 ms) polysynaptic inhibitory connections to cells 3 and 4 (all cells confirmed morphologically spiny). Dashed lines indicate (from left to right) time of presynaptic spike and PSP onset. Polysynaptic connections from L2/3 pyramidal cells were inferred by response latency >3 ms versus PSP amplitude. (E) Structure of intrinsic electrophysiology feature space. UMAP projection colored by cell subclass (left) and by depth of L2/3-type excitatory cells (right). (F) Variation in L2/3 intrinsic properties is strongly correlated with depth in human but not mouse. Example traces show superficial

A

B

C

D

E

F

G

and deep human cells [top, colors by L2/3 depth as in (E)]. Left is the phase plane representation of the first spike in a depolarizing step response. Right is the sag in response to hyperpolarization. Bottom is the regression of corresponding electrophysiology features versus depth by species, with bootstrapped 95% CI. (G) STP of L2/3 excitatory connections is structured by depth in human and not mouse. Top shows the PSP responses to spike trains for example cells from (F). The larger response to the first spike is quantified by paired-pulse STP, plotted below in relation to presynaptic depth (left) and AP up/down ratio (right) (postsynaptic relationships are shown in fig. S5E).

n = 72 mouse) connections. The paired-pulse STP showed a strong linear relationship with depth in the human data (P = 5.4 × 10⁻⁴), varying from weak facilitation for the most superficial cells (0.02 ± 0.03, mean ± SEM) to depression for the deepest (−0.32 ± 0.09). A strong correlation was also found with the AP upstroke/downstroke ratio of the presynaptic cell only (P = 1.2 × 10⁻⁴ versus P > 0.3 for post-

synaptic) (Fig. 5G). No corresponding trends were found in the mouse data (P > 0.2 for all regressions). Although lower sampling of L2/3 pyramidal cells in the mouse may contribute, regression coefficients for STP against the AP upstroke/downstroke ratio show a strong contrast (human −0.16 [−0.24, −0.08] CI; mouse 0.05 [−0.03, 0.13]), suggesting that there are real differences between these datasets in the

factors contributing to STP variability, whether explained by species, brain area, or other factors.

### Modeling STP of mouse and human connections

The STP metrics introduced above were chosen for their ease of interpretation but have several drawbacks: They are sensitive to noise, they require successfully repeated stimuli that are not available for all connections, and they are

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

8 of 13

Downloaded from https://www.science.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

difficult to use in a biophysical modeling context. Ideally, we would like a description that incorporates all of the data available for each connection and can predict synapse behavior in response to any arbitrary stimulus. We developed a generative model of stochastic vesicle release and STP with several adjustable parameters (see the materials and methods) and investigated which combinations of parameters could best explain the responses recorded for each connection. In this way, we captured and described more of the dynamic behavior of each connection with a small number of parameters. Our approach is similar to other recently developed models (52, 53) in that it does not depend on any particular stimulus structure (aside from having a diversity of interspike intervals). Accordingly, the model results for each connection make use of all presynaptic spikes and are robust to spike failures and early experiment termination.

Model performance was evaluated by using the maximum likelihood parameter set for each connection to simulate experimental data. These simulated data were then used to generate the same STP metrics that were previously collected from synaptic data. Both resting-state PSP amplitude and 90th percentile amplitude were almost perfectly correlated between recorded and simulated data (fig. S6, A and B), indicating that the model does exceptionally well at capturing synaptic strength. STP and variability (fig. S6, C to F) were also strongly correlated but with more scatter relative to strength metrics. STP measured from the second pulse in 50-Hz trains was only half as large as that measured in synaptic responses, indicating that the model as parameterized was not able to fully capture STP on this time scale. Release probability was more highly correlated with variance than with the number of release sites, suggesting that synapses may control variability primarily through their release probability.

Most discussions of short-term depression in the cortical literature begin with the assumption that depression is caused by the depletion of vesicles from the readily-releasable pool. However, recent evidence suggests that calcium channel inactivation may be a more prominent mechanism in cortical depression (54). We ran the model on two separate parameter spaces, one that uses vesicle depletion and another that uses a release-independent depression mechanism. In most cases, the model maximum likelihood value was found in the release-independent parameter space (fig. S7E, left). Release-dependent depression mechanisms should result in negative correlation between consecutive PSP amplitudes; however, we found little evidence for such negative correlations in our mouse data. Likewise, we found little relationship between paired correlation values and the model preference

for release-dependent mechanisms (fig. S7E, right). These results are consistent with the proposal that cortical synapses in mouse use release-independent depression mechanisms and that vesicle depletion plays a relatively minor role in depression. By comparison, our data from human connections had a modest preference for negative correlation between paired event amplitudes.

Organization in mouse and human synaptic dynamics

With a large dataset describing synaptic properties of connections, it becomes possible to ask what patterns emerge from the data. What synaptic features correlate with one another, and do connections naturally split into clusters based on these features or do they form a continuum? What aspects of the synaptic feature space are driven by presynaptic versus postsynaptic cell type? Excitatory connections onto Pvalb and Sst cells have distinctly different dynamics, suggesting a general rule that excitatory dynamics depend primarily on the postsynaptic cell type (16, 55). By contrast, inhibitory dynamics depend mainly on the presynaptic type, particularly when comparing Pvalb with Sst (47, 49). Although our data often followed these rules, we also found exceptions and suggest some refinements.

We used the stochastic release model described above to characterize 1196 connections (1035 mouse, 161 human). Although it would have been possible to use the maximum likelihood parameter set to describe each connection, parameters derived in this way are sensitive to experimental noise (56). Thus, we generated a more robust representation of synaptic behavior by exhaustively searching a large parameter space for each connection and then used sparse principal components analysis followed by UMAP dimensionality reduction to summarize these results (Fig. 6). This analysis groups connections based on the similarity of their model results; therefore, it has access to any synaptic strength and dynamical properties that the model could capture but does not have access to information about kinetics, cell subclass, or other cell properties. Connections in this analysis formed clear excitatory and inhibitory clusters (Fig. 6A), with a continuum of synaptic properties within each cluster. Perhaps the most prominent feature of this organization was that excitatory connections were strongly differentiated by the postsynaptic E/I class (66% classifier accuracy gain relative to shuffled; see the materials and methods), whereas inhibitory connection properties were mostly independent of postsynaptic cell class (7% gain).

Several synaptic properties were found to have a clear relationship to the UMAP organization (Fig. 6B). Which of these properties best explains the separation of excitatory synapses

by postsynaptic class? We generated a list of features from the maximum likelihood models and ranked these by the strength of their relationship to cell subclasses (see the materials and methods). When distinguishing between E→E and E→I connections, features that describe quantal release and synapse variance made up eight of the top 10 features (table S2). By contrast, STP parameters had relatively low importance in this classification, and the effect of synaptic strength was negligible. This relationship was also apparent when classifying excitatory connections among excitatory, Pvalb, Sst, and Vip postsynaptic subclasses (48% gain; Fig. 6B). Again, this relationship appears to be driven by quantal release parameters and variability metrics, whereas STP metrics were somewhat less important. Only minor relationships were found between excitatory synapse properties and postsynaptic excitatory subclasses (13% gain; Fig. 6C) or inhibitory synapse properties and postsynaptic subclasses (excitatory and inhibitory subclasses, 14% gain; excitatory subclasses, 5% gain; Fig. 6B). To simplify, excitatory synapses have low variability when the postsynaptic cell is excitatory, high variability when the postsynaptic cell is inhibitory, and highest variability onto Sst cells.

Inhibitory synapses, by contrast, did not follow this pattern. Instead, inhibitory synapse properties were more predictive of the presynaptic subclass (41% gain). The parameters that most contributed to this relationship include STP, variability, and to a lesser degree synaptic strength. Excitatory synapses were more weakly affected by presynaptic subclass (24% gain), but this relatively small effect appeared to be similarly driven by differences in STP and variability.

Human E→E synapses were strongly differentiated from mouse E→E synapses (43% gain). This was driven by a more diverse set of properties, including STP (especially recovery), variability, facilitation time constant, and strength (table S2, last column).

Discussion

By probing >20,000 possible connections across 28 mouse lines, we have explored a large fraction of the subclass-specific intralaminar connectivity in the mouse visual cortex. In addition, we sampled connections in the human cortex using similar methods for comparison. Past surveys near this scale have focused on the connectivity and strength of connections; a major advance provided by our study is the depth of characterization and analysis for each connection, in the context of transgenically identified cell subclasses and species.

Proposed standardized model of connectivity

The likelihood that two neurons are locally connected depends on multiple factors such as cell type, cortical region, species, and animal

Downloaded from https://www.sciencemag.org at Mayo Clinic on November 13, 2025

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

9 of 13

RESEARCH | RESEARCH ARTICLE

A

B

C

**Fig. 6. Dimensionality reduction on connection properties.** Relationships among connection properties and cell types revealed by dimensionality reduction. (A) All connections colored by postsynaptic E/I cell class. The UMAP output generates two clusters: inhibitory (left) and excitatory (right). (B) Four connection properties represented in reduced space, showing 90th percentile PSP amplitude (red, excitatory; blue, inhibitory); STP induced by 50 Hz trains (red, facilitating; blue, depressing), resting-state aCV during 50-Hz trains (purple, low variability; yellow, high variability), and the binomial CV derived from model parameters (release probability * number of release sites; purple, high CV; yellow, low CV). (C) Human and mouse connections colored by postsynaptic subclass. (D) Mouse connections colored by presynaptic subclass.

age. A longstanding goal has been to determine the governing principles of local circuit architecture. It is difficult to make direct comparisons among different studies because the observed rate of connectivity depends on several experimental details that are often inadequately reported or controlled for. Ideally, we would like a way to describe connectivity that allows more direct comparison between experiments regardless of their methodological differences. To facilitate such a comparison between connection subclasses in our own data, we developed a procedure for modeling connection probability as it relates to intersomatic

distance, axon truncation, cell depth, and signal detection power. With this model, we can estimate unbiased connection probabilities with confidence intervals that should be relatively robust to experimental bias. In principle, this approach is flexible enough to be replicated elsewhere in the field, and its adoption would substantially improve our ability to compare and reproduce results across studies.

**Conserved and canonical elements in the mouse intralaminar circuit**

As many previous studies have investigated the cortical circuit, a picture has emerged de-

scribing the relationships between the excitatory and inhibitory subclasses and their functional relevance. The details of this picture vary somewhat between descriptions, but a few key elements appear consistently, especially in the systems and theoretical neuroscience literature (Fig. 7A). Pvalb interneurons strongly inhibit nearby pyramidal cells and other Pvalb cells, Sst interneurons broadly inhibit nearby cells but avoid other Sst cells, and Vip cells selectively inhibit Sst cells and receive feedback excitation, forming a disinhibitory circuit (57, 58). We confirmed that each of these motifs is prominent across layers in the intralaminar cortical circuit.

Other circuit elements are equally prominent in our data but are more sparsely acknowledged in the literature. Many recent studies have focused on the importance of the Vip→Sst disinhibitory circuit. In the opposite direction, however, the connection from Sst to Vip has one of the highest connection probabilities, largest IPSP amplitudes, and strongest facilitation in our dataset. Although Sst→Vip connections have been described previously (4, 6, 7), they are often overlooked in consideration of the disinhibitory circuit. The synaptic features that we observed suggest an important functional ramification on the opposing Vip→Sst disinhibitory pathway (57, 59, 60).

Sst and Vip cells are often described as lacking recurrent connections (36, 37, 48, 61) despite some evidence to the contrary (6, 62). We confirmed that Sst and Vip do have the expected biases in their connectivity across all layers (Sst cells tend to avoid contacting other Sst cells, and Vip cells prefer to contact Sst cells). However, we also found sparse, recurrent connections within both interneuron populations approximately equal to the recurrent connectivity in excitatory populations. Furthermore, the strength of connections from Sst and Vip did not follow the same preferences, having roughly equal strength when connecting to preferred versus nonpreferred subclasses. Given their translaminar axon projection patterns, recurrent connections within these subclasses may be found more commonly across layer boundaries (63).

**Laminar variations on the cortical circuit**

Previous studies that sampled both L2/3 and L5 noted many similarities between the two layers (4, 6). Although we confirmed a consistent set of connectivity rules describing the intralaminar circuit, we also found variations on these rules that could support different modes of cortical function. Differences in intralaminar circuitry may contribute to laminar differences in receptive field properties (64) or visually mediated behaviors (58).

Layer 2/3 has strong interconnections between pyramidal and inhibitory cells (Fig. 7B).

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

10 of 13

Downloaded from https://www.science.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

**Fig. 7. Intralaminar circuit diagram.** The cortical intralayer circuit differs across layer and with activity. (A) Some commonly described elements of the intralaminar cortical circuit. Pvalb cells strongly inhibit pyramidal and other Pvalb cells. Sst cells provide broad inhibition, and Vip cells inhibit Sst cells to form a disinhibitory feedback pathway. (B and C) Circuit diagrams showing connections between major subclasses in mouse L2/3 (B) and L5 (C). The width of connecting lines roughly represents connection probability and PSP amplitude.

Connections that are prominent in each layer compared with the other are highlighted in orange, and green lines indicate connections that are less prevalent in that layer. For simplicity, connections from IT pyramidal to inhibitory cells in L5 (C) are omitted. (D) Two complementary circuits that activate at different times. Red connections are facilitating and will be stronger during sustained activity. Blue connections are depressing and are strongest during quiescent periods.

In deep layers, these connections are either absent or greatly reduced, and the relative sparsity of Vip cells in particular in deep layers should further enhance these differences (65, 66). Direct Vip inhibition of local pyramidal cells further complicates the prevailing view of the Vip→Sst disinhibitory pathway by, for example, allowing the possibility of feedback inhibition from higher cortical regions. Likewise, local excitatory inputs to Vip cells could be a source of feedforward disinhibition in layer 2/3.

Layer 5 excitatory subclasses differ in their visual responses and long-range projections, suggesting different functional roles in the circuit (67, 68). Accordingly, we saw differences in the intralaminar connections of L5 ET and IT neurons. ET pyramidal cells were generally more highly connected, receiving more local excitation and inhibition than IT cells. We confirmed a much higher rate of recurrent connections among L5 ET cells compared with IT cells and also found that connections between these two subclasses were unidirectional from IT→ET (Fig. 7C) (34). Layer 5 ET cells also received more frequent inhibition from Sst cells, as previously observed in the rat frontal cortex (35).

Could laminar differences in connectivity indicate cell type divisions within subclasses? Two recent studies investigated the correspondence among morphological, electrophysiological, and transcriptomic (MET) features of inhibitory neurons in primary visual cortex and motor cortex (63, 69). In visual cortex, different MET types had distinct patterns of local axonal innervation and dendritic morphologies (63), suggesting that their connec-

tivity will be different. MET types also exhibited layer localization, and thus some of the differences in connectivity that we observed as a function of layer may reflect differences in connectivity between different MET types.

#### **Dynamic flexibility in the cortical circuit**

Most studies in cortical synaptic physiology describe the circuit in its quiescent state. Ongoing activity in vivo, however, can transiently depress or facilitate synapses in a cell type-dependent way. In effect, the cortical circuit is modified from moment to moment based on the recent history of activity, so it may be misleading to conceptualize the cortex as a single, static circuit diagram. Although STP is variable across individual synapses, many subclass elements of the cortical circuit had a clear preference for either facilitation or depression (Fig. 7D). These patterns suggest an ability to dynamically switch between different network modalities in which intralaminar activity shifts from depressing interactions between pyramidal and Pvalb cells at the onset of activity to facilitating interactions among pyramidal, Sst, and Vip cells during sustained activity. Ultimately, activity patterns in vivo determine the extent of STP in each cell type, and the functional effects have yet to be determined.

#### **Excitatory cells receive weak local excitation compared with inhibitory cells**

Our previous study (13) observed low recurrent excitatory connectivity rates in all layers. We now find that recurrent excitatory connections have a set of features that distinguish them from excitatory inputs to inhibitory cells

and appear to limit their contribution to excitability. In addition to being sparse, they are relatively weak, get weaker with activity, and have slower PSP rise times. Slow PSP rise times limit excitability by raising the AP threshold (70). Slow rise times have been attributed to the PSPs being less synchronous at the multiple synapses that form E-E connections (71). Unlike the relatively high rate of intralaminar connectivity for E-I connections observed in all layers, recurrent excitatory connectivity is generally lower, and nearly absent in human L4 and in CT cells of mouse L6.

#### **Unidirectional disynaptic inhibition in human**

We observed disynaptic inhibition in human cortex between confirmed spiny pyramidal cells that is unidirectional, originating in L2 and targeting other L2 or L3 pyramidal cells. We did not see disynaptic inhibition in our mouse recordings, which may be due to the stronger excitation that we and others (9) have observed in human connections, particularly onto inhibitory cells. Disynaptic inhibition is often mediated by an interposed Sst cell (72, 73) because these cells have a low spiking threshold and receive facilitating inputs from excitatory cells. The latency of Sst-mediated disynaptic inhibition is often long (>100 ms), but we saw disynaptic IPSPs with a much shorter latency (3 to 6 ms). This suggests the disynaptic inhibition is driven by an intermediate fast-spiking Pvalb cell, which has been observed in humans and can be recruited by very large excitatory events (14). The unidirectional nature of this disynaptic inhibition from more superficial to deeper cortex further

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

11 of 13

Downloaded from https://www.sciencet.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

suggests a preferential routing of information by Pvalb cells in human cortex (14).

Connection types differ in variability

We used a new model of stochastic vesicle release and STP to estimate parameters that best describe each connection. A few governing principles emerge from this analysis. Cortical connection types could be organized into a two-dimensional feature space, with synaptic strength and variability forming two orthogonal gradients. Excitatory connections are grouped along the variability axis based on their postsynaptic subclass. These connections have lower variability when connecting to excitatory cells, higher variability onto Pvalb cells, and the highest variability onto Sst cells. By contrast, inhibitory connections partition more naturally based on their presynaptic subclass. This partition again separates connections based on variability, with Pvalb connections having overall lower variability than other inhibitory types. In most cases, variability correlates with STP, and the highest variance connections are more likely to be facilitating. The rules derived from these relationships are similar to rules described previously (47, 48, 74). However, we found synaptic variability to be a better predictor of cell subclass than STP or strength, particularly for excitatory connections.

Synaptic variability has been regarded as an undesirable consequence of signaling through metabolically expensive exocytosis. More recently, advances in machine learning that rely on stochasticity have supported the possibility that synaptic variability may offer computational benefits, such as a mechanism for regularization during learning (21). If that is the case, then it is further plausible that variability may be modulated by cell type and that these relationships are crucial features of cortical function. Indeed, our measurements of synaptic variability strongly differentiated cell type in a pattern that is largely, but not entirely, aligned with STP metrics, suggesting the possibility of cell type-specific tuning of variability. More broadly, analysis of our stochastic release model indicates that synaptic strength and variability form the two most significant parameters describing synapse behavior.

REFERENCES AND NOTES

1. A. M. Thomson, A. P. Bannister, Interlaminar connections in the neocortex. Cereb. Cortex 13, 5–14 (2003), doi: 10.1093/cercor/13.1.5; pmid: 12466210
2. H. Markram et al., Interneurons of the neocortical inhibitory system. Nat. Rev. Neurosci. 5, 793–807 (2004), doi: 10.1038/nrn1519; pmid: 15378039
3. R. J. Douglas, K. A. C. Martin, Neuronal circuits of the neocortex. Annu. Rev. Neurosci. 27, 419–451 (2004), pmid: 15217339
4. C. K. Pfeffer, M. Xue, M. He, Z. J. Huang, M. Scanziani, Inhibition of inhibition in visual cortex: The logic of connections between molecularly distinct interneurons. Nat. Neurosci. 16, 1068–1076 (2013), doi: 10.1038/nrn3446; pmid: 23817549

5. K. D. Harris, G. M. G. Shepherd, The neocortical circuit: Themes and variations. Nat. Neurosci. 18, 170–181 (2015), doi: 10.1038/nrn3917; pmid: 25622573
6. X. Jiang et al., Principles of connectivity among morphologically defined cell types in adult neocortex. Science 350, aac9462 (2015), doi: 10.1126/science.aac9462; pmid: 26612957
7. M. M. Karnani et al., Cooperative subnetworks of molecularly similar interneurons in mouse neocortex. Neuron 90, 86–100 (2016), doi: 10.1016/j.neuron.2016.02.037; pmid: 27021171
8. R. Tremblay, S. Lee, B. Rudy, GABAergic interneurons in the neocortex: From cellular properties to circuits. Neuron 91, 260–292 (2016), doi: 10.1016/j.neuron.2016.06.033; pmid: 27477017
9. G. Molnár et al., Complex events initiated by individual spikes in the human cerebral cortex. PLOS Biol. 6, e222 (2008), doi: 10.1371/journal.pbio.0060222; pmid: 18767905
10. G. Molnár et al., Human pyramidal to interneuron synapses are mediated by multi-vesicular release and multiple docked vesicles. eLife 8, e18167 (2016), doi: 10.7554/eLife.18167; pmid: 27536876
11. Y. Peng et al., High-throughput microcircuit analysis of individual human brains through next-generation multineuron patch-clamp. eLife 8, e48178 (2019), doi: 10.7554/eLife.48178; pmid: 31742558
12. G. Testa-Silva et al., High bandwidth synaptic communication and frequency tracking in human neocortex. PLOS Biol. 12, e1002007 (2014), doi: 10.1371/journal.pbio.1002007; pmid: 25422947
13. S. C. Seeman et al., Sparse recurrent excitatory connectivity in the microcircuit of the adult mouse and human cortex. eLife 7, e37349 (2018), doi: 10.7554/eLife.37349; pmid: 30256194
14. Y. Szeged et al., High-precision fast-spiking basket cell discharges during complex events in the human neocortex. eNeuro 4, ENEURO.0260-17.2017 (2017), doi: 10.1523/eneuro.0260-17.2017; pmid: 29034319
15. Y. N. Billeh et al., Systematic integration of structural and functional data into multi-scale models of mouse primary visual cortex. Neuron 106, 388–403.e18 (2020), doi: 10.1016/j.neuron.2020.01.040; pmid: 32142648
16. A. Reyes et al., Target-cell-specific facilitation and depression in neocortical circuits. Nat. Neurosci. 1, 279–285 (1998), doi: 10.1038/nrn1092; pmid: 10195160
17. A. V. Blackman, T. Abrahamsson, R. P. Costa, T. Lalanne, P. J. Sjöström, Target-cell-specific short-term plasticity in local circuits. Front. Synaptic Neurosci. 5, 11 (2013), doi: 10.3389/fnsyn.2013.00011; pmid: 24367330
18. R. S. Larsen, P. J. Sjöström, Synapse-type-specific plasticity in local circuits. Curr. Opin. Neurobiol. 35, 127–135 (2015), doi: 10.1016/j.conb.2015.08.001; pmid: 25310110
19. S. Lefort, C. C. H. Petersen, Layer-dependent short-term synaptic plasticity between excitatory neurons in the C2 barrel column of mouse primary somatosensory cortex. Cereb. Cortex 27, 3869–3876 (2017), doi: 10.1093/cercor/bhx094; pmid: 28444185
20. H. Markram et al., Reconstruction and simulation of neocortical microcircuit. Cell 163, 456–492 (2015), doi: 10.1016/j.cell.2015.09.029; pmid: 26451489
21. M. Llera-Montero, J. Sacramento, R. P. Costa, Computational roles of plastic probabilistic synapses. Curr. Opin. Neurobiol. 54, 90–97 (2019), doi: 10.1016/j.conb.2018.09.002; pmid: 30308457
22. L. F. Abbott, W. G. Regehr, Synaptic computation. Nature 431, 796–803 (2004), doi: 10.1038/nature03010; pmid: 15483601
23. D. Burnham, E. Shea-Brown, S. Mihalas, Learning to predict in networks with heterogeneous and dynamic synapses. bioRxiv 444107 [Preprint] (2021), doi: 10.1101/2021.05.18.444107
24. T. L. Dagle et al., A suite of transgenic driver and reporter mouse lines with enhanced brain-cell-type targeting and functionality. Cell 174, 465–480.e22 (2018), doi: 10.1016/j.cell.2018.06.035; pmid: 30007418
25. R. B. Levy, A. D. Reyes, Spatial profile of excitatory and inhibitory synaptic connectivity in mouse primary auditory cortex. J. Neurosci. 32, 5609–5619 (2012), doi: 10.1523/JNEUROSCI.5158-11.2012; pmid: 22514322
26. R. Perin, T. K. Berger, H. Markram, A synaptic organizing principle for cortical neuronal groups. Proc. Natl. Acad. Sci. U.S.A. 108, 5419–5424 (2011), doi: 10.1073/pnas.1016051108; pmid: 21383177
27. C. Holmgren, T. Harkany, B. Svennerfors, Y. Zilberter, Pyramidal cell communication within local networks in layer 2/3 of rat neocortex. J. Physiol. 551, 139–153 (2003), doi: 10.1113/jphysiol.2003.044784; pmid: 12813147

28. L. F. Rossi, K. D. Harris, M. Carandini, Spatial connectivity matches direction selectivity in visual cortex. Nature 588, 648–652 (2020), doi: 10.1038/s41586-020-2894-4; pmid: 33177719
29. S. Song, P. J. Sjöström, M. Riegl, S. Nelson, D. B. Chiklovski, Highly nonrandom features of synaptic connectivity in local cortical circuits. PLOS Biol. 3, e68 (2005), doi: 10.1371/journal.pbio.0030068; pmid: 15737062
30. F. Z. Hoffmann, J. Triesch, Nonrandom network connectivity comes in pairs. Netw. Neurosci. 1, 31–41 (2017), doi: 10.1162/NETN_a_00004; pmid: 29601066
31. A. Stepenyants, L. M. Martinez, A. S. Fereczko, Z. F. Kisvárday, The fractions of short- and long-range connections in the visual cortex. Proc. Natl. Acad. Sci. U.S.A. 106, 3555–3560 (2009), doi: 10.1073/pnas.0803930106; pmid: 19221032
32. S. Lefort, C. Tomm, J. C. Floyd, S. C. Petersen, The excitatory neuronal network of the C2 barrel column in mouse primary somatosensory cortex. Neuron 61, 301–316 (2009), doi: 10.1016/j.neuron.2008.12.020; pmid: 19186171
33. J.-S. Jouhanneau, J. Krenikov, A. L. Donn, J. F. Poulet, In vivo monosynaptic excitatory transmission between layer 2 cortical pyramidal neurons. Cell Rep. 13, 2098–2106 (2015), doi: 10.1016/j.celrep.2015.11.011; pmid: 26670044
34. S. P. Brown, S. Hestrin, Intracortical circuits of pyramidal neurons reflect their long-range axonal targets. Nature 457, 1133–1136 (2009), doi: 10.1038/nature07658; pmid: 19151698
35. M. Morishima, K. Kobayashi, S. Kato, K. Kobayashi, Y. Kawaguchi, Segregated excitatory–inhibitory recurrent subnetworks in layer 5 of the rat frontal cortex. Cereb. Cortex 27, 5846–5857 (2017), doi: 10.1093/cercor/bhx276; pmid: 29045559
36. A. Naka et al., Complementary networks of cortical somatostatin interneurons enforce layer specific control. eLife 8, e43696 (2019), doi: 10.7554/eLife.43696; pmid: 30883329
37. F. Scala et al., Layer 4 of mouse neocortex differs in cell types and circuit organization between sensory areas. Nat. Commun. 10, 4174 (2019), doi: 10.1038/s41467-019-12058-z; pmid: 31519874
38. H. Hu, J. Z. Cavendish, A. Agmon, Not all that differs in gold: Off-target recombination in the somatostatin-REG-CoA mouse line labels a subset of fast-spiking interneurons. Front. Neural Circuits 7, 195 (2013), doi: 10.3389/fncr.2013.00195; pmid: 24339803
39. B. R. Lee et al., Scaled, high fidelity electrophysiological, morphological, and transcriptomic cell characterization. eLife 10, e65482 (2021), doi: 10.7554/eLife.65482; pmid: 34387544
40. M. Galarreta, S. Hestrin, A network of fast-spiking cells in the neocortex connected by electrical synapses. Nature 402, 72–75 (1999), doi: 10.1038/47029; pmid: 10573418
41. M. Galarreta, S. Hestrin, Electrical and chemical synapses among parvalbumin fast-spiking GABAergic interneurons in adult mouse neocortex. Proc. Natl. Acad. Sci. U.S.A. 99, 12438–12443 (2002), doi: 10.1073/pnas.192159599; pmid: 12213962
42. P. Alcami, A. E. Pereda, Beyond plasticity: The dynamic impact of electrical synapses on neural circuits. Nat. Rev. Neurosci. 20, 253–271 (2019), doi: 10.1038/s41583-019-0133-5; pmid: 30824857
43. F. Walker et al., Parvalbumin- and vasoactive intestinal polypeptide-expressing neocortical interneurons impose differential inhibition on Martinotti cells. Nat. Commun. 7, 13664 (2016), doi: 10.1038/ncomms13664; pmid: 27897179
44. A. M. Thomson, D. C. West, Y. Wang, A. P. Bannister, Synaptic connections and small circuits involving excitatory and inhibitory neurons in layers 2-5 of adult rat and cat neocortex: Triple intracellular recordings and biocytin labelling in vitro. Cereb. Cortex 12, 936–953 (2002), doi: 10.1093/cercor/12.9.936; pmid: 12183393
45. A. Pala, C. C. H. Petersen, In vivo measurement of cell-type-specific synaptic connectivity and synaptic transmission in layer 2/3 mouse barrel cortex. Neuron 85, 68–75 (2015), doi: 10.1016/j.neuron.2014.11.025; pmid: 25543458
46. A. M. Thomson, Activity-dependent properties of synaptic transmission at two classes of connections made by rat neocortical pyramidal axons in vitro. J. Physiol. 502, 131–147 (1997), doi: 10.1111/j.1469-7793.1997.1131.x; pmid: 9134202
47. M. Beerlein, J. R. Gibson, B. W. Connors, Two dynamically distinct inhibitory networks in layer 4 of the neocortex. J. Neurophysiol. 90, 2987–3000 (2003), doi: 10.1152/jn.00283.2003; pmid: 12815025
48. Y. Ma, H. Hu, A. Agmon, Short-term plasticity of unitary inhibitory-to-inhibitory synapses depends on the presynaptic interneuron subtype. J. Neurosci. 32, 983–988 (2012), doi: 10.1523/JNEUROSCI.5007-11.2012; pmid: 22262896

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

12 of 13

Downloaded from https://www.scienc.org at Mayo Clinic on November 13, 2025

RESEARCH | RESEARCH ARTICLE

49. P. Balaram, J. H. Kaas, Towards a unified scheme of cortical lamination for primary visual cortex across primates: Insights from NeuN and VGLUT2 immunoreactivity. Front. Neuroanat. 8, 81 (2014), doi: 10.3389/fnana.2014.00081; pmid: 2517277
50. J. Berg et al., Human neocortical expansion involves glutamatergic neuron diversification. Nature 598, 151–158 (2021), doi: 10.1038/s41586-021-03813-8; pmid: 34616067
51. B. E. Kalmbach et al., n-channels contribute to divergent intrinsic membrane properties of supragranular pyramidal neurons in human versus mouse cerebral cortex. Neuron 100, 1194–1208.e5 (2018), doi: 10.1016/j.neuron.2018.10.012; pmid: 30392798
52. A. D. Bird, M. J. Wall, M. J. E. Richardson, Bayesian inference of synaptic quantal parameters from correlated vesicle release. Front. Comput. Neurosci. 10, 116 (2016), doi: 10.3389/fncom.2016.00116; pmid: 27932970
53. A. Barri, Y. Wang, D. Hansel, G. Mongillo, Quantifying repetitive transmission at chemical synapses: A generative-model approach. eNeuro 3, ENEURO.0113-15.2016 (2016), doi: 10.1523/ENEURO.0113-15.2016; pmid: 27200414
54. E. Nanou, W. A. Catterall, Calcium channels, synaptic plasticity, and neuropsychiatric disease. Neuron 98, 466–481 (2018), doi: 10.1016/j.neuron.2018.03.017; pmid: 29723500
55. H. Markram, Y. Wang, M. Tsodyks, Differential signaling via the same axon of neocortical pyramidal neurons. Proc. Natl. Acad. Sci. U.S.A. 95, 5323–5328 (1998), doi: 10.1073/pnas.95.5323; pmid: 9560274
56. O. Bykowska et al., Model-based inference of synaptic transmission. Front. Synaptic Neurosci. 11, 21 (2019), doi: 10.3389/fnsyn.2019.00021; pmid: 31481887
57. A. J. Keller et al., A Disinhibitory circuit for contextual modulation in primary visual cortex. Neuron 108, 1181–1193.e8 (2020), doi: 10.1016/j.neuron.2020.11.013; pmid: 33301712
58. Y. Fu et al., A cortical circuit for gain control by behavioral state. Cell 156, 1139–1152 (2014), doi: 10.1016/j.cell.2014.01.050; pmid: 24630718
59. D. J. Millman et al., VIP interneurons in mouse primary visual cortex selectively enhance responses to weak but specific stimuli. eLife 9, e55130 (2020), doi: 10.7554/eLife.55130; pmid: 33108272
60. M. Dipoppa et al., Vision and locomotion shape the interactions between neuron types in mouse visual cortex. Neuron 98, 602–615.e8 (2018), doi: 10.1016/j.neuron.2018.03.037; pmid: 29656873
61. H. Hu, Y. Ma, A. Agmon, Submillisecond firing synchrony between different subtypes of cortical interneurons connected chemically but not electrically. J. Neurosci. 31, 3351–3361 (2011), doi: 10.1523/JNEUROSCI.4881-10.2011; pmid: 21368047
62. J. R. Gibson, M. Beierlein, B. W. Connors, Two networks of electrically coupled inhibitory neurons in neocortex. Nature 402, 75–79 (1999), doi: 10.1038/47035; pmid: 10573419
63. N. W. Gouwens et al., Integrated morphoelectric and transcriptomic classification of cortical GABAergic cells. Cell 183, 935–953.e19 (2020), doi: 10.1016/j.cell.2020.09.057; pmid: 33186530
64. C. M. Neill, M. P. Stryker, Highly selective receptive fields in mouse visual cortex. J. Neurosci. 28, 7520–7536 (2008), doi: 10.1523/JNEUROSCI.0623-08.2008; pmid: 18650330
65. Z. Almási, C. Dávid, M. Witte, J. F. Staiger, Distribution patterns of three molecularly defined classes of GABAergic neurons across columnar compartments in mouse barrel cortex. Front. Neuroanat. 13, 45 (2019), doi: 10.3389/fnana.2019.00045; pmid: 31114486
66. Y. Gonchar, Q. Wang, A. Burkhalter, Multiple distinct subtypes of GABAergic neurons in mouse visual cortex identified by triple immunostaining. Front. Neuroanat. 1, 3 (2008), doi: 10.3389/fnana.2008.00043; pmid: 18958197
67. G. Lur, M. A. Vinck, L. Tang, J. A. Cardin, M. J. Higley, Projection-specific visual feature encoding by layer 5 cortical subnetworks. Cell Rep. 14, 2538–2545 (2016), doi: 10.1016/j.celrep.2016.02.050; pmid: 26972011
68. E. J. Kim, A. L. Juvinett, E. M. Kyubwa, M. W. Jacobs, E. M. Callaway, Three types of cortical layer 5 neurons that differ in brain-wide connectivity and function. Neuron 88, 1253–1267 (2015), doi: 10.1016/j.neuron.2015.11.002; pmid: 26671462
69. F. Scala et al., Phenotypic variation of transcriptomic cell types in mouse motor cortex. Nature 598, 144–150 (2021), doi: 10.1038/s41586-020-2907-3; pmid: 33184512
70. R. Azouz, C. M. Gray, Dynamic spike threshold reveals a mechanism for synaptic coincidence detection in cortical neurons in vivo. Proc. Natl. Acad. Sci. U.S.A. 97, 8110–8115 (2000), doi: 10.1073/pnas.130200797; pmid: 10859358
71. S. J. Barnes et al., Delayed and temporally imprecise neurotransmission in reorganizing cortical microcircuits. J. Neurosci. 35, 9024–9037 (2015), doi: 10.1523/JNEUROSCI.4583-14.2015; pmid: 26085628
72. J. Obermayer et al., Lateral inhibition by Martinotti interneurons is facilitated by cholinergic inputs in human and mouse neocortex. Nat. Commun. 9, 4101 (2018), doi: 10.1038/s41467-018-06628-w; pmid: 30291244
73. G. Silberberg, H. Markram, Disynaptic inhibition between neocortical pyramidal cells mediated by Martinotti cells. Neuron 53, 735–746 (2007), doi: 10.1016/j.neuron.2007.02.012; pmid: 17329212
74. A. Gupta, Y. Wang, H. Markram, Organizing principles for a diversity of GABAergic interneurons and synapses in the neocortex. Science 287, 273–278 (2000), doi: 10.1126/science.287.5451.273; pmid: 10634775

ACKNOWLEDGMENTS

We thank the founder of the Allen Institute, P. G. Allen, for his vision, encouragement, and support. Funding: This work was supported by the National Institutes of Health (grants R01DA036909 and RF1MH121274 to B.T., grant U19MH114830 to H.Z., and grant U01MH114812 to E.L. and T.C.). Author contributions: Conceptualization: T.J., C.K., H.Z.; Data curation: L.C., S.C.S., P.D., A.Ho., C.R., L.K., J.T., C.G., D.F., A.S., W.W., N.H., S.Sh., L.A., G.W., A.He., A.M., S.K., D.S., P.B., T.J.; Formal analysis: L.C., S.C.S., T.C., S.I., C.T., J.G.; Investigation: L.C., S.C.S., P.D., A.Ho., L.K., J.T., C.R., D.C., T.C., K.C., M.C., L.E., J.G., M., J.S., H.T., K.W., K.S., M.T., T.P., D.M., A.G., K.B., T.E., M.Ma., M.Mc., C.A.P., A.R., J.B., K.B., N.D., R.E., M.G., M.H., S.D.L., K.N., R.N., L.P., S.R.; Methodology: T.J., L.C.; Project administration: F.D., S.Su.; Resources: R.L., T.C., N.D., C.F., M.S., C.S., B.T., M.W., M.G., T.D., C.C., R.P.G., C.D.K., A.L.K., J.G.O., D.L.S.; Software: T.B., T.J., L.C.; Supervision: T.J., J.P., L.E., H.Z., E.L., G.M., M.Mc., N.D., R.N., S.So.; Visualization: L.C., S.C.S., T.C., S.I., C.G.; Writing – original draft: T.J., L.C., S.C.S., T.C., S.I., L.K., C.G.; Writing – review and editing: M.H.K., T.H., A.A. Competing interests: The authors declare no competing interests. Data and materials availability: All data and tools are available from our web portal at https://portal.bran-map.org/explore/connectivity/synaptic-physiology.

SUPPLEMENTARY MATERIALS

science.org/doi/10.1126/science.abj5861

Materials and Methods

Figs. S1 to S14

Tables S1 to S5

References (75–90)

MDAR Reproducibility Checklist

View/request a protocol for this paper from Bio-protocol.

28 May 2021; accepted 31 January 2022

10.1126/science.abj5861

Campagnola et al., Science 375, eabj5861 (2022) 11 March 2022

13 of 13

Downloaded from https://www.science.org at Mayo Clinic on November 13, 2025
