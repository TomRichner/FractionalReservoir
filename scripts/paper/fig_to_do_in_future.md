# These are user notes, not a to do list for agents.

- add more PSD figures comparing no adaptation, 1 TS, 3 TS

- fix the sensitivity plots (don't force the hists into range.  Make the darkest color dark gray, not black, fix the xlabel overlap.  Go wiht -50% to +50%

- median sensitivity plots need different colors

- memory capacity needs to be redone

- make the jacobian density plot use f(x) = log(1+log(1+x)) (doing now)

- make separate folders per run+config

- config.m files need to be independent

- time series example should show 1 synaptic output per neuron if all have same STD routes

- time series should show depression for for all if all have same STD routes

- time series should show a subset of neurons

- time series should show local lya

- local lya density histogram over time

- Eliminate SRNNModel2.m

- Eigenvalue density plot for different imbalanced networks or networks with E only SFA

- Verify precision of SRA1 vs ODE45 in noise-free (and maybe with-noise?) simulations and the similarity of LLE.  make a supplemental figure based on a 3 example time series plots with the base preset in the three conditions of the preset.  

- verify the benettin reshooting vs QR methods agree for SRNNCellTypePairs.m (it was done on the SRNNModel2.m).  this can also be done with the SRA1 vs ODE45 check.  

- only tau_a_E was swept for the tau sensitivity analysis.  it should be both E and I for runs in which both E and I have SFA.  We need to fix this.  Can this be done considering how the param space analysis class works? 

- could we go to a single cell type model (but with a dale's law weight matrix)?  this would reduce the std routes and possibly clean up the number of eignevalues in the jacobian.  results should be the same.  might make the tau_a sweep work without modification.  Is single cell type supported?  We ran into that during a refactor previously. 