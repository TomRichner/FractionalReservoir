pp = PairedPulseMIAnalysis('n_networks', 5, 'note', 'test');
pp.model_defaults.T_range = [0, 100];  % Shorter for quick test
pp.run();
pp.plot();
pp.plot_histograms();
pp.plot_mi_imagesc();
pp.plot_mi_vs_lle();