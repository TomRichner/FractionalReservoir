function cfg = sfaEI_mu7_grouped_figures_config()
% Figure-only replay: NO calls to analysis or legacy inline simulations.
cfg=struct();
cfg.preset_name='celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7';
cfg.run_mode='fast';
cfg.run_dir='data/sfaEI_mu7_fast';
cfg.fig_root='figs/sfaEI_mu7_fast_grouped';
cfg.visible_figures=false;
cfg.verbose='minimal';
root=fileparts(which('setup_paths'));
old_learning=fullfile(fileparts(root),'StochasticPlasticDynamicalSystemPaper','figs_pytorch','plot_for_Brian_seeds.png');
cfg.figures=grouped_figure_registry(cfg.preset_name,'figs/sfaEI_mu7_fast',old_learning);
end
