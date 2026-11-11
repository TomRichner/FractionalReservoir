% FUTURE full analysis: run on the stronger computer. Not part of a replot.
% Fast is the default; set run_mode to 'medium' for the longer/statistical run.
setup_paths();
run_mode = 'fast';
cfg = mu7revised_config(run_mode);
run_dir = run_all_paper_analyses(cfg);
results = make_all_paper_figures(cfg);
