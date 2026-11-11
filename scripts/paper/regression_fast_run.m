function [run_dir, results, minutes] = regression_fast_run(tag, opts)
% REGRESSION_FAST_RUN Run the regression bundle: every analysis stage and every figure, on the small network.
%
%   regression_fast_run()                       % tag 'reference'
%   [run_dir, results, minutes] = regression_fast_run('after_change')
%   regression_fast_run('abc', 'figures', false)   % analyses only
%
% Runs both entry points on regression_fast_config(tag):
%   run_all_paper_analyses  ->  data/regression_fast/<tag>
%   make_all_paper_figures  ->  figs/regression_fast/<tag>
%
% A FUNCTION, unlike the other *_run.m files, because it takes the tag and is
% called by test_regression_bundle. It can still be launched cold: setup_paths
% is its first line.
%
% THE RUN DIRECTORY MUST NOT EXIST OR MUST BE EMPTY (run_all_paper_analyses
% refuses otherwise). Use a new tag, or delete the old directory deliberately.
%
% Outputs: the run directory, the per-figure results of
% make_all_paper_figures (empty when figures are off), and the wall time in
% minutes as [analyses, figures].
%
% See also: regression_fast_config, regression_digest, test_regression_bundle
arguments
    tag (1,:) char = 'reference'
    opts.figures (1,1) logical = true
end
setup_paths();
cfg = regression_fast_config(tag);

t0 = tic;
run_dir = run_all_paper_analyses(cfg);
minutes = [toc(t0) / 60, 0];

results = [];
if opts.figures
    t1 = tic;
    results = make_all_paper_figures(cfg);
    minutes(2) = toc(t1) / 60;
end

fprintf('\n========================================================\n');
fprintf('REGRESSION BUNDLE COMPLETE (tag %s)\n', tag);
fprintf('  preset   : %s (%s)\n', cfg.preset_name, cfg.run_mode);
fprintf('  run_dir  : %s\n', run_dir);
fprintf('  analyses : %.1f min\n', minutes(1));
if opts.figures
    fprintf('  figures  : %d of %d succeeded, %.1f min\n', ...
        sum([results.ok]), numel(results), minutes(2));
end
fprintf('========================================================\n');
end
