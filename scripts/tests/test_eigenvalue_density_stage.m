%% test_eigenvalue_density_stage.m
% The opt-in stage 'eigenvalue_density' (FD-042): its theory, its plan, the
% stage on contract-conformant outlier_spectra files, and its five figures.
%
% Checks:
%   1. THE THEORY against sampled matrices (eigdens_check_theory, the five
%      synthetic cases of the exploration's edens_check_theory) at n = 200,
%      10 draws: sampled and predicted cumulative radial distributions
%      within 3.5 sampling errors (1/sqrt(n draws)); the column solver and
%      the block solver within 1e-3; Harris's closed form against the column
%      solver within 1e-3 in density relative to uniform. Plus the closed
%      forms: uniform gains give (r/R)^2 and edge_rel 1; half the gains zero
%      gives an atom of 1/2 and a uniform rest; the measures of a uniform
%      disk.
%   2. THE PLAN: a bundle naming the stage and outlier_spectra is planned
%      (folder, input folder, file, grids); refused (BadConfig) without
%      outlier_spectra, with an unknown option, or with an unknown run mode.
%   3. THE STAGE on contract-conformant files (make_synthetic_outlier_spectra:
%      an 80-neuron mu7chain network, reference and E-to-E x 1.5, three
%      conditions, two seeds, three states each, the two layouts of R.cal):
%      every index row becomes a run, in order, with the shapes the figures
%      read; the theoretical radius from the saved gains equals the one the
%      fixture computed (R_eff, 1e-9); the frozen-gain eigenvalue counts are
%      the spectrum's own; the loader refuses a missing input and a
%      mismatched network (errors, no fallback).
%   4. THEORY WITH THE SAVED GAINS against the measured cumulative radial
%      distribution, at the test size: the largest gap inside the bulk at a
%      single state is below 0.15 (one 80-eigenvalue spectrum: a
%      Kolmogorov-Smirnov distance of 1.36/sqrt(80) = 0.15 is the 5% level),
%      and the mean of the measured curves over the six states of a network
%      and condition is within 0.08 of the mean theory (480 eigenvalues).
%   5. THE FIGURES draw from the stage's output for both networks (PNG and
%      _table.md written; tags with _EEx1p5 for that network), refuse a
%      network the file lacks, and error on an empty run directory;
%      fig_eigdens_theory_check draws without a run directory.
%
% Asserts; prints one PASS line. About a minute. Assumes setup_paths has run
% and a pool allocation is set (pool_set(0) will do).
%
% See also: eigenvalue_density_stage_plan, run_eigenvalue_density_stage,
%           eigdens_load_runs, make_synthetic_outlier_spectra

t_all = tic;

%% 1. The theory
C = eigdens_check_theory('n', 200, 'n_draws', 10, 'seed', 3);
for c = 1:numel(C)
    assert(C(c).max_diff < 3.5 * C(c).sampling_error, 'theory vs sampled, case %s: %.4f against error %.4f', ...
        C(c).name, C(c).max_diff, C(c).sampling_error);
    if c <= 4; assert(C(c).max_diff_solvers < 1e-3, 'solvers disagree, case %s', C(c).name); end
end
assert(C(4).max_diff_harris < 1e-3, 'Harris closed form vs column solver: %.2g', C(4).max_diff_harris);
T = eigdens_radial_theory(ones(100, 1), linspace(0, 1, 51));
assert(max(abs(T.n_less - linspace(0, 1, 51).^2)) < 1e-9 && abs(T.edge_rel - 1) < 1e-12, 'uniform gains: the circular law');
T = eigdens_radial_theory([zeros(50, 1); ones(50, 1)], [0 0.5 sqrt(0.5) * 0.5]);
assert(abs(T.p0 - 0.5) < 1e-12 && abs(T.n_less(1) - 0.5) < 1e-12 && abs(T.n_less(3) - (0.5 + 0.5 * 0.25)) < 1e-6, ...
    'half the gains zero: an atom of 1/2 and a uniform disk of radius sqrt(1/2)');
mu = eigdens_measures('uniform', 1, 100);
assert(abs(mu.shell - 0.36) < 1e-12 && abs(mu.r50 - sqrt(0.5)) < 1e-12 && mu.n_unstable == 0, 'the uniform disk''s measures');

%% 2. The plan
cfg = struct('run_mode', 'essential', 'stages', {{'outlier_spectra', 'eigenvalue_density'}}, 'verbose', 'near-none');
pl = eigenvalue_density_stage_plan(cfg);
assert(strcmp(pl.folder, 'eigenvalue_density') && strcmp(pl.input_folder, 'outlier_spectra') && ...
    strcmp(pl.data_file, 'eigenvalue_density_data.mat') && isequal(pl.grids.u, 0:0.01:1.4));
for m = run_mode_names()
    c = cfg; c.run_mode = m{1};
    eigenvalue_density_stage_plan(c);                % every mode is accepted
end
refused = {setfield(cfg, 'stages', {'eigenvalue_density'}), setfield(cfg, 'eigenvalue_density', struct('n', 3)), ...
    setfield(cfg, 'run_mode', 'huge'), rmfield(cfg, 'run_mode')}; %#ok<SFLD>
for k = 1:numel(refused)
    assert(throws(@() eigenvalue_density_stage_plan(refused{k}), 'eigenvalue_density_stage_plan:BadConfig'), 'refusal %d', k);
end

%% 3. The stage on contract-conformant files
run_dir = tempname; empty_run = tempname; fig_dir = tempname;
cleanup = onCleanup(@() cellfun(@rmdir_if, {run_dir, empty_run, fig_dir}));
mkdir(run_dir); mkdir(empty_run);
make_synthetic_outlier_spectra(run_dir);
index = eigdens_load_runs(run_dir);
assert(height(index) == 12);
t0 = tic;
detail = run_eigenvalue_density_stage(cfg, run_dir, pl);
t_stage = toc(t0);
mat = fullfile(run_dir, 'eigenvalue_density', 'eigenvalue_density_data.mat');
assert(isfile(mat) && contains(detail, mat));
D = load(mat);
assert(numel(D.runs) == 12 && isequal({D.runs.network}, index.network') && isequal({D.runs.condition}, index.condition'));
assert(isequal({D.networks.name}, {'reference', 'EEx1p5'}) && abs(D.networks(1).harris.g) < 1e-12, ...
    'two networks; W alone has a uniform density (Harris g = 0) in the reference network');
for k = 1:numel(D.runs)
    r = D.runs(k);
    K = numel(r.R_eff);
    assert(K == 3 && isequal(size(r.cdf_M), [3, numel(pl.grids.u)]) && isequal(size(r.dens_M), [3, numel(pl.grids.u2_edges) - 1]) ...
        && isequal(size(r.m_M), [4 3]) && isequal(size(r.gains), [80 3]) && numel(r.ev_example) == 80, 'run %d: shapes', k);
    assert(max(abs(r.R_eff - r.R_eff_saved)) < 1e-9 * max(r.R_eff), 'run %d: R_eff from the saved gains', k);
    assert(all(isfinite(r.full.n_band_x)), 'run %d: the dense pass was read', k);
    X = eigdens_load_runs(run_dir, index, k);
    for s = 1:K                                   % the counts are the spectrum's own
        ev = X.ev_M(:, s);
        assert(r.m_M(4, s) == nnz(real(ev) > 1) && r.m_M(3, s) == nnz(real(ev) > 0.8 & real(ev) <= 1));
        assert(abs(r.p0(s) - mean(X.gain(:, s) == 0)) < 1e-12);
    end
end
% No input, a wrong network: errors, never a fallback.
assert(throws(@() eigdens_load_runs(empty_run), 'eigdens_load_runs:NoInput'));
assert(throws(@() eigdens_load_runs(fullfile(empty_run, 'nope')), 'eigdens_load_runs:NoInput'));
assert(throws(@() run_eigenvalue_density_stage(cfg, empty_run, pl), 'eigdens_load_runs:NoInput'));
f1 = fullfile(run_dir, 'outlier_spectra', index.file{1});
L = load(f1, 'R'); R = L.R; R.rng_seeds = [5 6];           % another draw of W than the saved one
bad_dir = fullfile(empty_run, 'bad'); mkdir(fullfile(bad_dir, 'outlier_spectra'));
save(fullfile(bad_dir, 'outlier_spectra', index.file{1}), 'R');
runs = index(1, :); save(fullfile(bad_dir, 'outlier_spectra', 'ospec_index.mat'), 'runs');
assert(throws(@() eigdens_load_runs(bad_dir, eigdens_load_runs(bad_dir), 1), 'eigdens_load_runs:WrongNetwork'));

%% 4. Theory with the saved gains reproduces the measured distribution
inside = pl.grids.u <= 1;
worst_state = max([D.runs.theory_gap]);
assert(worst_state < 0.15, 'largest single-state gap %.3f', worst_state);
worst_mean = 0;
for w = {'reference', 'EEx1p5'}
    for cnd = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}
        rs = D.runs(strcmp({D.runs.network}, w{1}) & strcmp({D.runs.condition}, cnd{1}));
        gap = max(abs(mean(cat(1, rs.cdf_M), 1, 'omitnan') - mean(cat(1, rs.cdf_theory), 1, 'omitnan')) .* inside);
        worst_mean = max(worst_mean, gap);
    end
end
assert(worst_mean < 0.08, 'mean measured vs mean theory: %.3f', worst_mean);

%% 5. The figures
figs = {@fig_eigdens_radial, 'Fig_Eigdens_Radial'; @fig_eigdens_gains, 'Fig_Eigdens_Gains'; ...
    @fig_eigdens_measures, 'Fig_Eigdens_Measures'; @fig_eigdens_plane, 'Fig_Eigdens_Plane'};
t0 = tic;
for i = 1:size(figs, 1)
    for w = {'reference', 'EEx1p5'}
        o = fullfile(fig_dir, sprintf('%s_%s', figs{i, 2}, w{1}));
        out = figs{i, 1}('run_dir', run_dir, 'out_dir', o, 'network', w{1}, 'visible', false, 'verbose', 'near-none');
        close(out.figs);
        tag = figs{i, 2};
        if ~strcmp(w{1}, 'reference'); tag = [tag '_' w{1}]; end %#ok<AGROW>
        assert(isfile(fullfile(o, [tag '.png'])) && isfile(fullfile(o, [tag '_table.md'])), '%s: PNG and table', tag);
        assert(strcmp(out.source, mat));
    end
    assert(throws_prefix(@() figs{i, 1}('run_dir', empty_run, 'save', false, 'visible', false), 'resolve_data_file:'), ...
        '%s on an empty run directory', figs{i, 2});
    assert(throws(@() figs{i, 1}('run_dir', run_dir, 'network', 'EEx9p9', 'save', false, 'visible', false), ...
        [func2str(figs{i, 1}) ':NoNetwork']), '%s: an absent network', figs{i, 2});
end
out = fig_eigdens_radial('run_dir', run_dir, 'conditions', {'sfa3_stdc2', 'no_adaptation'}, 'save', false, 'visible', false);
close(out.figs);
o = fullfile(fig_dir, 'theory');
out = fig_eigdens_theory_check('out_dir', o, 'n', 100, 'n_draws', 4, 'visible', false);
close(out.figs);
assert(isfile(fullfile(o, 'Fig_Eigdens_Theory_Check.png')) && isfile(fullfile(o, 'Fig_Eigdens_Theory_Check_table.md')));
t_fig = toc(t0);

fprintf(['PASS: test_eigenvalue_density_stage -- theory vs sampled matrices, the plan, the stage on 12 runs (%.1f s), ' ...
    'theory with the saved gains (worst state %.3f, worst mean %.3f), the figures (%.0f s); %.0f s in all.\n'], ...
    t_stage, worst_state, worst_mean, t_fig, toc(t_all));

%% ------------------------------------------------------------------------
function tf = throws(fn, id)
try
    fn();
    tf = false;
catch err
    tf = strcmp(err.identifier, id);
    if ~tf; fprintf(2, 'expected %s, got %s: %s\n', id, err.identifier, err.message); end
end
end

function tf = throws_prefix(fn, prefix)
try
    fn();
    tf = false;
catch err
    tf = startsWith(err.identifier, prefix);
end
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end
