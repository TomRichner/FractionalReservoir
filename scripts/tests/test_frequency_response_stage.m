%% test_frequency_response_stage.m
% The network frequency response (FD-041): the opt-in stage
% 'frequency_response' (frequency_response_stage_plan,
% run_frequency_response_stage), its computation (src/analysis/frequency_response/)
% and its four figures, on 60-neuron networks.
%
% The stage reads the calibration states of the stage outlier_spectra
% (FD-040). That stage was built in parallel, so the input here is made by
% the fixture fresp_test_make_ospec: files in the contract's layout from
% tiny networks (two networks x three conditions x two seeds, three states
% each).
%
% Checks:
%   1. The plan in EVERY run mode: the folder, file and input, the two
%      networks and the default grids, the same in every mode. Overrides
%      are taken. An unknown field, a bad grid, a bad rate, a missing run
%      mode, and a bundle whose cfg.stages lacks outlier_spectra are
%      refused (and that last one accepted with input_exists).
%   2. The stage on the fixture's files: one entry per (network, condition,
%      seed) of the index, the curves of the right size and finite, the
%      matched operating point per adapting condition and seed. A condition
%      without slow variables does not depend on frequency and has no
%      checks and no single-neuron responses.
%   3. CLOSED FORM AGAINST THE JACOBIAN, to round-off: the gains of
%      fresp_gain against fresp_schur (the model's own Jacobian blocks) at
%      real and complex frequencies, with both mechanisms free and each
%      alone, on both adapting conditions; M(z) is the weights times one
%      gain per neuron and route (off-diagonal remainder at round-off);
%      the single-neuron responses against fresp_neuron_schur; and the same
%      numbers as the stage saved them.
%   4. EXACTNESS: eigenvalues of the dense full Jacobian that are not also
%      eigenvalues of the slow block satisfy the reduced problem, 1 + tau_d z
%      an eigenvalue of M(z) (smallest singular value at round-off); and
%      the frozen limit: M at a very high frequency is the frozen-gain
%      matrix I + tau_d J_xx.
%   5. The loader: the stage errors on a run directory without the
%      outlier_spectra index and on a network the index lacks; a run file
%      that lacks a contract field is refused by name.
%   6. The four figures draw from the stage's output (PNG and _table.md
%      each), take 'conditions', and error on an empty run directory.
%
% Asserts; prints one PASS line. About a minute, serial. Assumes
% setup_paths has run and a pool allocation is set (pool_set(0) will do).
%
% See also: frequency_response_stage_plan, run_frequency_response_stage,
%           fresp_run_one, fresp_load_states, fresp_test_make_ospec

%% 1. The plan in every run mode
modes = run_mode_names();
for m = 1:numel(modes)
    pl = frequency_response_stage_plan(struct('run_mode', modes{m}));
    assert(strcmp(pl.folder, 'frequency_response') && strcmp(pl.file, 'frequency_response_data.mat') && ...
        strcmp(pl.input_folder, 'outlier_spectra') && isequal(pl.networks, {'reference', 'EEx1p5'}), '%s: the plan', modes{m});
    o = pl.opts;
    assert(numel(o.f) == 61 && abs(o.f(1) - 1e-3) < 1e-15 && abs(o.f(end) - 100) < 1e-9 && numel(o.f_single) == 201 && ...
        isempty(o.exact_states) && isempty(o.shown_state) && o.n_lead == 20 && o.n_sample == 200 && ...
        isequal(o.exact_f, [0.001 0.03 0.5 3 30]) && numel(o.sigma) == 201 && numel(o.w) == 121 && ...
        o.plane_all_states && o.r0 == 0.25, '%s: the settings', modes{m});
    assert(~isfield(o, 'networks') && ~isfield(o, 'input_exists'));
end
pl = frequency_response_stage_plan(struct('run_mode', 'fast', 'frequency_response', ...
    struct('networks', {{'reference'}}, 'f', [0.1 1 10], 'exact_states', [1 3], 'shown_state', 3, 'r0', 0.5)));
assert(isequal(pl.networks, {'reference'}) && isequal(pl.opts.f, [0.1 1 10]) && ...
    isequal(pl.opts.exact_states, [1 3]) && pl.opts.shown_state == 3 && pl.opts.r0 == 0.5);
refused = { ...
    struct('run_mode', 'fast', 'frequency_response', struct('freqs', 1)), ...            an unknown field
    struct('run_mode', 'fast', 'frequency_response', struct('f', [1 0.5])), ...          not ascending
    struct('run_mode', 'fast', 'frequency_response', struct('f', [0 1])), ...            zero frequency
    struct('run_mode', 'fast', 'frequency_response', struct('w', [-1 0 1])), ...         lower half plane
    struct('run_mode', 'fast', 'frequency_response', struct('r0', 1.5)), ...             not a rate
    struct('run_mode', 'fast', 'frequency_response', struct('networks', {{}})), ...      no network
    struct('run_mode', 'fast', 'frequency_response', struct('exact_states', 0)), ...     not a state
    struct('run_mode', 'nonesuch'), ...                                                  no such mode
    struct('run_mode', 'fast', 'stages', {{'eig_heatmap', 'frequency_response'}})};      % the input stage is not run
for k = 1:numel(refused)
    assert(throws(@() frequency_response_stage_plan(refused{k}), 'frequency_response_stage_plan:BadConfig'), ...
        'setting %d was not refused', k);
end
frequency_response_stage_plan(struct('run_mode', 'fast', 'stages', {{'outlier_spectra', 'frequency_response'}}));
frequency_response_stage_plan(struct('run_mode', 'fast', 'stages', {{'frequency_response'}}, ...
    'frequency_response', struct('input_exists', true)));

%% 2. The stage on contract-conformant files
tmp = tempname; tmpf = tempname; empty_run = tempname; broken = tempname;
cleanup = onCleanup(@() cellfun(@rmdir_if, {tmp, tmpf, empty_run, broken}));
networks = {'reference', 'EEx1p5'};
conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
seeds = [1 2; 2 3];
t0 = tic;
fresp_test_make_ospec(tmp, 'networks', networks, 'conditions', conditions, 'seeds', seeds, 'n', 60, 'n_cal', 3);
t_fix = toc(t0);

cfg = struct('run_mode', 'essential', 'verbose', 'near-none');
cfg.frequency_response = struct('f', logspace(-3, 2, 21), 'f_single', logspace(-3, 2, 41), ...
    'sigma', -12:0.5:8, 'w', 0:1:30, 'n_sample', 60);
plan = frequency_response_stage_plan(cfg);
t0 = tic;
detail = run_frequency_response_stage(cfg, tmp, plan);
t_run = toc(t0);
mat = fullfile(tmp, 'frequency_response', 'frequency_response_data.mat');
assert(isfile(mat) && contains(detail, mat) && contains(detail, '12 runs'), 'the stage writes its own folder');
D = load(mat);
S = D.settings;
assert(isequal(S.networks, networks) && isequal(S.conditions, conditions) && isequal(S.seeds, seeds) && ...
    isequal(S.f, plan.opts.f) && S.r0 == 0.25);
assert(numel(D.runs) == 12 && numel(D.matched) == 4, 'one entry per run; matched per adapting condition and seed');
nf = numel(S.f);
for k = 1:numel(D.runs)
    A = D.runs(k);
    tag = sprintf('%s / %s / seed [%d %d]', A.network, A.condition, A.seed);
    assert(ismember(A.network, networks) && ismember(A.condition, conditions) && ismember(A.seed, seeds, 'rows'), tag);
    assert(A.n == 60 && numel(A.t_cal) == 3 && isequal(A.exact_states, [1 2]) && A.shown_state == 2, tag);
    assert(isequal(size(A.R.both), [3 nf]) && isequal(size(A.max_re), [3 nf]) && isequal(size(A.K_eig), [2 nf 3]), tag);
    assert(all(isfinite([A.R.both(:); A.R.sfa(:); A.R.std(:); A.max_re(:); A.R_inf(:)])) && all(A.R.both(:) >= 0), tag);
    assert(numel(A.plane) == 2 && numel(A.plane(1).lam) == A.N && isequal(size(A.plane(2).R), [numel(S.sigma), numel(S.w)]), tag);
    assert(abs(A.plane(1).abscissa - max(real(double(A.plane(1).lam)))) < 1e-5 * max(1, abs(A.plane(1).abscissa)), tag);
    assert(A.has_slow == ~strcmp(A.condition, 'no_adaptation'), tag);
    if ~A.has_slow
        % Nothing to eliminate: no dependence on frequency, no checks.
        assert(A.N == A.n && max(abs(A.R.both - A.R_inf), [], 'all') < 1e-12 && ...
            max(abs(A.max_re - A.max_re_inf), [], 'all') < 1e-12 && isempty(A.single) && isempty(A.exact(1).sigma), tag);
    else
        % The frozen radius is the high-frequency limit; adaptation and
        % depression lower the radius at zero frequency.
        assert(A.N > A.n && all(A.R_zero.both < A.R_inf) && all(abs(A.R.both(:, end) - A.R_inf) < 0.02 * A.R_inf), tag);
        assert(numel(A.single) == 3 && isequal({A.single.label}, {'E median', 'I median', 'E driven'}), tag);
    end
end
ee = D.runs(find(strcmp({D.runs.network}, 'EEx1p5'), 1)); rf = D.runs(find(strcmp({D.runs.network}, 'reference'), 1));
assert(ee.structural.has_outlier && ~rf.structural.has_outlier, 'E-to-E x 1.5 has a structural outlier; the reference network has none');
for k = 1:numel(D.matched)
    mt = D.matched(k);
    assert(strcmp(mt.network, 'reference') && ismember(mt.condition, conditions(2:3)) && ...
        size(mt.G.both, 2) == numel(S.f_single) && all(isfinite(mt.G.both(:))) && all(mt.F > 0));
    % At high frequency the slow variables do nothing: the gain is the frozen gain.
    assert(max(abs(double(mt.G.both(:, end)) - 1)) < 0.02 && all(mt.G0.both < 1));
end

%% 3. Closed form against the Jacobian, and 4. exactness
T = fresp_load_states(tmp);
assert(height(T) == 12 && isequal(unique(T.network, 'stable')', networks));
worst = struct('closed', 0, 'offdiag', 0, 'single', 0, 'sigma', 0, 'frozen', 0);
zs = [2i * pi * [0.001 0.05 1 20], 0.4 + 3i, -0.3 + 0.7i, 0];
frees = {struct('sfa', true, 'std', true), struct('sfa', true, 'std', false), struct('sfa', false, 'std', true)};
for c = 2:3
    row = T(strcmp(T.network, 'EEx1p5') & strcmp(T.condition, conditions{c}) & T.seed1 == 1, :);
    st = fresp_load_states(tmp, row);
    assert(isequal(size(st.S), [st.params.state_layout.N, 3]) && numel(st.cont) == 3 && all(imag(st.cont) >= 0));
    ws = fresp_weight_stats(st.params);
    Sk = st.S(:, 2);
    J = SRNNCellTypePairs.compute_Jacobian_fast(Sk, st.params);
    ok = ws.col_sq > 0;
    for q = 1:numel(zs)
        for fr = 1:3
            [~, dn, off] = fresp_schur(J, st.params, ws, zs(q), frees{fr});
            dc = fresp_gain(Sk, st.params, zs(q), frees{fr});
            worst.closed = max(worst.closed, max(abs(dc(ok) - dn(ok))) / max(abs(dn(ok))));
            worst.offdiag = max(worst.offdiag, off);
        end
    end
    % The frozen limit: at a very high frequency M(z) is I + tau_d J_xx.
    x = st.params.state_layout.x;
    Mf = eye(st.params.n) + st.params.tau_d * full(J(x, x));
    [d_hi, d_inf] = fresp_gain(Sk, st.params, 2i * pi * 1e9);
    worst.frozen = max([worst.frozen, norm(fresp_matrix(d_hi, ws) - Mf, 'fro') / norm(Mf, 'fro'), ...
        norm(fresp_matrix(d_inf, ws) - Mf, 'fro') / norm(Mf, 'fro')]);
end
for k = find([D.runs.has_slow])
    A = D.runs(k);
    E = A.exact;
    worst.closed = max([worst.closed, E.closed_diff]);
    worst.offdiag = max([worst.offdiag, E.offdiag]);
    worst.single = max([worst.single, A.single.max_diff]);
    for e = 1:numel(E)
        off_slow = E(e).dslow > 1e-6;
        assert(nnz(off_slow) >= 20, 'too few eigenvalues off the slow block to judge exactness');
        worst.sigma = max(worst.sigma, max(E(e).sigma(off_slow)));
        assert(any(E(e).kind == 1) && any(E(e).kind == 2) && any(E(e).kind == 3));
    end
    assert(A.pop.route_diff < 1e-12, 'in this preset a neuron has one gain for all its routes');
end
assert(worst.closed < 1e-10, 'closed form against the Jacobian''s blocks: %.2g', worst.closed);
assert(worst.offdiag < 1e-10, 'M(z) is not the weights times one gain per neuron and route: %.2g', worst.offdiag);
assert(worst.single < 1e-9, 'single neuron, closed form against its Jacobian blocks: %.2g', worst.single);
assert(worst.sigma < 1e-8, 'eigenvalues of the full Jacobian do not satisfy the reduced problem: %.2g', worst.sigma);
assert(worst.frozen < 1e-8, 'the high-frequency limit is not the frozen-gain matrix: %.2g', worst.frozen);

%% 5. The loader and the stage refuse what is not there
mkdir(empty_run);
assert(throws(@() run_frequency_response_stage(cfg, empty_run, plan), 'fresp_load_states:NoIndex'), ...
    'a run directory without the outlier_spectra index is an error');
assert(~isfolder(fullfile(empty_run, 'frequency_response')), 'a refused stage writes nothing');
plan_x = frequency_response_stage_plan(struct('run_mode', 'essential', ...
    'frequency_response', struct('networks', {{'reference', 'EEx2p0'}})));
assert(throws(@() run_frequency_response_stage(cfg, tmp, plan_x), 'run_frequency_response_stage:NoSuchNetwork'));
copyfile(fullfile(tmp, 'outlier_spectra'), fullfile(broken, 'outlier_spectra'));
L = load(fullfile(broken, 'outlier_spectra', T.file{1}));
R = rmfield(L.R, 'cal');
save(fullfile(broken, 'outlier_spectra', T.file{1}), 'R');
assert(throws(@() fresp_load_states(broken, T(1, :)), 'fresp_load_states:Contract'), ...
    'a run file without the calibration states is refused');

%% 6. The figures
t0 = tic;
figs = {@fig_freqresp_network, 'Fig_Freqresp_Network'; @fig_freqresp_plane, 'Fig_Freqresp_Plane'; ...
    @fig_freqresp_phase, 'Fig_Freqresp_Phase'; @fig_freqresp_single_neuron, 'Fig_Freqresp_Single_Neuron'};
for k = 1:size(figs, 1)
    name = func2str(figs{k, 1});
    out_dir = fullfile(tmpf, name);
    out = figs{k, 1}('run_dir', tmp, 'out_dir', out_dir, 'visible', false, 'verbose', 'near-none');
    assert(isscalar(out.figs) && isgraphics(out.figs) && strcmp(out.source, mat), name);
    close(out.figs);
    assert(isfile(fullfile(out_dir, [figs{k, 2} '.png'])) && isfile(fullfile(out_dir, [figs{k, 2} '_table.md'])), ...
        '%s: the PNG and the table', name);
    tbl = fileread(fullfile(out_dir, [figs{k, 2} '_table.md']));
    assert(contains(tbl, 'sfa3_stdc2') && contains(tbl, '| condition |') && contains(tbl, 'frequency_response_data.mat'), ...
        '%s: the table names the conditions and its data', name);
    % The bundle's condition list is taken, in its order.
    out = figs{k, 1}('run_dir', tmp, 'save', false, 'visible', false, 'verbose', 'near-none', ...
        'conditions', {'sfa3_stdc2', 'no_adaptation'});
    close(out.figs);
    assert(throws_prefix(@() figs{k, 1}('run_dir', tmp, 'save', false, 'visible', false, 'conditions', {'nonesuch'}), ...
        'select_conditions:'), '%s: an unknown condition is an error', name);
    assert(throws_prefix(@() figs{k, 1}('run_dir', empty_run, 'save', false, 'visible', false), ...
        'resolve_data_file:'), '%s: a run without the stage is an error', name);
end
tbl = fileread(fullfile(tmpf, 'fig_freqresp_plane', 'Fig_Freqresp_Plane_table.md'));
assert(contains(tbl, 'seed [1 2]') && contains(tbl, 'calibration state 2'), 'the plane table states the seed and state drawn');
out = fig_freqresp_plane('run_dir', tmp, 'save', false, 'visible', false, 'seed', [2 3]); close(out.figs);
assert(throws(@() fig_freqresp_plane('run_dir', tmp, 'save', false, 'visible', false, 'seed', [9 9]), ...
    'fig_freqresp_plane:NoSuchSeed'));
t_fig = toc(t0);

fprintf(['PASS: test_frequency_response_stage -- the plan in every mode, refusals, the stage on contract files, ' ...
    'closed form against the Jacobian %.1e (off-diagonal %.1e, single neuron %.1e), exactness %.1e, frozen limit %.1e, ' ...
    'the loader, four figures (fixture %.0f s, stage %.0f s, figures %.0f s).\n'], ...
    worst.closed, worst.offdiag, worst.single, worst.sigma, worst.frozen, t_fix, t_run, t_fig);

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
    if ~tf; fprintf(2, 'expected %s..., got %s: %s\n', prefix, err.identifier, err.message); end
end
end

function rmdir_if(d)
if isfolder(d); rmdir(d, 's'); end
end
