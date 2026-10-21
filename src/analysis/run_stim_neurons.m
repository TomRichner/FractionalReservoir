function out_dir = run_stim_neurons(cfg, opts)
% RUN_STIM_NEURONS Rerun one sine-train simulation and keep single-neuron traces.
%
%   out_dir = RUN_STIM_NEURONS(cfg)          % cfg from sine_paired_stim_config()
%   out_dir = RUN_STIM_NEURONS(cfg, 'condition', 'sfa3_std2', 'period', 2)
%
% run_stim_evoked saves population means only. For the neuron-level review
% figure (fig_sine_stim_neurons, TR 2026-10-20) this repeats ONE of its runs
% -- same preset, condition, seeds, stimulus, fs and T_range, so the same
% trajectory -- and keeps, over a display window, every quantity for a
% fixed random sample of neurons plus the means of the three classes
% (driven E, undriven E, I):
%
%   u      external input u_i (the pulses; zero on undriven neurons)
%   noise  the Wiener term as an EQUIVALENT INPUT, sigma_u * xi_1 / sqrt(dt):
%          per step the integrator adds (sigma_u/tau_d) * sqrt(dt) * xi_1 to
%          x (sde_fixed_step, sra1; its second normal xi_2 only shapes the
%          drift's intermediate stage), which is what an input of that value
%          held for the step would add. Regenerated from the model's
%          noise_seed exactly as build_noise draws it.
%   x, r, syn (r prod_m b_m, route 1), sfa ((c/K) sum_k a_k), std (prod_m b_m)
%
% at the full rate fs. The rerun is checked against the saved population
% mean of run_stim_evoked (the undriven E neurons' x) before anything is
% written. Output: <run_dir>/stim_neurons/stim_neurons_data.mat; that folder
% must not exist yet.
%
% See also: run_stim_evoked, fig_sine_stim_neurons, sine_paired_stim_config

arguments
    cfg struct
    opts.condition (1,:) char   = 'sfa3_std2'
    opts.period    (1,1) double = 2
    opts.window    (1,2) double = [28 38]      % s of stimulation time kept
    opts.n_sample  (1,1) double = 25
    opts.sample_seed (1,1) double = 20261020
end
setup_paths();
project_root = fileparts(which('setup_paths'));
run_dir = cfg.run_dir;
if isempty(regexp(run_dir, '^([A-Za-z]:[\\/]|[\\/])', 'once')), run_dir = fullfile(project_root, run_dir); end
out_dir = fullfile(run_dir, 'stim_neurons');
if isfolder(out_dir)
    error('run_stim_neurons:Exists', '%s exists; delete it to rerun.', out_dir);
end
t_run = tic;

% The stimulus exactly as run_stim_evoked builds it for a sine protocol.
s = cfg.stim;
k_per = find(cfg.sine_periods == opts.period, 1);
input = struct('intrinsic_drive', 0, 'amp', s.amp, 'neuron_frac', s.neuron_frac, ...
    'target_types', {s.target_types}, 'generator', @generate_am_pulse_train, ...
    'width_s', s.width_s, 'pulse_period', s.pulse_period, 'mod_period', opts.period, ...
    'mod_depth', s.mod_depth, 't_first', 0);
stim_T = opts.period * cfg.sine_cycles(k_per);
extra = {};
if ~isempty(cfg.sigma_u_noise), extra = {'sigma_u_noise', cfg.sigma_u_noise}; end
model = build_from_preset(cfg.preset_name, opts.condition, ...
    'T_range', [-cfg.settle_s, stim_T], 'fs', cfg.fs, 'rng_seeds', cfg.rng_seeds, ...
    'input_config', input, 'lya_method', 'none', 'store_full_state', true, ...
    'store_decimated_state', false, 'verbose', cfg.verbose, extra{:});
params = model.get_params();
model.run();
[~, ~, info] = generate_am_pulse_train(params, stim_T, cfg.fs, cfg.rng_seeds(2), input);

% Classes and the sample
n = model.n;
is_E = false(n, 1); is_E(params.type_indices{1}) = true;
is_drv = false(n, 1); is_drv(info.stim_idx) = true;
cls = 3 * ones(n, 1); cls(is_E & is_drv) = 1; cls(is_E & ~is_drv) = 2;   % 1 driven E, 2 undriven E, 3 I
class_names = {'driven E', 'undriven E', 'I'};
stream = RandStream('twister', 'Seed', opts.sample_seed);
sample = sort(randperm(stream, n, opts.n_sample));

% The window, full rate
t_all = model.t_out(:)';
win = t_all >= opts.window(1) - 1e-9 & t_all <= opts.window(2) + 1e-9;
iw = find(win);
[x, a, b, r, so] = SRNNCellTypePairs.unpack_and_compute_states(model.S_out(iw, :), params);
full = struct('u', model.u_ex(:, iw), 'x', zeros(n, numel(iw)), 'r', zeros(n, numel(iw)), ...
    'syn', zeros(n, numel(iw)), 'sfa', zeros(n, numel(iw)), 'std', ones(n, numel(iw)));
for q = 1:numel(params.cell_type_names)
    nm = params.cell_type_names{q}; ix = params.type_indices{q}; nq = numel(ix);
    posts = fieldnames(so.(nm));
    full.x(ix, :) = x.(nm); full.r(ix, :) = r.(nm);
    full.syn(ix, :) = so.(nm).(posts{1});
    if ~isempty(a.(nm)), full.sfa(ix, :) = params.c_eff(q) * reshape(sum(a.(nm), 2), nq, []); end
    B = b.(nm).(posts{1});
    if ~isempty(B), full.std(ix, :) = reshape(prod(B, 2), nq, []); end
end

% The Wiener term, regenerated as build_noise draws it (xi1 first, then xi2).
seed = model.noise_seed; if isempty(seed), seed = model.rng_seeds(1) + 224737; end
state = rng; rng(seed);
xi1 = randn(n, numel(model.t_ex) - 1);
rng(state);
step = iw; step(step > size(xi1, 2)) = size(xi1, 2);     % step j runs from t_ex(j) to t_ex(j+1)
full.noise = model.sigma_u_noise * xi1(:, step) * sqrt(cfg.fs);
clear xi1

% Check against the saved population mean before writing anything.
S = load(fullfile(run_dir, 'stim_evoked', 'stim_evoked_data.mat'), 'results');
kk = find(strcmp({S.results.condition}, opts.condition) & ...
    strcmp({S.results.protocol}, sprintf('sine_%gs', opts.period)), 1);
Rk = S.results(kk);
t_mid = t_all(iw(1:cfg.readout_deci:end));
saved = interp1(Rk.t, double(Rk.readouts.x.E_nonstim), t_mid);
mine = mean(full.x(cls == 2, 1:cfg.readout_deci:end), 1);
mismatch = max(abs(saved - mine));
assert(mismatch < 1e-4, 'run_stim_neurons:Mismatch', ...
    'The rerun differs from the saved run (max |dx| %.3g).', mismatch);

vars = {'u', 'noise', 'x', 'r', 'syn', 'sfa', 'std'};
traces = struct(); means = struct();
for v = 1:numel(vars)
    traces.(vars{v}) = single(full.(vars{v})(sample, :));
    for c = 1:3
        means.(vars{v})(c, :) = single(mean(full.(vars{v})(cls == c, :), 1));
    end
end
counts = arrayfun(@(c) nnz(cls == c), 1:3);
sample_counts = arrayfun(@(c) nnz(cls(sample) == c), 1:3);
data = struct('t', t_all(iw), 'traces', traces, 'means', means, 'sample', sample, ...
    'sample_class', cls(sample)', 'class_names', {class_names}, 'counts', counts, ...
    'sample_counts', sample_counts, 'pulse_times', info.pulse_times, ...
    'pulse_amps', info.pulse_amps, 'mod_period', opts.period, 'condition', opts.condition, ...
    'mismatch_vs_saved', mismatch, 'noise_definition', ...
    'sigma_u * xi_1 / sqrt(dt): the per-step Wiener increment of x expressed as an equivalent input (input units)');
mkdir(out_dir);
settings = struct('cfg', cfg, 'opts', opts, 'stim_T', stim_T, 'runtime_s', toc(t_run), ...
    'created', char(datetime('now')), 'git', capture_git_provenance(out_dir, project_root));
save(fullfile(out_dir, 'stim_neurons_data.mat'), 'data', 'settings', '-v7.3');
vprintf(cfg.verbose, 'minimal', '  run_stim_neurons: %s %g s, %.0f s, rerun matches saved run to %.2g\n', ...
    opts.condition, opts.period, settings.runtime_s, mismatch);
end
