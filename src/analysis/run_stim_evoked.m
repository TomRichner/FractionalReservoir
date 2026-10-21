function out_dir = run_stim_evoked(cfg)
% RUN_STIM_EVOKED Simulate the human sine and paired-pulse stimulation protocols.
%
%   out_dir = RUN_STIM_EVOKED(cfg)      % cfg from sine_paired_stim_config()
%
% The model counterpart of the SEEG experiment in brain-stim-analysis: one
% network, driven at a fixed set of neurons (the "stimulated site") by
%
%   - 5 Hz pulse trains whose amplitude is sinusoidally modulated with a
%     2, 5 or 10 s period (generate_am_pulse_train; the first pulse at the
%     envelope's trough, 70 % of the peak), and
%   - the paired-pulse protocol, a single pulse then pairs at 2, 10 and
%     100 ms, one event every 2.5 s (generate_paired_pulse_input),
%
% under every adaptation condition of the preset, on the same W, setpoints,
% SFA ladders, noise and driven subset (cfg.rng_seeds). Each run starts
% cfg.settle_s before the first pulse so the network settles on noise alone.
%
% WHAT IS SAVED. Not the state: population means only, which is all the
% figure needs and keeps a 70-s run at 2 kHz to a few MB. For every
% condition x protocol, at the readout rate fs / cfg.readout_deci, the mean
% over neurons of
%
%   x     dendritic potential                  (the "CCEP-like" signal)
%   r     rate, phi(x - (c/K) sum_k a_k)       (pre-depression)
%   syn   synaptic output r * prod b, route 1
%   sfa   SFA feedback (c/K) sum_k a_k
%   std   depression prod_m b_m, route 1
%
% in four groups, <type>_stim and <type>_nonstim for each cell type (a
% group with no neurons is NaN). Route 1 is every route here: the run
% asserts SRNNCellTypePairs.routes_identical. The equations are in
% docs/EquationsParametersDocs/Equations_stability_paper.md.
%
% The state is unpacked in chunks from S_out (store_full_state on,
% store_decimated_state off) rather than through plot_data, which would hold
% a second full copy of the trajectory: the 70-s runs peak near 6 GB this
% way instead of about 12.
%
% Output: <cfg.run_dir>/stim_evoked/stim_evoked_data.mat with
%   results  - struct array, one per condition x protocol (fields below)
%   settings - cfg, the resolved stimulus configs, git provenance, timing
% The figure (fig_sine_paired_stim) and stim_evoked_measures read it.
%
% Lyapunov estimation is off (lya_method 'none'); no parallel pool is used.
%
% See also: sine_paired_stim_config, stim_evoked_measures,
%           fig_sine_paired_stim, generate_am_pulse_train,
%           generate_paired_pulse_input, run_paper_illustrations

arguments
    cfg struct
end
setup_paths();
project_root = fileparts(which('setup_paths'));
run_dir = cfg.run_dir;
if ~is_absolute(run_dir), run_dir = fullfile(project_root, run_dir); end
if isfolder(run_dir) && numel(dir(run_dir)) > 2
    error('run_stim_evoked:RunDirNotEmpty', ...
        'run_dir %s is not empty; delete it or name another.', run_dir);
end
out_dir = fullfile(run_dir, 'stim_evoked');
mkdir(out_dir);
log_file = fullfile(run_dir, 'command_window.log');
log_cleanup = vlog('open', log_file); %#ok<NASGU>

[~, ~, conds] = srnn_param_preset(cfg.preset_name);
cond_names = cellfun(@(c) c.name, conds, 'UniformOutput', false);
if ~isempty(cfg.conditions), cond_names = cfg.conditions; end
titles = srnn_condition_titles();

% Protocols: one per sine period, then paired pulse.
protocols = struct('name', {}, 'kind', {}, 'mod_period', {}, 'n_cycles', {}, ...
    'cycle_len', {}, 'stim_T', {});
for k = 1:numel(cfg.sine_periods)
    Tm = cfg.sine_periods(k);
    protocols(end + 1) = struct('name', sprintf('sine_%gs', Tm), 'kind', 'sine', ... %#ok<AGROW>
        'mod_period', Tm, 'n_cycles', cfg.sine_cycles(k), 'cycle_len', Tm, ...
        'stim_T', Tm * cfg.sine_cycles(k));
end
if cfg.paired_cycles > 0
    protocols(end + 1) = struct('name', 'paired_pulse', 'kind', 'paired', ...
        'mod_period', NaN, 'n_cycles', cfg.paired_cycles, 'cycle_len', 10, ...
        'stim_T', 10 * cfg.paired_cycles);
end

extra = {};
if ~isempty(cfg.sigma_u_noise), extra = {'sigma_u_noise', cfg.sigma_u_noise}; end

results = [];
t_total = tic;
for p = 1:numel(protocols)
    pr = protocols(p);
    input = stimulus_config(cfg, pr);
    for c = 1:numel(cond_names)
        t_run = tic;
        model = build_from_preset(cfg.preset_name, cond_names{c}, ...
            'T_range', [-cfg.settle_s, pr.stim_T], 'fs', cfg.fs, ...
            'rng_seeds', cfg.rng_seeds, 'input_config', input, ...
            'lya_method', 'none', 'store_full_state', true, ...
            'store_decimated_state', false, 'verbose', cfg.verbose, extra{:});
        params = model.get_params();
        assert(SRNNCellTypePairs.routes_identical(params), ...
            'run_stim_evoked:RoutesDiffer', 'The readouts take route 1 and need identical routes.');
        model.run();
        t_sim = toc(t_run);

        % The stimulus as the analysis needs it (pulse times, driven set),
        % from the same generator and seed the model used.
        if strcmp(pr.kind, 'sine')
            [~, ~, info] = generate_am_pulse_train(params, pr.stim_T, cfg.fs, cfg.rng_seeds(2), input);
        else
            [~, ~, info] = generate_paired_pulse_input(params, pr.stim_T, cfg.fs, cfg.rng_seeds(2), input);
        end
        stim_mask = false(model.n, 1); stim_mask(info.stim_idx) = true;
        assert(isequal(find(any(model.u_ex ~= 0, 2)), info.stim_idx(:)), ...
            'run_stim_evoked:StimSet', 'Driven neurons differ from the generator''s.');

        [t, ro] = population_readouts(model, params, stim_mask, cfg.readout_deci);
        r = struct('condition', cond_names{c}, 'title', titles(cond_names{c}), ...
            'protocol', pr.name, 'kind', pr.kind, 'mod_period', pr.mod_period, ...
            'n_cycles', pr.n_cycles, 'cycle_len', pr.cycle_len, 'stim_T', pr.stim_T, ...
            't', t, 'readouts', ro, 'info', info, ...
            'n_stim', nnz(stim_mask), 'n', model.n, ...
            'cell_type_names', {model.cell_type_names}, ...
            'ode_solver', model.ode_solver, 'runtime_s', toc(t_run), 'sim_s', t_sim);
        if isempty(results), results = r; else, results(end + 1) = r; end %#ok<AGROW>
        vprintf(cfg.verbose, 'minimal', '  %-13s %-14s %5.0f s simulated in %6.1f s\n', ...
            pr.name, cond_names{c}, pr.stim_T + cfg.settle_s, toc(t_run));
        clear model
    end
end

settings = struct('cfg', cfg, 'protocols', protocols, 'cond_names', {cond_names}, ...
    'total_runtime_s', toc(t_total), 'created', char(datetime('now')), ...
    'git', capture_git_provenance(run_dir, project_root));
mat_file = fullfile(out_dir, 'stim_evoked_data.mat');
save(mat_file, 'results', 'settings', '-v7.3');
vprintf(cfg.verbose, 'minimal', '  saved %s (%.1f MB, %.1f min)\n', mat_file, ...
    dir(mat_file).bytes / 2^20, settings.total_runtime_s / 60);
end

function input = stimulus_config(cfg, pr)
% The model's input_config for one protocol. intrinsic_drive 0 as required.
s = cfg.stim;
common = struct('intrinsic_drive', 0, 'amp', s.amp, 'neuron_frac', s.neuron_frac, ...
    'target_types', {s.target_types});
if strcmp(pr.kind, 'sine')
    input = common;
    input.generator    = @generate_am_pulse_train;
    input.width_s      = s.width_s;
    input.pulse_period = s.pulse_period;
    input.mod_period   = pr.mod_period;
    input.mod_depth    = s.mod_depth;
    input.t_first      = 0;
else
    input = common;
    input.generator = @generate_paired_pulse_input;
    input.width_ms  = 1000 * s.width_s;
end
end

function [t, ro] = population_readouts(model, params, stim_mask, deci)
% Population means of x, r, syn, sfa, std per type x {stim, nonstim}, from
% S_out in chunks. Route 1 of each presynaptic type (identical routes).
idx_t = 1:deci:numel(model.t_out);
t = model.t_out(idx_t)';
names = params.cell_type_names;
vars = {'x', 'r', 'syn', 'sfa', 'std'};
groups = {};
for q = 1:numel(names), groups = [groups, {[names{q} '_stim'], [names{q} '_nonstim']}]; end %#ok<AGROW>
ro = struct();
for v = 1:numel(vars)
    for g = 1:numel(groups), ro.(vars{v}).(groups{g}) = nan(1, numel(idx_t), 'single'); end
end
ro.count = struct();
for q = 1:numel(names)
    ix = params.type_indices{q}(:);
    ro.count.([names{q} '_stim']) = nnz(stim_mask(ix));
    ro.count.([names{q} '_nonstim']) = nnz(~stim_mask(ix));
end

chunk = 5000;
for i0 = 1:chunk:numel(idx_t)
    sel = i0:min(i0 + chunk - 1, numel(idx_t));
    [x, a, b, r, so] = SRNNCellTypePairs.unpack_and_compute_states(model.S_out(idx_t(sel), :), params);
    for q = 1:numel(names)
        nm = names{q};
        local = stim_mask(params.type_indices{q}(:));
        nq = numel(local); nts = numel(sel);
        posts = fieldnames(so.(nm));
        sfa = zeros(nq, nts);
        if ~isempty(a.(nm)), sfa = params.c_eff(q) * reshape(sum(a.(nm), 2), nq, nts); end
        dep = ones(nq, nts);
        B = b.(nm).(posts{1});
        if ~isempty(B), dep = reshape(prod(B, 2), nq, nts); end
        data = struct('x', x.(nm), 'r', r.(nm), 'syn', so.(nm).(posts{1}), 'sfa', sfa, 'std', dep);
        masks = {local, ~local};
        tags = {'_stim', '_nonstim'};
        for g = 1:2
            if ~any(masks{g}), continue; end
            for v = 1:numel(vars)
                ro.(vars{v}).([nm tags{g}])(sel) = single(mean(data.(vars{v})(masks{g}, :), 1));
            end
        end
    end
end
end

function tf = is_absolute(p)
tf = ~isempty(regexp(p, '^([A-Za-z]:[\\/]|[\\/])', 'once'));
end
