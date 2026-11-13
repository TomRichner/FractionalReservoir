function out = p1a_run_pair(condition, seed, amplitudes, duration, late_start, keep_traces)
%P1A_RUN_PAIR One condition and network seed of pilot 1a: the OFF arm and one ON arm per amplitude.
%
%   out = p1a_run_pair(condition, seed, amplitudes, duration, late_start, keep_traces)
%
% condition   'original' | 'chain_mts' (run13_condition_model)
% seed        network seed; the noise seed follows it (run13_condition_model's
%             default, the offset rule of ied_stim_experiment)
% amplitudes  ON-arm amplitudes of the uniform DC input
% duration    s observed after the 10-s unstimulated settling interval; every
%             arm has this same length, so all arms share the network, the
%             setpoints, the SFA ladders, x0 and the Wiener path
% late_start  s; a secondary window [late_start duration] that drops the
%             stimulus onset (ied_stim_experiment discarded 10 s per half)
% keep_traces true: return population traces of every arm (for the example)
%
% The input is ied_stim_experiment's: dc_staircase_stimulus, uniform over all
% 350 neurons, no ramp, no added input noise; zero during the 10-s settling.
% OFF arm: level 0 for the whole run. ON arm: level A from t = 0 to the end.
%
% out.rows   table, one row per arm: condition, seed, amplitude (0 = OFF),
%            the window measures of p1a_window_measures over [0 duration]
%            (no suffix) and over [late_start duration] (suffix _late), and
%            events_per_min_fixed: the OFF arm's excursion criterion
%            re-centred on the arm's own median (ied_stim_experiment's
%            fixed-criterion control; equal to events_per_min for OFF)
% out.x0_diff  max |x(t=0)| difference of each ON arm from the OFF arm (should
%            be 0: the arms differ only from t = 0 on)
% out.traces   struct array (keep_traces), one per arm
% out.seconds  wall time

% The model seeds the Wiener path and the per-neuron setpoints with
% rng(seed), which keeps the session's generator TYPE: 'twister' in a client
% MATLAB, 'threefry' on a pool worker, so a worker would build a different
% network and noise path from the same seeds. Pin the type to the client's,
% as the documented experiment (run in the client) had it.
prior_rng = rng; rng(0, 'twister'); restore_rng = onCleanup(@() rng(prior_rng)); %#ok<NASGU>
started = tic;
levels = [0 amplitudes(:)'];
rows = table(); traces = struct([]); x0_off = []; base_excess = []; x0_diff = zeros(1, numel(amplitudes));
for k = 1:numel(levels)
    A = levels(k);
    ic = struct('dc_levels', A, 'hold_dur', duration, 'ramp_dur', 0, ...
        'noise_intensity', 0, 'intrinsic_drive', [], 'generator', @dc_staircase_stimulus);
    [m, groups, info] = run13_condition_model(condition, 'duration', duration, ...
        'seed', seed, 'input_config', ic);
    m.run();
    s = p1a_signals(m, info, groups);
    assert(all(s.u == A), 'The input is not the uniform level %g over the observed window.', A);
    [w, det] = p1a_window_measures(s, [0 duration]);
    wl = p1a_window_measures(s, [late_start duration]);
    if k == 1
        x0_off = s.x(:, 1); base_excess = det.threshold - det.baseline;
        w.events_per_min_fixed = w.events_per_min;
    else
        x0_diff(k - 1) = max(abs(s.x(:, 1) - x0_off));
        wf = p1a_window_measures(s, [0 duration], det.baseline + base_excess);
        w.events_per_min_fixed = wf.events_per_min;
    end
    row = struct2table(w);
    late = struct2table(wl);
    late.Properties.VariableNames = strcat(late.Properties.VariableNames, '_late');
    row = [table(string(condition), seed, info.cfg.noise_seed, A, ...
        'VariableNames', {'condition', 'seed', 'noise_seed', 'amplitude'}) row late]; %#ok<AGROW>
    rows = [rows; row]; %#ok<AGROW>
    if keep_traces
        [~, ~, events] = p1a_window_measures(s, [0 duration]);
        det_full = ied_detect_events(s.t, s.x, s.r, s.groups, struct());
        tr = struct('condition', condition, 'seed', seed, 'amplitude', A, 't', s.t, ...
            'u', single(s.u), 'module_x', single(s.module_x), ...
            'rate_E_mean', single(s.rate_E), 'rate_I_mean', single(s.rate_I), ...
            'sfa_E_mean', single(s.sfa_E), 'sfa_I_mean', single(s.sfa_I), ...
            'depression_EE_mean', single(s.dep), ...
            'group_recruitment', single(det_full.group_recruitment), ...
            'events', events, 'params', info.params);
        traces = [traces tr]; %#ok<AGROW>
    end
end
out = struct('rows', rows, 'x0_diff', x0_diff, 'traces', traces, 'seconds', toc(started));
end
