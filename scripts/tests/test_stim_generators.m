% TEST_STIM_GENERATORS The AM pulse train and the paired-pulse protocol produce what the experiment specifies.
%
% Assert-based and silent until it throws. No simulation: a minimal params
% struct (two cell types) stands in for a model.
%
%   generate_am_pulse_train: pulse count and period, envelope range 0.7-1.0
%       of the peak with the first pulse at the trough and the peak half a
%       modulation period later, width in seconds, charge per pulse under
%       the model's linear interpolant, the driven set E-only and fixed.
%   generate_paired_pulse_input: event times 1.25/3.75/6.25/8.75 s of each
%       10-s cycle with ISIs 0/2/10/100 ms, the default draw unchanged
%       (same neurons as rng(seed); randperm(n, round(n*frac))), and the
%       same driven set as the AM train for one seed and config.

n = 40; fs = 2000;
params = struct('n', n, 'cell_type_names', {{'E', 'I'}}, ...
    'type_indices', {{1:20, 21:40}});

%% AM pulse train
cfg = struct('amp', 40, 'width_s', 1e-3, 'pulse_period', 0.2, 'mod_period', 2, ...
    'mod_depth', 0.3, 'neuron_frac', 0.5, 'target_types', {{'E'}});
T = 10;
[u, t, info] = generate_am_pulse_train(params, T, fs, 7, cfg);
assert(isequal(size(u), [n, numel(t)]) && abs(t(end) - T) < 1e-12 && t(1) == 0);
% 5 Hz over [0, 10] s: pulses at 0, 0.2, ..., 9.8 (the one at 10 s does not fit its width)
assert(numel(info.pulse_times) == 50, 'pulse count');
assert(all(abs(diff(info.pulse_times) - 0.2) < 1e-12), 'pulse period');
rel = info.pulse_amps / cfg.amp;
assert(abs(min(rel) - 0.7) < 1e-12 && abs(max(rel) - 1.0) < 1e-12, 'envelope range');
assert(abs(rel(1) - 0.7) < 1e-12, 'first pulse at the trough');
[~, ipk] = max(rel);
assert(abs(info.pulse_times(ipk) - 1.0) < 1e-12, 'peak half a period after the first pulse');
assert(abs(mean(rel) - 0.85) < 1e-12, 'mean 0.85 of the peak over whole cycles');
% Driven set: half of E, nobody in I, identical rows
stim = find(any(u ~= 0, 2))';
assert(isequal(stim, info.stim_idx) && numel(stim) == 10 && all(stim <= 20), 'driven set');
assert(all(all(u(stim, :) == u(stim(1), :))), 'same input on every driven neuron');
% Width: 1 ms = 2 samples; charge per pulse = amp * width under linear interpolation
assert(info.width_samples == 2);
row = u(stim(1), :);
k1 = round(info.pulse_times(1) * fs) + 1;
assert(isequal(row(k1:k1 + 1), [28 28]) && row(k1 + 2) == 0);
tt = 0.2 + linspace(-0.01, 0.01, 200001);   % the second pulse: the first, at t = 0, has no ramp before it
q = trapz(tt, interp1(t, row, tt, 'linear', 0));
assert(abs(q - info.pulse_amps(2) * 1e-3) < 1e-6, 'pulse charge under linear interpolation');
% Total count of nonzero runs equals the pulse count
assert(nnz(diff([0, row ~= 0]) == 1) == 50);
% Envelope convention with a 5-s period: trough at 0, peak at 2.5 s, which falls
% between pulses (2.4 and 2.6 s carry equal amplitudes, the largest)
cfg5 = cfg; cfg5.mod_period = 5;
[~, ~, i5] = generate_am_pulse_train(params, 20, fs, 7, cfg5);
[~, ipk] = max(i5.pulse_amps);
assert(abs(i5.pulse_times(ipk) - 2.4) < 1e-12 && abs(i5.pulse_amps(ipk) - i5.pulse_amps(ipk + 1)) < 1e-9 && abs(i5.pulse_amps(1) - 28) < 1e-12);
assert(isequal(i5.stim_idx, info.stim_idx), 'driven set independent of the period');
% The global RNG is left alone
rng(123); a1 = rand; rng(123); generate_am_pulse_train(params, 1, fs, 7, cfg); a2 = rand;
assert(a1 == a2, 'global RNG untouched');

%% Paired pulse
pcfg = struct('amp', 40, 'width_ms', 1, 'neuron_frac', 0.5, 'target_types', {{'E'}});
[up, tp, pinfo] = generate_paired_pulse_input(params, 30, fs, 7, pcfg);
assert(isequal(pinfo.stim_idx, info.stim_idx), 'same driven set as the AM train');
row = up(pinfo.stim_idx(1), :);
starts = tp(diff([0, row ~= 0]) == 1)';
expect = [];
ev = [1.25 3.75 6.25 8.75]; isis = [0 0.002 0.010 0.100];
for c = 0:2
    for e = 1:4
        expect = [expect, c * 10 + ev(e)]; %#ok<AGROW>
        if isis(e) > 0, expect = [expect, c * 10 + ev(e) + isis(e)]; end %#ok<AGROW>
    end
end
assert(numel(starts) == numel(expect) && max(abs(starts - expect)) < 0.5 / fs, 'paired-pulse event times');
assert(isequal(pinfo.onset_times, reshape(ev' + (0:2) * 10, 1, [])), 'onset times');
assert(isequal(pinfo.onset_isi, repmat(isis, 1, 3)), 'onset ISIs');
assert(all(row(row ~= 0) == 40), 'paired-pulse amplitude');
% Default draw (no target_types) unchanged from the pre-2026-10-20 generator
[ud, ~, dinfo] = generate_paired_pulse_input(params, 10, fs, 5, struct('neuron_frac', 0.25));
rng(5); old = sort(randperm(n, round(n * 0.25)));
assert(isequal(dinfo.stim_idx, old) && isequal(find(any(ud ~= 0, 2))', old), 'default draw unchanged');
