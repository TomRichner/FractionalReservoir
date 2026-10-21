function M = stim_evoked_measures(results, opts)
% STIM_EVOKED_MEASURES Evoked-response size, sine-fit lead and paired-pulse ratios.
%
%   M = STIM_EVOKED_MEASURES(results)
%   M = STIM_EVOKED_MEASURES(results, 'discard_s', 10, 'readouts', R)
%
% Reads the population readouts run_stim_evoked saved and applies the human
% analysis of brain-stim-analysis to them, per condition x protocol:
%
% EPOCHS. Around every pulse onset, -20 to +190 ms of a readout, minus its
% -20 to -2 ms mean. Only pulses at or after discard_s (the build-up of slow
% adaptation once the train starts) and whose epoch fits are used.
%
% SINE TRAINS. The epochs are averaged over modulation cycles at each pulse
% position in the cycle; the response SIZE at a position is the RMS of that
% average over 10 to 60 ms. A sinusoid at the modulation frequency is fitted
% to size against stimulus phase phi (0 = envelope trough, the first pulse;
% pi = envelope peak) by least squares, size = a cos(phi) + b sin(phi) + c.
% LEAD = pi - atan2(b, a) in degrees, wrapped to (-180, 180]: positive means
% the response peaks BEFORE the envelope. Also amplitude hypot(a, b), offset
% c, relative modulation hypot(a,b)/c and r2. Cycle-averaged continuous
% traces of every readout over one modulation cycle go with it.
%
% PAIRED PULSES. Epochs (-50 to +400 ms here, same baseline) averaged per
% event type over cycles. Ratio for ISI d: RMS of the pair's average over
% [d + 10, d + 60] ms divided by RMS of the single pulse's over [10, 60] ms
% (the human definition, so the first pulse's tail is in the numerator at
% short ISIs). The LTI prediction is the same ratio computed on
% s(t) + s(t - d), s the single-pulse average.
%
% Options:
%   discard_s - stimulation time excluded from the averages (default 10)
%   readouts  - N x 3 cell {name, var, groups}: var a field of
%               results(k).readouts ('x','r','syn','sfa','std'), groups a cell
%               of group names whose count-weighted mean is the signal.
%               The FIRST row is the primary readout (M(k).primary).
%               Default: the CCEP-like readouts listed in default_readouts.
%
% Output M: struct array parallel to results, with condition, protocol,
% kind, mod_period, and per readout (M(k).fits(i), M(k).sizes{i}, ...) for
% sine; M(k).paired for paired pulse. See the field assignments below.
%
% See also: run_stim_evoked, fig_sine_paired_stim

arguments
    results struct
    opts.discard_s (1,1) double = 10
    opts.readouts cell = default_readouts()
end
M = [];
for k = 1:numel(results)
    R = results(k);
    m = struct('condition', R.condition, 'title', R.title, 'protocol', R.protocol, ...
        'kind', R.kind, 'mod_period', R.mod_period, 'readout_names', {opts.readouts(:, 1)'}, ...
        'fits', [], 'sizes', {{}}, 'phase', [], 't_in_cycle', [], 'n_cycles_used', NaN, ...
        'epoch_t', [], 'epochs', {{}}, 'cycle_t', [], 'cycle_traces', struct(), ...
        'envelope_t', [], 'envelope', [], 'paired', struct());
    fs_r = 1 / median(diff(R.t));
    if strcmp(R.kind, 'sine')
        m = sine_measures(m, R, fs_r, opts);
    else
        m = paired_measures(m, R, fs_r, opts);
    end
    if isempty(M), M = m; else, M(end + 1) = m; end %#ok<AGROW>
end
end

function R = default_readouts()
% name, var, groups. The first is the primary "CCEP-like" readout: the mean
% dendritic potential of the UNDRIVEN neurons, i.e. the response at a
% recording site that receives the pulse only through the network.
R = { ...
    'x_nonstim',   'x',   {'E_nonstim', 'I_nonstim'}; ...
    'x_all',       'x',   {'E_stim', 'E_nonstim', 'I_stim', 'I_nonstim'}; ...
    'x_E_nonstim', 'x',   {'E_nonstim'}; ...
    'x_I_nonstim', 'x',   {'I_nonstim'}; ...
    'r_nonstim',   'r',   {'E_nonstim', 'I_nonstim'}; ...
    'r_stim',      'r',   {'E_stim'}; ...
    'syn_stim',    'syn', {'E_stim'}};
end

function y = signal(R, var, groups)
% Count-weighted mean of the listed groups (groups with no neurons skipped).
num = 0; den = 0;
for g = 1:numel(groups)
    n = R.readouts.count.(groups{g});
    if n == 0, continue; end
    num = num + n * double(R.readouts.(var).(groups{g}));
    den = den + n;
end
y = num / den;
end

function [E, tl] = epochs(y, t, onsets, fs_r, win, base)
lag = round(win(1) * fs_r):round(win(2) * fs_r);
tl = lag / fs_r;
E = nan(numel(onsets), numel(lag));
for i = 1:numel(onsets)
    i0 = round((onsets(i) - t(1)) * fs_r) + 1;
    ix = i0 + lag;
    if ix(1) >= 1 && ix(end) <= numel(y)
        E(i, :) = y(ix);
    end
end
bm = tl >= base(1) & tl <= base(2);
E = E - mean(E(:, bm), 2);
end

function v = rms_win(e, tl, w)
v = sqrt(mean(e(:, tl >= w(1) & tl <= w(2)) .^ 2, 2));
end

function m = sine_measures(m, R, fs_r, opts)
Tm = R.mod_period;
pt = R.info.pulse_times;
cyc = floor(pt / Tm + 1e-9);
pos_t = pt - cyc * Tm;                      % time in cycle
use = pt >= opts.discard_s - 1e-9 & pt + 0.19 <= R.t(end);
% Keep complete cycles only: every position present.
n_pos = round(Tm / R.info.pulse_period);
cyc_ok = unique(cyc(use));
cyc_ok = cyc_ok(arrayfun(@(c) nnz(use & cyc == c) == n_pos, cyc_ok));
use = use & ismember(cyc, cyc_ok);
pos_vals = (0:n_pos - 1) * R.info.pulse_period;
m.t_in_cycle = pos_vals;
m.phase = 2 * pi * pos_vals / Tm;
m.n_cycles_used = numel(cyc_ok);
m.envelope_t = linspace(0, Tm, 400);
m.envelope = R.info.amp * (1 - (R.info.mod_depth / 2) * (1 + cos(2 * pi * m.envelope_t / Tm)));

fits = [];
for i = 1:size(opts.readouts, 1)
    y = signal(R, opts.readouts{i, 2}, opts.readouts{i, 3});
    [E, tl] = epochs(y, R.t, pt(use), fs_r, [-0.02 0.19], [-0.02 -0.002]);
    p_use = pos_t(use);
    avg = nan(n_pos, numel(tl));
    for j = 1:n_pos
        avg(j, :) = mean(E(abs(p_use - pos_vals(j)) < 1e-6, :), 1);
    end
    sz = rms_win(avg, tl, [0.01 0.06]);
    f = sine_fit(m.phase(:), sz(:));
    if isempty(fits), fits = f; else, fits(end + 1) = f; end %#ok<AGROW>
    m.sizes{i} = sz(:)';
    m.epochs{i} = avg;
    m.epoch_t = tl;
end
m.fits = fits;

% Cycle-averaged continuous traces of every readout, over the cycles used.
n_s = round(Tm * fs_r);
m.cycle_t = (0:n_s - 1) / fs_r;
vars = fieldnames(R.readouts); vars = setdiff(vars, {'count'}, 'stable');
for v = 1:numel(vars)
    groups = fieldnames(R.readouts.(vars{v}));
    for g = 1:numel(groups)
        y = double(R.readouts.(vars{v}).(groups{g}));
        acc = zeros(1, n_s); n = 0;
        for c = cyc_ok(:)'
            i0 = round((c * Tm - R.t(1)) * fs_r) + 1;
            if i0 >= 1 && i0 + n_s - 1 <= numel(y)
                acc = acc + y(i0:i0 + n_s - 1); n = n + 1;
            end
        end
        m.cycle_traces.(vars{v}).(groups{g}) = acc / max(n, 1);
    end
end
end

function f = sine_fit(phi, y)
X = [cos(phi), sin(phi), ones(size(phi))];
beta = X \ y;
yh = X * beta;
r2 = 1 - sum((y - yh) .^ 2) / sum((y - mean(y)) .^ 2);
peak = atan2(beta(2), beta(1));
lead = mod(rad2deg(pi - peak) + 180, 360) - 180;
f = struct('a', beta(1), 'b', beta(2), 'offset', beta(3), ...
    'amplitude', hypot(beta(1), beta(2)), 'rel_mod', hypot(beta(1), beta(2)) / beta(3), ...
    'peak_phase', mod(peak, 2 * pi), 'lead_deg', lead, 'r2', r2);
end

function m = paired_measures(m, R, fs_r, opts)
on = R.info.onset_times; isi = R.info.onset_isi;
use = on >= opts.discard_s - 1e-9;
isis = R.info.isis;
P = struct('isis', isis, 'readouts', {opts.readouts(:, 1)'});
P.ratio = nan(size(opts.readouts, 1), numel(isis));
P.ratio_lti = P.ratio;
P.n_events = arrayfun(@(d) nnz(use & abs(isi - d) < 1e-9), isis);
for i = 1:size(opts.readouts, 1)
    y = signal(R, opts.readouts{i, 2}, opts.readouts{i, 3});
    [E, tl] = epochs(y, R.t, on(use), fs_r, [-0.05 0.4], [-0.02 -0.002]);
    iu = isi(use);
    avg = nan(numel(isis), numel(tl));
    for j = 1:numel(isis)
        avg(j, :) = mean(E(abs(iu - isis(j)) < 1e-9, :), 1);
    end
    s = avg(1, :);
    s_ref = rms_win(s, tl, [0.01 0.06]);
    lti = nan(numel(isis), numel(tl));
    for j = 1:numel(isis)
        d = isis(j);
        sh = round(d * fs_r);
        shifted = [zeros(1, sh), s(1:end - sh)];
        lti(j, :) = s + (d > 0) * shifted;
        if d == 0, lti(j, :) = s; end
        P.ratio(i, j) = rms_win(avg(j, :), tl, [d + 0.01, d + 0.06]) / s_ref;
        P.ratio_lti(i, j) = rms_win(lti(j, :), tl, [d + 0.01, d + 0.06]) / s_ref;
    end
    m.epochs{i} = avg;
    P.lti{i} = lti;
end
m.epoch_t = tl;
m.paired = P;
m.n_cycles_used = numel(unique(floor(on(use) / R.info.cycle_len)));
end
