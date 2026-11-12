% P1B_STIM_BURSTING Block B pilot 1b: does tonic stimulation suppress the bursts under chained MTS adaptation?
%
% The bursting network of pilot 0b (preset bursting_pairs, 50 neurons, 35 E
% and 15 I) under its original adaptation and under chained
% multiple-timescale (MTS) adaptation, both built by bursting_condition_model
% (scripts/examples/bursting_chain_mts_example; the one definition of the
% conditions). Seed [42 42], 10 s of settling, then 100 s: no stimulus from
% 0 to 40 s (the baseline), stimulus ON from 40 s to 100 s. Measured before
% (10-40 s) and during (50-100 s, leaving 10 s for the transient).
%
% STIMULUS CASES (the settings block below states each exactly)
%   ref       the manuscript's stimulation figure (fig_stim_engages_adaptation)
%             run as it runs it: dc_staircase_stimulus with dc_levels [0 0.4],
%             i.e. a DC step of 0.4 to EVERY neuron, plus independent white
%             noise per neuron of intensity 0.001 (sd 0.001*sqrt(400) = 0.02
%             per 400 Hz sample) on [0, 100] s, u_ex_scale 1, rk4 at 400 Hz.
%             The figure holds each level 50 s from t = 0 with no settling;
%             here the zero level is held 40 s after 10 s of settling.
%   dc0p4     the same step with no noise (p1b_dc_step)
%   dc0p2     half the amplitude, every neuron, no noise
%   dc0p1     a quarter, every neuron, no noise
%   part0p4   0.4 to 7 of the 35 E neurons (20 %), E neurons 1:5:35
%             (indices 1 6 11 16 21 26 31), no noise, the same in every
%             condition
% Conditions 'original' and 'chain_mts' under every case; 'product_mts'
% under ref only (it does not burst; pilot 0b).
%
% A BURST is pilot 0b's: population mean E rate above 0.2 for at least
% 0.1 s, dips under 0.2 s merged (p1b_detect_bursts: the same arithmetic
% as the local function in bursting_chain_mts_figure.m, which is not
% callable). One rule added here: an above-threshold stretch longer than
% 5 s (ten times the original's median burst length, 0.46 s, rounded) is
% SUSTAINED high firing, not a burst, and is not counted as one; without
% it, a step that lifts the population above 0.2 is one 60 s "burst".
%
% Run setup_paths once in the session first (cwd at the repository root).
% No pool: eleven 110 s simulations, about 3 minutes in all in a MATLAB
% with no desktop. Figures are drawn invisible and closed after saving
% when save_figs is set.
%
% WHAT TO EXPECT with the settings below (block B pilot 1b, on main at
% b051af1). Before the stimulus both conditions burst at 0.33 per s
% (largest population mean E rate 0.59 original, 0.39 chain MTS). During
% (50-100 s), mean E rate / largest mean E rate / peak-to-peak:
%   original   ref 0.297/0.309/0.018  dc0p4 0.297/0.305/0.010
%              dc0p2 0.217/0.222/0.007  dc0p1 0.170/0.175/0.005
%              part0p4 0.185/0.236/0.092: NOT steady, a small 0.88 Hz
%              oscillation crossing the 0.2 threshold (0.88 "bursts"/s)
%   chain_mts  ref 0.200/0.208/0.011  dc0p4 0.200/0.205/0.005
%              dc0p2 0.133/0.136/0.003  dc0p1 0.097/0.098/0.001
%              part0p4 0.120/0.122/0.003
%   product_mts ref 0.170/0.176/0.009 (no bursts before either)
% Every all-neuron step ends the bursts in both conditions and leaves a
% STEADY rate above baseline, well below saturation (no neuron's mean rate
% above 0.9), with the adaptation term 1.6 to 3.6 times its baseline mean and the
% depression factor lower. Chain MTS under ref sits at a mean E rate of
% 0.200, exactly the burst threshold, so the input noise makes it cross
% the threshold 38 times (0.76 "bursts"/s, each peaking at most 0.204):
% read the peak-to-peak, not the bare count, there. The noise-free run
% (dc0p4) counts 0. The model is deterministic (sigma_u_noise 0, rk4;
% ref's input noise is seeded), so a rerun gives these numbers exactly.
%
% SAVED with save_figs, in out_dir:
%   p1b_original.png/.fig, p1b_chain_mts.png/.fig   one condition each,
%       columns = cases, rows = rates, mean rates with bursts, adaptation,
%       depression, input
%   p1b_summary.png/.fig    before against during, per case and condition
%   traces.mat   R.<condition>.<case>: t, mean_r_E, mean_r_I, r_E, r_I
%       (single), mean_adaptation, mean_depression, u_mean, u_target,
%       target, bursts, measures (p1b_measures), values (every value the
%       model was built with); plus S (settings) and stims (the cases).
%
% See also: bursting_condition_model, p1b_run_case, p1b_measures,
%           p1b_dc_step, fig_stim_engages_adaptation

%% Settings
S = struct();
S.rng_seeds  = [42 42];      % the stimulation figure's seed pair
S.T_settle   = 10;           % s simulated before t = 0, not shown
S.T_show     = 100;          % s shown and analysed
S.plot_deci  = 4;            % 400 Hz integration -> 100 Hz traces
S.t_stim     = 40;           % s, stimulus onset
S.win_before = [10 40];      % s, baseline window
S.win_during = [50 100];     % s, stimulus window (10 s transient skipped)
S.pinned_thresh = 0.9;       % a rate above this is near saturation
S.burst = struct('thresh', 0.2, 'min_dur', 0.1, 'merge', 0.2, ...
                 'sustained_max', 5);   % pilot 0b's, plus the 5 s rule above
amp_ref   = 0.4;             % fig_stim_engages_adaptation: dc_levels [0 0.4]
noise_ref = 0.001;           % fig_stim_engages_adaptation: noise_intensity
part_E    = 1:5:35;          % the partially stimulated E neurons (7 of 35)
conds       = {'original', 'chain_mts', 'product_mts'};
cond_labels = {'original (sfa3\_std1)', 'chain MTS', 'product MTS'};
product_cases = {'ref'};     % product MTS runs these cases only
save_figs = true;
out_dir   = 'C:\Users\m218089\Desktop\github_repos\FractionalReservoir\figs\explore_B\p1b_stim_bursting';
visible   = false;

%% The stimulus cases
ic_ref = struct('dc_levels', [0 amp_ref], 'hold_dur', S.t_stim, 'ramp_dur', 0, ...
    'noise_intensity', noise_ref, 'intrinsic_drive', [], 'generator', @dc_staircase_stimulus);
step = @(amp, target) struct('amp', amp, 't_on', S.t_stim, 'target', target, ...
    'intrinsic_drive', [], 'generator', @p1b_dc_step);
stims = struct( ...
    'name',  {'ref', 'dc0p4', 'dc0p2', 'dc0p1', 'part0p4'}, ...
    'short', {'ref (0.4+noise)', 'DC 0.4', 'DC 0.2', 'DC 0.1', '0.4 to 20% E'}, ...
    'label', {sprintf('reference: DC %g all + noise %g', amp_ref, noise_ref), ...
              sprintf('DC %g, all neurons', amp_ref), ...
              sprintf('DC %g, all neurons', amp_ref / 2), ...
              sprintf('DC %g, all neurons', amp_ref / 4), ...
              sprintf('DC %g, E %s only', amp_ref, mat2str(part_E))}, ...
    'target', {[], [], [], [], part_E}, ...
    'input_config', {ic_ref, step(amp_ref, []), step(amp_ref / 2, []), ...
                     step(amp_ref / 4, []), step(amp_ref, part_E)});

%% Run
t0 = tic;
R = struct();
for c = 1:numel(conds)
    R.(conds{c}) = cell(1, numel(stims));
    for k = 1:numel(stims)
        if strcmp(conds{c}, 'product_mts') && ~ismember(stims(k).name, product_cases)
            continue
        end
        R.(conds{c}){k} = p1b_run_case(conds{c}, stims(k), S);
    end
end
run_time = toc(t0);

%% The table of measures
fprintf(['%-12s %-8s | before: bursts/s rE rI maxE p2p fHz ad dep | during: bursts/s rE rI maxE ' ...
    'p2p fHz ad dep | during: above-thr pinned(t,n) n_sust t_last\n'], 'condition', 'case');
for c = 1:numel(conds)
    for k = 1:numel(stims)
        C = R.(conds{c}){k};
        if isempty(C), continue; end
        b = C.measures.before; d = C.measures.during;
        fprintf(['%-12s %-8s | %.2f %.3f %.3f %.3f %.3f %.2f %.3f %.3f | %.2f %.3f %.3f %.3f %.3f %.2f %.3f %.3f | ' ...
            '%.2f %.3f %.2f %d %.2f\n'], conds{c}, stims(k).name, ...
            b.burst_rate, b.mean_r_E, b.mean_r_I, b.max_mean_r_E, b.p2p_mean_r_E, b.dominant_freq, ...
            b.mean_adaptation, b.mean_depression, ...
            d.burst_rate, d.mean_r_E, d.mean_r_I, d.max_mean_r_E, d.p2p_mean_r_E, d.dominant_freq, ...
            d.mean_adaptation, d.mean_depression, ...
            d.frac_above, d.frac_pinned_t, d.frac_pinned_n, d.n_sustained, C.measures.t_last_burst);
    end
end
fprintf('run time %.1f s\n', run_time);

%% Figures
F = struct();
for c = 1:2
    F.(conds{c}) = p1b_condition_figure(R.(conds{c}), S, sprintf(['%s: bursting_pairs, seed %s, ' ...
        'stimulus from %g s (red dashed). Burst: mean E rate > %.2f for >= %.2f s (gaps < %.2f s ' ...
        'merged), grey; above threshold > %g s = sustained, orange.'], ...
        strrep(cond_labels{c}, '\_', '_'), mat2str(S.rng_seeds), S.t_stim, S.burst.thresh, ...
        S.burst.min_dur, S.burst.merge, S.burst.sustained_max), visible);
end
F.summary = p1b_summary_figure(R, stims, conds, cond_labels, S, visible);

%% Save
if save_figs
    if ~isfolder(out_dir), mkdir(out_dir); end
    tags = {'original', 'p1b_original'; 'chain_mts', 'p1b_chain_mts'; 'summary', 'p1b_summary'};
    for i = 1:size(tags, 1)
        fig = F.(tags{i, 1});
        exportgraphics(fig, fullfile(out_dir, [tags{i, 2} '.png']), 'Resolution', 110);
        old_fcn = fig.CreateFcn;
        fig.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';   % reopen visible
        savefig(fig, fullfile(out_dir, [tags{i, 2} '.fig']));
        fig.CreateFcn = old_fcn;
    end
    Rs = struct();
    for c = 1:numel(conds)
        for k = 1:numel(stims)
            if ~isempty(R.(conds{c}){k})
                Rs.(conds{c}).(stims(k).name) = R.(conds{c}){k};
            end
        end
    end
    traces = struct('R', Rs, 'S', S, 'stims', stims, 'run_time_s', run_time); %#ok<NASGU>
    save(fullfile(out_dir, 'traces.mat'), '-struct', 'traces');
    if ~visible, cellfun(@close, struct2cell(F)); end
end
