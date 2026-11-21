% RUN_STIM_OUTLIER Does the population mode of the frozen-gain connectivity cross Re = 1 at each burst, and stay inside under stimulation?
%
% A pilot (2026-11-21) connecting two results: sustained input suppresses
% the bursts of the bursting network (fig_stim_engages_adaptation; block B
% pilot 1b), and adaptation pulls a structural outlier eigenvalue of the
% connectivity back toward the bulk (the outlier_spectra stage). Hypothesis
% (StochasticPlasticDynamicalSystemPaper, docs/future_directions/
% adaptation_reins_in_outlier_eigenvalues, section 1 "Bursting network" and
% prediction 5): a population burst is the population-coherent mode of
% M(t) = W diag(d(t)) crossing Re = 1; under stimulation it stays inside.
%
% THE NETWORK, STIMULUS AND CONDITIONS are pilot 1b's (p1b_stim_bursting),
% case dc0p4: bursting_pairs (n = 50, 35 E / 15 I, indegree 10), condition
% 'original' (the preset's sfa3_std1: SFA on E, tau_a 0.25/1.58/10 s,
% c = 0.5; STD on E->E and E->I, tau_rec 1 s, tau_rel 0.25 s), seed
% [42 42] (fig_stim_engages_adaptation's), 10 s settling, then a DC step
% of 0.4 to every neuron from 40 s to 100 s, no noise, rk4 at 400 Hz. Also
% 'chain_mts' on the same seed. Further weight seeds [s s], s = 1..8,
% 'original'. Baseline window 10-40 s, stimulus window 50-100 s (pilot
% 1b's: 10 s of transient skipped). Bursts: pilot 1b's rule (mean E rate
% > 0.2 for >= 0.1 s, gaps < 0.2 s merged).
%
% FROZEN CONTROL: from the last sample before the step, every state but x
% held (p2_frozen_rhs, pilot 2's method), the same step applied.
%
% Parts (set so_parts before running to do only some):
%   'main'     seed [42 42], original (with full Jacobian and frozen arm)
%              and chain_mts (with frozen arm)
%   'seeds'    seeds 1..8, original, with frozen arm
%   'fullJ'    the full Jacobian's spectrum at five states during the
%              step, seeds 42, 1, 6, 3 (about 1 min)
%   'figures'  so_main_figure and so_seeds_figure from the saved data,
%              and the sidecar (so_write_sidecar)
% Data: data/explorations/stim_outlier/; figures: figs/narrative_stim_outlier/.
% Run setup_paths first (cwd at the repo root). No pool. On the PC
% (2026-11-21): main 157 s, seeds 472 s.
%
% See also: so_run_network, so_crossings, so_main_figure, so_seeds_figure,
%           so_write_sidecar, p1b_stim_bursting, run_p2_frozen_state

if ~exist('so_parts', 'var'); so_parts = {'main', 'seeds', 'fullJ', 'figures'}; end
root     = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
data_dir = fullfile(root, 'data', 'explorations', 'stim_outlier');
fig_dir  = fullfile(root, 'figs', 'narrative_stim_outlier');
if ~isfolder(data_dir); mkdir(data_dir); end
if ~isfolder(fig_dir);  mkdir(fig_dir);  end

S = struct();
S.T_settle   = 10;
S.T_show     = 100;
S.t_stim     = 40;
S.amp        = 0.4;          % fig_stim_engages_adaptation's DC level
S.burst      = struct('thresh', 0.2, 'min_dur', 0.1, 'merge', 0.2, 'sustained_max', 5);
S.win_before = [10 40];
S.win_during = [50 100];
S.spec_every = 4;            % 400 Hz -> a spectrum every 10 ms
S.full_every = 40;           % full Jacobian every 100 ms
S.frozen_every = 8;          % frozen arm every 20 ms
S.tol_pre    = 0.6;          % s: a crossing up to this long before a burst's onset coincides with it
S.tol_post   = 0.1;
main_seed  = [42 42];
extra_seeds = 1:8;

if any(strcmp(so_parts, 'main'))
    t0 = tic;
    S.do_full = true;  S.do_frozen = true;
    R_orig = so_run_network('original', main_seed, S);
    S.do_full = false;
    R_chain = so_run_network('chain_mts', main_seed, S);
    main_s = toc(t0);
    save(fullfile(data_dir, 'stim_outlier_main.mat'), 'R_orig', 'R_chain', 'S', 'main_s', '-v7.3');
    fprintf('main done in %.0f s\n', main_s);
end

if any(strcmp(so_parts, 'seeds'))
    t0 = tic;
    S.do_full = false;  S.do_frozen = true;
    R_seeds = cell(1, numel(extra_seeds));
    for i = 1:numel(extra_seeds)
        Ri = so_run_network('original', extra_seeds(i) * [1 1], S);
        Ri.intact = rmfield(Ri.intact, 'ev');    % keep the file small
        Ri.frozen = rmfield(Ri.frozen, 'ev');
        R_seeds{i} = Ri;
        fprintf('seed %d: %.0f s\n', extra_seeds(i), toc(t0));
    end
    seeds_s = toc(t0);
    save(fullfile(data_dir, 'stim_outlier_seeds.mat'), 'R_seeds', 'extra_seeds', 'S', 'seeds_s', '-v7.3');
end

if any(strcmp(so_parts, 'fullJ'))
    % All eigenvalues of the full Jacobian at five states during the step,
    % on the figure's seed and on the seeds whose frozen-gain mode sits at
    % or beyond Re = 1 under the step. Columns of FJ.s<seed>: t | Re, Im
    % of the rightmost network mode (1/s) | its coherence | abscissa |
    % Re lambda_1(M) | its coherence.
    FJ = struct();
    ic = struct('amp', S.amp, 't_on', S.t_stim, 'target', [], 'intrinsic_drive', [], 'generator', @p1b_dc_step);
    for sd = [42 1 6 3]
        m = bursting_condition_model('original', 'rng_seeds', [sd sd], 'T_settle', S.T_settle, ...
            'T_show', S.T_show, 'input_config', ic, 'plot_deci', 4, 'model_args', {'store_full_state', true});
        m.run();
        params = m.get_params(); ws = ospec_weight_stats(params); tq = m.t_out;
        tt = [45 60 75 90 99]; X = zeros(numel(tt), 7);
        for q = 1:numel(tt)
            Sx = m.S_out(find(tq >= tt(q), 1), :)';
            f = ospec_full_spectrum(Sx, params, ws, 223);
            s = ospec_snapshot(Sx, params, ws, 3, 50);
            X(q, :) = [tt(q), real(f.net), imag(f.net), f.coh(find(f.xfrac >= 0.1, 1)), ...
                f.abscissa, real(s.ev_top(1)), s.coh(1)];
        end
        FJ.(sprintf('s%d', sd)) = X;
    end
    save(fullfile(data_dir, 'stim_outlier_fullJ_check.mat'), 'FJ');
end

if any(strcmp(so_parts, 'figures'))
    D = load(fullfile(data_dir, 'stim_outlier_main.mat'));
    so_main_figure(D.R_orig, D.S, fullfile(fig_dir, 'stim_outlier_main'));
    E = load(fullfile(data_dir, 'stim_outlier_seeds.mat'));
    so_seeds_figure([{D.R_orig}, E.R_seeds], D.S, fullfile(fig_dir, 'stim_outlier_seeds'));
    G = load(fullfile(data_dir, 'stim_outlier_fullJ_check.mat'));
    so_write_sidecar(D, E, fullfile(fig_dir, 'stim_outlier.md'), G.FJ);
end
