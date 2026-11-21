function soied_write_sidecar(seeds, fig_seed)
% SOIED_WRITE_SIDECAR Measures for every network run, the main figure for one, and the .md sidecar with the numbers.
%
%   soied_write_sidecar([42 73 101], 42)
%
% Reads data/explorations/stim_outlier_ied/soied_seed<seed>.mat (soied_run),
% computes soied_measures, draws soied_main_figure for fig_seed, writes
% figs/narrative_stim_outlier_ied/soied_main.md and soied_measures.mat.
%
% See also: soied_run, soied_measures, soied_main_figure

root = fileparts(which('setup_paths'));
dd = fullfile(root, 'data', 'explorations', 'stim_outlier_ied');
fd = fullfile(root, 'figs', 'narrative_stim_outlier_ied');
if ~isfolder(fd), mkdir(fd); end
Ms = cell(1, numel(seeds)); Rs = cell(1, numel(seeds));
for k = 1:numel(seeds)
    R = load(fullfile(dd, sprintf('soied_seed%03d.mat', seeds(k))), 'R'); R = R.R;
    Ms{k} = soied_measures(R);
    if seeds(k) == fig_seed
        info = soied_main_figure(R, Ms{k}, fd);
        Sset = R.settings;
    end
    Rs{k} = struct('seed', R.seed, 'struct', R.struct, 'check', R.check, ...
        'n_events', [numel(R.block(1).events.time) numel(R.block(2).events.time)], ...
        'groups_per_event', {cellfun(@numel, R.block(1).events.groups)});
end
save(fullfile(fd, 'soied_measures.mat'), 'Ms', 'Rs', 'info');

fid = fopen(fullfile(fd, 'soied_main.md'), 'w');
c = onCleanup(@() fclose(fid));
p = @(varargin) fprintf(fid, varargin{:});
p('# Discharges, stimulation and the frozen-gain coherent mode in the IED network of Figure 8 (pilot, %s)\n\n', char(datetime('today', 'Format', 'yyyy-MM-dd')));
p('Code: `scripts/explorations/stim_outlier/soied_*.m` (soied_run, then soied_write_sidecar). Data: `data/explorations/stim_outlier_ied/`. Figure: `soied_main.png` (network seed %d).\n\n', fig_seed);
p('## What was run\n\n');
p(['The ied_stim bundle''s network (ied_base_config, run13: n = 350, 245 E / 105 I, five subnets of 70, one SFA variable per neuron, ' ...
   'STD on E outputs, noise 0.025, sra1) at fs = %g Hz, network seeds %s (the bundle''s first three). Per network two blocks, as the bundle''s ' ...
   '''essential'' numbering: block 1 input 0, block 2 constant uniform input %.2f, noise seed base + (seed - 42) + 1000 x block; ' ...
   '10-s model warmup, %g-s settle, then **%g s analysed** (the bundle analyses 300 s; shortened to fit a 30-min budget). ' ...
   'Discharges: ied_detect_events on plot_data (250 Hz) over the analysed window, exactly as ied_stim_blocks. ' ...
   'Event time = the detector''s, the dendritic peak.\n\n'], Sset.fs, mat2str(seeds), Sset.amplitude, Sset.settle, Sset.block_len);
p(['Frozen-gain M = W diag(d) = I + tau_d J_xx (ospec_effective_matrix) every %g s, and every %g s in [%g, %g] s around each discharge; ' ...
   'dense eig with vectors. coh10 = fraction of a unit eigenvector''s norm^2 in the span of the ten population indicators (E, I of each subnet; ' ...
   'random vector 0.029). **Coherent mode** = rightmost mode with coh10 >= %.1f; rc its real part. R_eff = bulk radius predicted from the gains ' ...
   '(2-type blocks, ospec); R10 the same on the 10 population blocks. Full Jacobian (1190 x 1190) every %g s by dense eig; network mode = ' ...
   'eigenvector >= 10%% in the dendritic rows. Around the first 6 baseline discharges and 6 control times (>= 1 s from any discharge) the full ' ...
   'Jacobian every %g s in [-0.6, 0.2] s (soied_event_full; re-simulation reproduces the block: checked).\n\n'], ...
   Sset.coarse_dt, Sset.dense_dt, Sset.dense_win, Ms{1}.thr, Sset.full_dt, 0.05);
p(['ospec_full_spectrum accepts this model class (an SRNNCellTypePairs subclass) but with N = 1190 > 1000 it uses eigs ' ...
   '(''largestreal'', k = 12), which converged for %s of 12 eigenvalues at the probed baseline states, and among those none was a network mode ' ...
   '(they are single-neuron slow modes at the kinetic edge). The network mode was therefore taken from a dense eig (about 0.5 s per state).\n\n'], ...
   strjoin(arrayfun(@(k) num2str(Ms{k}.ospec_full_probe.n_converged), 1:numel(Ms), 'UniformOutput', false), ', '));
p('## Check of M (block 1, t = %.2f s, seed %d)\n\n', Rs{1}.check.t, Rs{1}.seed);
for k = 1:numel(Rs)
    ck = Rs{k}.check;
    p('- seed %d: max|M| = %.3f; max|M - (I + tau_d J_xx)| = %.3g (the route ospec builds M by); max|M - W diag(F phi'')| with F the route output factor and phi'' the slope formed from the route models = %.3g; max|M - (I + tau_d J_xx,FD)| with a central finite difference of dynamics_fast (h = 1e-7) = %.3g.\n', ...
        Rs{k}.seed, ck.max_abs_M, ck.M_vs_I_plus_tau_Jxx, ck.M_vs_W_diag_F_rp, ck.M_vs_I_plus_tau_Jxx_FD);
end
p('\n## Per network and block\n\n');
p('b = baseline (input 0), s = stimulation (0.2). Fractions on the uniform grid. Gains d: mean over neurons (post <- pre). Full J: rightmost network mode, 1-s grid.\n\n');
p('| seed | block | discharges (per min) | frac Re l1(M) > 1 | frac rc > 1 | frac most-coherent > 1 | mean Re l1 | mean rc (p5, p50, p95) | coh10 of l1 (median) | R_eff | R10 | Re K10_1 | d(E<-E) | d(I<-E) | d(E<-I) | rate E | rate I | full J: frac net Re > 0 | net Re median (mean) 1/s | net abs(Im) median | net coh10 median |\n');
p('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(Ms)
    for b = 1:2
        w = Ms{k}.block(b); lab = 'bs';
        p('| %d | %s | %d (%.1f) | %.3f | %.3f | %.3f | %.3f | %.3f (%.2f, %.2f, %.2f) | %.2f | %.3f | %.3f | %.2f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.2f (%.2f) | %.2f | %.2f |\n', ...
            Ms{k}.seed, lab(b), w.n_events, w.events_per_min, w.frac_l1, w.frac_coh, w.frac_mc, w.re_l1, w.re_coh, w.re_coh_p, ...
            w.coh_l1, w.R_eff, w.R10, w.K10_1, w.dbar(1, 1), w.dbar(2, 1), w.dbar(1, 2), w.rate(1), w.rate(2), ...
            w.frac_net_pos, w.net_re, w.net_re_mean, w.net_im, w.net_coh);
    end
end
p('\n## Discharges against the coherent mode (baseline blocks)\n\n');
p(['Strict crossing = upward passage of rc through 1. Excursion = rc above theta, the 95th percentile of rc at baseline quiet times (>= 1 s from any discharge); ' ...
   'passages < 0.3 s apart merged. A crossing or excursion is "free" when no discharge falls in [onset - 0.2, onset + 0.6] s. Lags: onset minus discharge peak. ' ...
   'Attribution: at the max of rc in [-0.3, +0.05] s, the subnet holding the largest share of the coherent eigenvector, and whether the detector marked that subnet.\n\n']);
p('| seed | block | discharges | with rc > 1 in [-0.6, 0.2] | with a strict crossing in [-0.6, 0] | median lag (s) | strict crossings (free) | theta | quiet rc p5/p50/p95 | with excursion | median lag (s) | excursions (free; free per min) | dominant subnet marked | median share on marked subnets | subnets per discharge (median) |\n');
p('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(Ms)
    for b = 1:2
        w = Ms{k}.block(b); lab = 'bs';
        p('| %d | %s | %d | %d | %d | %.3f | %d (%d) | %.3f | %.2f / %.2f / %.2f | %d | %.3f | %d (%d; %.1f) | %d / %d | %.2f | %.0f |\n', ...
            Ms{k}.seed, lab(b), w.n_events, sum(w.ev_above), sum(w.ev_cross), median(w.lag_cross, 'omitnan'), w.n_cross, w.n_cross_free, ...
            Ms{k}.theta, Ms{k}.rc_quiet_p, sum(w.ev_exc), median(w.lag_exc, 'omitnan'), w.n_exc, w.n_exc_free, w.exc_rate_free, ...
            sum(w.dom_in_event == 1), sum(isfinite(w.dom_in_event)), median(w.share_event, 'omitnan'), median(w.n_groups, 'omitnan'));
    end
end
p('\nEvent-triggered median of rc (all baseline discharges, dense grid), lag -0.6 to 0.2 s in steps of 0.05 s:\n\n');
for k = 1:numel(Ms)
    w = Ms{k}.block(1); sel = 1:5:numel(w.trig_lags);
    p('- seed %d: %s\n', Ms{k}.seed, strjoin(arrayfun(@(i) sprintf('%.2f', w.trig_rc(i)), sel, 'UniformOutput', false), ', '));
end
p('\n## Full Jacobian around discharges (6 discharges, 6 control times per network, every 0.05 s)\n\n');
p('| seed | discharges with net Re > 0 in [-0.3, 0] | controls with it | max net Re in [-0.3, 0], discharges (1/s) | same, controls | median profile peaks at (s) | median first > 0 from -0.4 s at (s) | median net coh10 in [-0.3, 0] |\n|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(Ms)
    e = Ms{k}.ef;
    p('| %d | %.2f | %.2f | %s | %s | %.2f | %.2f | %.2f |\n', Ms{k}.seed, e.ev_pos_pre, e.ctrl_pos_pre, ...
        mat2str(e.ev_max_pre, 2), mat2str(e.ctrl_max_pre, 2), e.lag_peak, e.lag_first_pos, e.ev_net_coh_pre);
end
p('\nMedian profiles (lag -0.6 to 0.2 s, step 0.05 s), discharges / controls:\n\n');
for k = 1:numel(Ms)
    e = Ms{k}.ef;
    p('- seed %d: %s\n  / %s\n', Ms{k}.seed, strjoin(compose('%.2f', e.ev_net_med), ', '), strjoin(compose('%.2f', e.ctrl_net_med), ', '));
end
p('\n## Structure\n\n');
for k = 1:numel(Rs)
    s = Rs{k}.struct; evw = sort(real(s.ev_W), 'descend'); kk = sort(real(s.K10_eig), 'descend');
    p('- seed %d: rightmost Re eig(W) %.2f, %.2f, %.2f; structural bulk radius (2-type) %.2f; Re eig of the 10 x 10 block-mean matrix %s.\n', ...
        Rs{k}.seed, evw(1:3), s.R_struct, mat2str(kk', 3));
end
p('\n## Figure states (seed %d)\n\nBaseline window [%g %g] s, stimulation window [%g %g] s. Zoomed discharge at %.3f s (subnet %s). Complex planes: quiet t = %.2f s, before discharge t = %.2f s, stimulation t = %.2f s. Eigenvalues as rates (lambda - 1)/tau_d; ringed: the coherent mode; dashed circle: R_eff around -1/tau_d.\n', ...
    fig_seed, info.window(1, :), info.window(2, :), info.zoom_event, mat2str(info.zoom_group), info.t_quiet, info.t_pre, info.t_stim);
end
