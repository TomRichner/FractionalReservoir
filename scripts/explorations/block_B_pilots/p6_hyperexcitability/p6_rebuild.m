function p6_rebuild(mode, idx_list)
% P6_REBUILD Block B pilot 6, question 3: rebuild sampled networks and run them.
%
%   p6_rebuild('verify', idx_list)  the sweep's own settings (20 s, top-K
%       Lyapunov), to check that a rebuilt network reproduces its stored
%       lambda_1 and mean rate. Writes traces/verify_<idx>_<cond>.mat.
%   p6_rebuild('show', idx_list)    30 s, no Lyapunov estimate, rates kept at
%       50 Hz, for the time-series figure. Writes traces/show_<idx>_<cond>.mat.
%
% The network is reconstructed by ParamSpaceAnalysis2.rebuild_model from the
% result's own network_seed, model_defaults, condition and grid values, which
% is the call the sweep made. Serial in the client; one call per one or two
% networks keeps each MATLAB call to a few minutes.

%
% THE GENERATOR. The sweep ran on parallel-pool workers, whose default random
% generator is Threefry; a client session's is the Mersenne twister. The model
% seeds W, the noise and the stimulus with rng(seed), which keeps whatever
% generator is current, so the same seed draws a DIFFERENT network in the
% client. Found here: a client rebuild missed the stored values (rate 0.115
% stored against 0.127, lambda_1 -0.10 against +0.04, network 761 chained
% MTS); with the global generator set to Threefry first, all 12 rebuilt runs
% match the stored lambda_1 and rate to round-off. So this function switches
% the generator and restores it on exit.

P = p6_config();
s0 = rng;
restore = onCleanup(@() rng(s0));
rng(0, 'threefry');
tr_dir = fullfile(P.out_dir, 'traces');
if ~exist(tr_dir, 'dir'), mkdir(tr_dir); end
d = dir(fullfile(P.run_dir, 'param_space_*'));
psa = ParamSpaceAnalysis2.from_dir(fullfile(d.folder, d.name));
psa.load_condition_results();

for idx = idx_list(:)'
    for c = 1:numel(P.rebuild_conds)
        cond = P.rebuild_conds{c};
        res = psa.results.(cond){idx};
        m = psa.rebuild_model(res);
        m.verbose = 'near-none';
        t0 = tic;
        switch mode
            case 'verify'
                m.build(); m.run();
                S = m.lya_summary();
                out = struct('config_idx', idx, 'cond', cond, ...
                    'LLE_stored', res.LLE, 'LLE_rebuilt', S.LLE, ...
                    'rate_stored', res.mean_rate, ...
                    'rate_rebuilt', mean([m.plot_data.r.E(:); m.plot_data.r.I(:)], 'omitnan'), ... % the sweep's pooled mean
                    'npos_stored', res.n_positive, 'npos_rebuilt', S.n_positive, ...
                    'hks_stored', res.h_KS_bits, 'hks_rebuilt', S.h_KS_bits, ...
                    'seconds', toc(t0));
                fprintf('%5d %-14s lambda_1 %8.4f -> %8.4f (diff %.2g)   rate %.5f -> %.5f (diff %.2g)   %.0f s\n', ...
                    idx, cond, out.LLE_stored, out.LLE_rebuilt, out.LLE_rebuilt - out.LLE_stored, ...
                    out.rate_stored, out.rate_rebuilt, out.rate_rebuilt - out.rate_stored, out.seconds);
            case 'show'
                m.T_range = [0 30];
                m.lya_method = 'none';
                m.plot_freq = 50;
                m.build(); m.run();
                rE = single(m.plot_data.r.E); rI = single(m.plot_data.r.I);
                out = struct('config_idx', idx, 'cond', cond, 't', m.plot_data.t(:)', ...
                    'rE', rE, 'rI', rI, 'mean_rE', mean(rE, 1), 'mean_rI', mean(rI, 1), ...
                    'LLE_stored', res.LLE, 'rate_stored', res.mean_rate, 'seconds', toc(t0));
                fprintf('%5d %-14s 30 s run, %d E neurons, mean E rate %.3f, %.0f s\n', ...
                    idx, cond, size(rE, 1), mean(rE(:)), out.seconds);
            otherwise
                error('p6_rebuild:mode', 'mode is verify or show');
        end
        save(fullfile(tr_dir, sprintf('%s_%d_%s.mat', mode, idx, cond)), '-struct', 'out');
    end
end
end
