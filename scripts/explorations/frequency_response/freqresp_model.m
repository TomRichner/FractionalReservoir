function [model, params, ws] = freqresp_model(network, condition, seed, T_end)
%FREQRESP_MODEL The network of an outlier run, built and not simulated.
%
%   [model, params, ws] = freqresp_model(network, condition, seed, T_end)
%
% The same preset, network change and seeds as outlier_run (the preset
% name and the network rules are repeated here because outlier_run holds
% them inline; if one changes, change the other). The Lyapunov pass is
% off: it does not touch the trajectory, and freqresp_states checks the
% states it gets against the spectra outlier_run saved.
%
%   network    'reference' or 'EEx<a>p<b>' (E-to-E mean times a.b)
%   condition  'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'
%   seed       1 x 2, default [7 8]
%   T_end      length of the run in seconds, default 80 (60 to settle and
%              20 measured; a run of another length is another noise
%              realization)
%
%   params = model.get_params(), ws = outlier_weight_stats(params)
%
% See also: outlier_run, freqresp_states

arguments
    network   (1,:) char
    condition (1,:) char {mustBeMember(condition, {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'})}
    seed      (1,2) double = [7 8]
    T_end     (1,1) double = 80
end

preset = ['celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_' ...
    'stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain'];
tok = regexp(network, '^EEx(\d+)p(\d+)$', 'tokens', 'once');
if strcmp(network, 'reference')
    net_args = {};
elseif ~isempty(tok)
    net_args = {'mu_EE_relative', str2double([tok{1} '.' tok{2}]) * 7};
else
    error('freqresp_model:Network', 'Unknown network ''%s'' (reference, EEx<a>p<b>).', network);
end
model = build_from_preset(preset, condition, net_args{:}, ...
    'rng_seeds', seed, 'T_range', [0 T_end], 'fs', 500, ...
    'store_full_state', true, 'verbose', 'near-none', 'lya_method', 'none');
if nargout >= 2
    params = model.get_params();
    if any(params.release_rate ~= 0)
        error('freqresp_model:Release', 'The closed forms of freqresp_gain assume no spontaneous release.');
    end
end
if nargout >= 3
    ws = outlier_weight_stats(params);
end
end
