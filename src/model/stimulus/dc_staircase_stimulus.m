function [u_ex, t_ex] = dc_staircase_stimulus(params, T, fs, rng_seed, input_config)
% DC_STAIRCASE_STIMULUS Uniform tonic DC stepped through a sequence of levels,
% plus independent per-neuron white noise.
%
% Each level dc_levels(k) is held for hold_dur seconds, applied identically
% to every neuron over [0, T]. The first hold ramps in linearly from 0 over
% the first ramp_dur seconds. fs-invariant white noise
% (input_config.noise_intensity) is added per neuron on top. Signature
% matches the SRNNModel2 generator hook:
%   [u_ex, t_ex] = generator(params, T, fs, rng_seed, input_config)
%
% Required input_config fields:
%   dc_levels       - row vector of absolute DC levels (applied to all neurons)
%   hold_dur        - seconds each level is held
%   ramp_dur        - seconds to ramp 0 -> dc_levels(1) at the start
% Optional:
%   noise_intensity - fs-invariant white-noise intensity (0 or absent = off)
%
% Promoted from a local function in the bursting figure's own script so it lives
% on the path (added by setup_paths) and can be serialized to parfor workers by
% the multi-seed analyses.
    dt   = 1 / fs;
    t_ex = (0:dt:T)';          % nt x 1, matches built-in generator
    nt   = numel(t_ex);

    dc_levels = input_config.dc_levels;
    hold_dur  = input_config.hold_dur;
    ramp_dur  = input_config.ramp_dur;
    nL        = numel(dc_levels);

    % Staircase: level k over [(k-1)*hold_dur, k*hold_dur)
    dc_profile = zeros(nt, 1);
    for k = 1:nL
        seg = t_ex >= (k-1)*hold_dur & t_ex < k*hold_dur;
        dc_profile(seg) = dc_levels(k);
    end
    dc_profile(t_ex >= nL*hold_dur) = dc_levels(nL);   % final boundary sample

    % Ramp the first hold in linearly: 0 -> dc_levels(1) over the first ramp_dur s
    ramp_idx = t_ex <= ramp_dur;
    dc_profile(ramp_idx) = linspace(0, dc_levels(1), nnz(ramp_idx))';

    % Same drive to every neuron: n x nt
    u_ex = repmat(dc_profile', params.n, 1);

    % Add independent white noise per neuron over [0, T] (probes the network's
    % filtering). The sqrt(fs) factor makes this an fs-invariant white noise: the
    % continuous-time PSD ~ noise_intensity^2 is independent of fs (standard
    % Euler-Maruyama 1/sqrt(dt) scaling). The model's linear interpolant
    % band-limits it to ~Nyquist (fs/2), flat over our <100 Hz band. Seeded for
    % reproducibility.
    if isfield(input_config, 'noise_intensity') && input_config.noise_intensity > 0
        rng(rng_seed);
        u_ex = u_ex + input_config.noise_intensity * sqrt(fs) * randn(params.n, numel(t_ex));
    end
end
