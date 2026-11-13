function [u_ex, t_ex] = pulse_train_input(params, T, fs, rng_seed, config) %#ok<INUSD>
% PULSE_TRAIN_INPUT A DC step or a rectangular pulse train of a given mean, the same on every neuron.
%
%   [u_ex, t_ex] = PULSE_TRAIN_INPUT(params, T, fs, rng_seed, config)
%
% An input generator for SRNNCellTypePairs (input_config.generator; see
% SRNNCellTypePairs.generate_external_input). Called on t = 0:1/fs:T; the
% model prepends zeros for t < 0 itself. Draws no random numbers, so
% rng_seed is unused.
%
% config fields (beside `generator`)
%   pulse_mean   time-averaged input while the stimulus is on (input units)
%   pulse_freq   pulse frequency in Hz; 0 = DC (the mean, constant)
%   pulse_duty   fraction of each period the pulse is on (0 < duty <= 1)
%   t_on, t_off  the stimulus is on for t_on <= t < t_off (s)
%
% Pulses are rectangular, positive only, in phase on every neuron, with
% amplitude pulse_mean / pulse_duty, the first starting at t_on. They are
% laid on the sample grid: the period must be a whole number of samples
% (fs / pulse_freq), the width too (duty * period), and the on time a whole
% number of periods, so the mean of the samples over [t_on, t_off) is
% pulse_mean exactly; the function errors otherwise rather than round.
% The model interpolates u linearly between samples, and the integral of
% that interpolant over a whole number of periods equals the sample mean,
% so the input the integrator sees has the same mean too.
%
% See also: p5_dc_vs_pulses, SRNNCellTypePairs.generate_external_input

t_ex = (0:1/fs:T)';
k = round(t_ex * fs);                       % sample index, exact integers
k_on  = round(config.t_on  * fs);
k_off = round(config.t_off * fs);
if abs(k_on - config.t_on * fs) > 1e-9 || abs(k_off - config.t_off * fs) > 1e-9
    error('pulse_train_input:Grid', 't_on and t_off must fall on the sample grid.');
end
on_window = k >= k_on & k < k_off;

if config.pulse_freq == 0
    u = config.pulse_mean * double(on_window);
else
    P  = fs / config.pulse_freq;            % period, samples
    Wd = config.pulse_duty * P;             % width, samples
    if abs(P - round(P)) > 1e-9 || abs(Wd - round(Wd)) > 1e-9 || round(Wd) < 1
        error('pulse_train_input:Grid', ...
            'Period %.4g and width %.4g samples must be whole numbers (fs %g Hz).', P, Wd, fs);
    end
    P = round(P); Wd = round(Wd);
    if mod(k_off - k_on, P) ~= 0
        error('pulse_train_input:Grid', 'The on time must be a whole number of periods.');
    end
    amp = config.pulse_mean / config.pulse_duty;
    u = amp * double(on_window & mod(k - k_on, P) < Wd);
end
u_ex = repmat(u(:)', params.n, 1);
end
