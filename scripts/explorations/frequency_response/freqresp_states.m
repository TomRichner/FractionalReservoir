function [file, Q] = freqresp_states(network, condition, seed, out_dir)
%FREQRESP_STATES The calibration states of an outlier run, simulated again and saved.
%
%   [file, Q] = freqresp_states(network, condition, seed, out_dir)
%
% WHY. outlier_run keeps the spectra at its 12 calibration states and not
% the states themselves, and the frequency response needs the Jacobian at
% those states. This simulates the same network on the same seeds for the
% same 80 s (serially, in the client, as outlier_run does) and takes the
% same rows: every 10 ms from 60 s, and of those the 12 that
% unique(round(linspace(1, n_s, 12))) picks.
%
% THE CHECK. When the outlier run is on disk, the rightmost eigenvalue of
% the frozen-gain matrix at each state is compared with the one saved
% there (R.dense.lam1_M, double precision; the condition without
% adaptation is compared through the rightmost eigenvalue of the full
% Jacobian, saved in single precision). Q.check_max is the largest
% absolute difference; the states are the run's own when it is at
% round-off.
%
% Q fields: network, condition, seed, t_cal (1 x 12), S_cal (N x 12),
% check_max (NaN when there is no outlier run to compare with), seconds.
%
% See also: outlier_run, freqresp_model, explore_frequency_response

arguments
    network   (1,:) char
    condition (1,:) char
    seed      (1,2) double = [7 8]
    out_dir   (1,:) char = ''
end
T_settle = 60; T_end = 80; dt_sample = 0.01; n_cal = 12;

root = fileparts(which('setup_paths'));
if isempty(out_dir)
    out_dir = fullfile(root, 'data', 'explorations', 'frequency_response');
end
if ~isfolder(out_dir); mkdir(out_dir); end
file = fullfile(out_dir, sprintf('states_%s_%s_s%d-%d.mat', network, condition, seed(1), seed(2)));
if isfile(file)
    L = load(file);
    Q = L.Q;
    return;
end

t0 = tic;
[model, params, ws] = freqresp_model(network, condition, seed, T_end);
model.run();
step = max(1, round(dt_sample * model.fs));
first = find(model.t_out >= T_settle, 1, 'first');
idx = first:step:numel(model.t_out);
cal_idx = unique(round(linspace(1, numel(idx), n_cal)));
t_cal = model.t_out(idx(cal_idx));
S_cal = model.S_out(idx(cal_idx), :)';
model = [];                                   %#ok<NASGU>  release the trajectory

check_max = NaN;
ref_file = fullfile(root, 'data', 'explorations', 'outliers', ...
    outlier_run_name(network, condition, seed, ''));
if isfile(ref_file)
    L = load(ref_file, 'R');
    lam_ref = [L.R.dense.lam1_M];
    lam_now = complex(zeros(1, numel(cal_idx)));
    for k = 1:numel(cal_idx)
        ev = eig(outlier_effective_matrix(S_cal(:, k), params, ws));
        [~, i1] = max(real(ev));
        lam_now(k) = ev(i1);
    end
    % A conjugate pair's rightmost member is either one.
    check_max = max(abs(complex(real(lam_now), abs(imag(lam_now))) ...
        - complex(real(lam_ref), abs(imag(lam_ref)))));
end

Q = struct('network', network, 'condition', condition, 'seed', seed, ...
    't_cal', t_cal(:)', 'S_cal', S_cal, 'check_max', check_max, 'seconds', toc(t0));
save(file, 'Q');
end
