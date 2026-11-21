function disk = frozen_gain_disk_at_states(model, idx)
% FROZEN_GAIN_DISK_AT_STATES The frozen gains at stored states of a run, for the uniform-disk reference.
%
%   disk = FROZEN_GAIN_DISK_AT_STATES(model, idx)
%
% model has run with store_full_state; idx are rows of its S_out. Shared by
% run_eig_heatmap (the baseline run; its disk_by_cond) and run_local_lyapunov
% (the switching run, FD-034), so the two stages save one definition.
%
% disk has the network size n, tau_d, the theoretical bulk radius R of W
% (the class's R), the times t of the states, and per state rms_gain: the
% root mean square over neurons of the frozen gain d_j (the class's
% frozen_gain: the slope of the nonlinearity times the route's output
% factor). The dendritic block of the Jacobian is (-I + W diag(d)) / tau_d,
% so a uniform disk of center -1/tau_d and radius R rms(d) / tau_d is what
% a random W would give (fig_eig_real_histogram draws its projection).
%
% See also: run_eig_heatmap, run_local_lyapunov, fig_eig_real_histogram,
%           jacobian_spectra_at_states
params = model.get_params();
t_used = model.t_out(idx);
rms_gain = zeros(numel(idx), 1);
for k = 1:numel(idx)
    d = feval([class(model) '.frozen_gain'], model.S_out(idx(k), :)', params);
    rms_gain(k) = sqrt(mean(d .^ 2));
end
disk = struct('n', model.n, 'tau_d', params.tau_d, 'R', model.R, 'rms_gain', rms_gain, 't', t_used(:));
end
