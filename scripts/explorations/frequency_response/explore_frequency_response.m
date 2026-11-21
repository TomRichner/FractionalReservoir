% EXPLORE_FREQUENCY_RESPONSE How much do adaptation and depression attenuate a mode, as a function of how fast it rotates?
%
% THE QUESTION (Tom, 2026-11-20). "If there's a mode that's rotating
% quickly, adaptation can't necessarily keep up with it. ... how much does
% it attenuate modes with different oscillatory frequencies? ... could we
% make something of a frequency response for this? Or is this not well
% defined once you're at the network level?"
%
% THE ANSWER IN ONE PARAGRAPH. It is well defined, for the linearization
% about a state. Adaptation and depression variables are driven by their
% own neuron's rate only, so they can be eliminated neuron by neuron: a
% mode going as exp(z t) sees the weights with each presynaptic neuron's
% column scaled by a complex gain d_j(z) (freqresp_gain), and z is an
% eigenvalue of the full Jacobian exactly when 1 + tau_d z is an
% eigenvalue of M(z) = W .* d(z). The frozen-gain matrix of the outlier
% exploration is z -> infinity; z = 0 is the fully adapted network. The
% bulk of M(z) depends on the moduli of the gains only, so the bulk holds
% a growing mode rotating at omega when the LOOP GAIN
% R(i omega) / |1 + i omega tau_d| exceeds 1 (R the bulk radius from the
% gains). A structural outlier follows the complex mean gain, so for it
% the phase matters.
%
% WHAT IS RUN (seed [7 8], 500 neurons, the two networks and three
% conditions of the outlier exploration; serial, no pool):
%   freqresp_states    the 12 calibration states of each outlier run,
%                      simulated again (about 20 s a run) and checked
%                      against the spectra saved there
%   freqresp_network   per run: gains, predicted radius and mean-part
%                      eigenvalues over 0.001 to 100 Hz, all eigenvalues
%                      of M(i omega), the exactness check against the
%                      dense full Jacobian, the radius over the complex
%                      plane (1 to 5 minutes a run)
%   freqresp_single    single-neuron responses, numeric against closed
%                      form, at operating points from the runs and at
%                      matched ones
%   fig_frequency_response   the four figures and their sidecars
% Every step keeps its result on disk and is skipped when it is there
% (delete a file under data/explorations/frequency_response/ to redo it).
% Then the checks and the numbers behind the sidecars are printed.
%
% See also: freqresp_gain, freqresp_network, freqresp_single, fig_frequency_response

setup_paths();
addpath(fileparts(mfilename('fullpath')));
if ~exist('seed', 'var'); seed = [7 8]; end
networks = {'reference', 'EEx1p5'};
conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
data_dir = fullfile(fileparts(which('setup_paths')), 'data', 'explorations', 'frequency_response');

%% Runs
for w = 1:numel(networks)
    for c = 1:numel(conditions)
        [~, Q] = freqresp_states(networks{w}, conditions{c}, seed);
        fprintf('[states] %-9s %-13s rightmost frozen eigenvalue against the outlier run: max difference %.2g\n', ...
            networks{w}, conditions{c}, Q.check_max);
        freqresp_network(networks{w}, conditions{c}, seed);
    end
end
[~, B] = freqresp_single(seed);
fig_frequency_response('seed', seed);

%% Checks: is the elimination exact and local?
fprintf('\nEXACTNESS (z_k eigenvalues of the dense full Jacobian; sigma = smallest singular value of M(z_k) - (1 + tau_d z_k) I over norm(M(z_k)))\n');
for w = 1:numel(networks)
    for c = 2:numel(conditions)
        L = load(fullfile(data_dir, sprintf('network_%s_%s_s%d-%d.mat', networks{w}, conditions{c}, seed(1), seed(2))));
        A = L.A;
        for E = A.exact
            off = E.dslow > 1e-6;                    % not an eigenvalue of the slow block
            lead = E.kind <= 2;
            fprintf(['  %-9s %-11s state %2d: %d eigenvalues tested; %d also eigenvalues of the slow block (reduced problem undefined there); ' ...
                'of the other %d: sigma median %.1e, max %.1e (leading and continuation: max %.1e); ' ...
                'closed form against Schur %.1e; off-diagonal remainder %.1e; routes of one neuron differ by %.1e\n'], ...
                networks{w}, conditions{c}, E.state, numel(E.z), nnz(~off), nnz(off), ...
                median(E.sigma(off)), max(E.sigma(off)), max(E.sigma(off & lead)), ...
                max(E.closed_diff), max(E.offdiag), A.pop.route_diff);
        end
    end
end
fprintf('  single neuron, closed form against the Jacobian''s blocks: largest relative difference %.1e\n', max([B.own.max_diff]));

%% Checks: the bulk does not see the phase of the gains
fprintf('\nBULK AND PHASE (largest eigenvalue modulus of W .* d against W .* |d| at 0.7 Hz, state 7; predicted radius)\n');
for w = 1:numel(networks)
    for c = 2:numel(conditions)
        [~, params, ws] = freqresp_model(networks{w}, conditions{c}, seed);
        [~, Q] = freqresp_states(networks{w}, conditions{c}, seed);
        d = freqresp_gain(Q.S_cal(:, 7), params, 2i * pi * 0.7);
        e1 = eig(freqresp_matrix(d, ws)); e2 = eig(freqresp_matrix(abs(d), ws));
        fprintf('  %-9s %-11s complex gains: 95th percentile modulus %.3f, max %.3f; moduli only: %.3f, %.3f; predicted R %.3f\n', ...
            networks{w}, conditions{c}, prctile(abs(e1), 95), max(abs(e1)), prctile(abs(e2), 95), max(abs(e2)), freqresp_radius(d, ws));
    end
end

%% One against multiple timescales: what the difference is made of
fprintf('\nONE AGAINST MULTIPLE TIMESCALES, loop gain R(f) / |1 + 2 pi i f tau_d| of the bulk (median over 12 states)\n');
for w = 1:numel(networks)
    L1 = load(fullfile(data_dir, sprintf('network_%s_sfa1_std1_s%d-%d.mat', networks{w}, seed(1), seed(2)))); A1 = L1.A;
    L3 = load(fullfile(data_dir, sprintf('network_%s_sfa3_stdc2_s%d-%d.mat', networks{w}, seed(1), seed(2)))); A3 = L3.A;
    den = abs(1 + 2i * pi * A1.f * A1.tau_d);
    n1 = median(A1.R.both, 1) / median(A1.R_inf); n3 = median(A3.R.both, 1) / median(A3.R_inf);
    pk = @(v) max(v ./ den);
    fprintf(['  %-9s frozen radius %.3f (one) %.3f (multiple). Peak loop gain: one %.3f, multiple %.3f; ' ...
        'multiple''s filter on one''s frozen radius %.3f; one''s filter on multiple''s frozen radius %.3f\n'], ...
        networks{w}, median(A1.R_inf), median(A3.R_inf), pk(n1 * median(A1.R_inf)), pk(n3 * median(A3.R_inf)), ...
        pk(n3 * median(A1.R_inf)), pk(n1 * median(A3.R_inf)));
    fprintf('            frozen radius is made of: fraction with nonzero slope %.2f / %.2f, mean slope %.2f / %.2f, mean output factor %.3f / %.3f, rms frozen gain %.3f / %.3f\n', ...
        mean(A1.pop.slope(:) > 0), mean(A3.pop.slope(:) > 0), mean(A1.pop.slope(:)), mean(A3.pop.slope(:)), ...
        mean(A1.pop.d_inf(:)) / mean(A1.pop.slope(:)), mean(A3.pop.d_inf(:)) / mean(A3.pop.slope(:)), ...
        sqrt(mean(A1.pop.d_inf(:).^2)), sqrt(mean(A3.pop.d_inf(:).^2)));
end
m = B.matched;
for r0 = unique([m.r0])
    a1 = m(strcmp({m.condition}, 'sfa1_std1') & [m.r0] == r0);
    a3 = m(strcmp({m.condition}, 'sfa3_stdc2') & [m.r0] == r0);
    ratio = abs(a1.G.sfa) ./ abs(a3.G.sfa);
    fprintf(['  matched, rate %.2f: adaptation alone, gain one / multiple: largest over neurons and frequencies %.6f, smallest %.3f; ' ...
        'output factor (median) %.3f / %.3f; G(0) both %.3f / %.3f\n'], r0, max(ratio(:)), min(ratio(:)), ...
        median(a1.F), median(a3.F), median(a1.G0.both), median(a3.G0.both));
end
