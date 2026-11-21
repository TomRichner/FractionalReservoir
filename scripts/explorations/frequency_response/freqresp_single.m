function [file, B] = freqresp_single(seed, opts)
%FREQRESP_SINGLE Single-neuron frequency responses: at operating points taken from the runs, and at matched operating points.
%
%   [file, B] = freqresp_single(seed, Name, Value, ...)
%
% THE RESPONSE. The gain from a neuron's dendritic input x_j to its
% synaptic output on one route, with its slow variables free, divided by
% the frozen gain (the limit of high frequency):
%   G(f) = d(a, j; 2 pi i f) / d_inf(a, j) ,  a = the excitatory target.
% |G| below 1 is attenuation by adaptation and depression; G -> 1 where
% they cannot follow.
%
% TWO ROUTES TO IT, COMPARED.
%   numeric  from the model's own Jacobian, for the one neuron: with s the
%            neuron's slow variables and i a postsynaptic neuron of type a
%            it has a synapse onto,
%            d = (tau_d / W(i, j)) [ J(x_i, x_j) + J(x_i, s) (z I - J(s, s))^(-1) J(s, x_j) ]
%   closed   freqresp_gain
% for adaptation and depression free ('both'), adaptation alone ('sfa')
% and depression alone ('std').
%
% OPERATING POINTS FROM THE RUNS (B.own), per adapting condition, at
% calibration state opts.state (default 7, t = 70.9 s) of seed [7 8]:
%   'E median'  reference network: the excitatory neuron of median rate
%               among those with nonzero slope
%   'I median'  the same for the inhibitory neurons
%   'E driven'  E-to-E x 1.5: the excitatory neuron at the 90th percentile
%               of rate among those with nonzero slope
%
% MATCHED OPERATING POINTS (B.matched), per adapting condition and rate
% r0 of opts.r0: every neuron at rate r0, slope 1, slow variables at
% steady state (freqresp_matched_state). Kept: the response of every
% excitatory neuron of the reference network (E-to-E route).
%
% B fields
%   f           1 x nf Hz (default 0.001 to 100, 501 points)
%   own         struct array: network, condition, label, neuron, slope,
%               rate, F, a (1 x K), s (route state), tau_a (1 x K), and
%               G.both / .sfa / .std (1 x nf complex, closed form),
%               Gnum (same, numeric), G0 (the three at f = 0),
%               max_diff (largest |closed - numeric| over |numeric|),
%               summary (freqresp_summary of each)
%   matched     struct array: condition, r0, F (n_E x 1, the frozen gain:
%               the route's output factor at this rate), and
%               G.both / .sfa / .std (n_E x nf complex), G0 (n_E x 1 each)
%
% See also: freqresp_gain, freqresp_summary, explore_frequency_response

arguments
    seed (1,2) double = [7 8]
    opts.f (1,:) double = logspace(-3, 2, 501)
    opts.state (1,1) double = 7
    opts.r0 (1,:) double = [0.15 0.25 0.5]
    opts.out_dir (1,:) char = ''
end

root = fileparts(which('setup_paths'));
out_dir = opts.out_dir;
if isempty(out_dir)
    out_dir = fullfile(root, 'data', 'explorations', 'frequency_response');
end
file = fullfile(out_dir, sprintf('single_s%d-%d.mat', seed(1), seed(2)));
if isfile(file)
    if nargout >= 2; L = load(file); B = L.B; end
    return;
end

f = opts.f;
z = [2i * pi * f, 0];
modes = struct('both', struct('sfa', true, 'std', true), ...
    'sfa', struct('sfa', true, 'std', false), 'std', struct('sfa', false, 'std', true));
names = fieldnames(modes);
conditions = {'sfa1_std1', 'sfa3_stdc2'};
picks = {'reference', 'E median', 1, 0.5; 'reference', 'I median', 2, 0.5; 'EEx1p5', 'E driven', 1, 0.9};

own = [];
for c = 1:numel(conditions)
    for p = 1:size(picks, 1)
        network = picks{p, 1};
        [~, params, ws] = freqresp_model(network, conditions{c}, seed);
        [~, Q] = freqresp_states(network, conditions{c}, seed, out_dir);
        S = Q.S_cal(:, opts.state);
        [~, d_inf, parts] = freqresp_gain(S, params, 0);
        idx = params.type_indices{picks{p, 3}};
        active = idx(parts.slope(idx) > 0);
        [~, order] = sort(parts.rate(active));
        j = active(order(max(1, round(picks{p, 4} * numel(active)))));

        o = struct('network', network, 'condition', conditions{c}, 'label', picks{p, 2}, ...
            'neuron', j, 'slope', parts.slope(j), 'rate', parts.rate(j), 'F', parts.F(1, j));
        q = picks{p, 3};
        o.a = S(params.state_layout.a{q}(ws.neuron_of_state(params.state_layout.a{q}) == j))';
        o.tau_a = params.tau_a_matrix{q}(idx == j, :);
        rows = params.route{q, 1}.idx;
        o.s = S(rows(ws.neuron_of_state(rows) == j))';

        J = SRNNCellTypePairs.compute_Jacobian_fast(S, params);
        o.max_diff = 0;
        for m = 1:numel(names)
            d = freqresp_gain(S, params, z, modes.(names{m}));
            G = reshape(d(1, j, :), 1, []) ./ d_inf(1, j);
            Gn = one_neuron_schur(J, params, ws, j, z, modes.(names{m})) ./ d_inf(1, j);
            o.G.(names{m}) = G(1:end - 1);
            o.Gnum.(names{m}) = Gn(1:end - 1);
            o.G0.(names{m}) = real(G(end));
            o.max_diff = max(o.max_diff, max(abs(G - Gn) ./ abs(Gn)));
            o.summary.(names{m}) = freqresp_summary(f, G(1:end - 1), real(G(end)));
        end
        own = [own, o]; %#ok<AGROW>
    end
end

matched = [];
for c = 1:numel(conditions)
    [~, params] = freqresp_model('reference', conditions{c}, seed);
    iE = params.type_indices{1};
    for r0 = opts.r0
        S = freqresp_matched_state(params, r0);
        mt = struct('condition', conditions{c}, 'r0', r0);
        [~, d_frozen] = freqresp_gain(S, params, 0);
        mt.F = d_frozen(1, iE)';                    % slope is 1: the frozen gain is the output factor
        for m = 1:numel(names)
            [d, d_inf] = freqresp_gain(S, params, z, modes.(names{m}));
            G = reshape(d(1, iE, :), numel(iE), []) ./ d_inf(1, iE)';
            mt.G.(names{m}) = G(:, 1:end - 1);
            mt.G0.(names{m}) = real(G(:, end));
        end
        matched = [matched, mt]; %#ok<AGROW>
    end
end

B = struct('seed', seed, 'f', f, 'state', opts.state, 'own', own, 'matched', matched);
if ~isfolder(out_dir); mkdir(out_dir); end
save(file, 'B');
end

function d = one_neuron_schur(J, params, ws, j, z, free)
% The gain of neuron j on its route to the excitatory type, from the
% Jacobian's blocks for that neuron alone.
L = params.state_layout;
mine = find(ws.neuron_of_state == j & ~ws.is_x)';
keep = false(size(mine));
if free.sfa
    keep = keep | ismember(mine, [L.a{:}]);
end
if free.std
    route_rows = [];
    for pre = 1:params.n_cellTypes
        for post = 1:params.n_cellTypes
            route_rows = [route_rows, params.route{pre, post}.idx]; %#ok<AGROW>
        end
    end
    keep = keep | ismember(mine, route_rows);
end
s = mine(keep);
post = params.type_indices{1};
post = post(post ~= j);
[~, best] = max(abs(ws.Wf(post, j)));
i = post(best);
xi = L.x(i); xj = L.x(j);
Jss = full(J(s, s)); Jxs = full(J(xi, s)); Jsx = full(J(s, xj));
d = complex(zeros(1, numel(z)));
for q = 1:numel(z)
    v = full(J(xi, xj));
    if ~isempty(s)
        v = v + Jxs * ((z(q) * eye(numel(s)) - Jss) \ Jsx);
    end
    d(q) = params.tau_d * v / ws.Wf(i, j);
end
end
