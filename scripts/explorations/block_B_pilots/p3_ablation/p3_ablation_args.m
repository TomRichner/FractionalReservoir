function spec = p3_ablation_args(ablation, ref)
% P3_ABLATION_ARGS The model arguments that turn an intact condition into one of pilot 3's ablations.
%
%   spec = P3_ABLATION_ARGS(ablation, ref)
%
% Everything an ablation changes is expressed as name-value model arguments
% for the condition builders' 'model_args' (applied last, so they override
% the condition), plus, for the mean-matched ablations, a constant
% per-neuron input. No model code is involved.
%
% ablation   one of P3_ABLATION_ARGS('names')
% ref        the INTACT run of the same network and condition, as
%            p3_run_one returns it (its .ref field): the condition's own
%            tau_a, synapse_config and c, and the intact run's time
%            averages. Not needed for 'intact'.
%
% spec.model_args   name-value cell for the builder's 'model_args'
% spec.drive        [] or n x 1: a constant input u_i for the whole run
%                   (settling included; p3_run_one delivers it)
% spec.col_scale    [] or n x 1: per-presynaptic-neuron weight factor,
%                   applied after the build with replace_connectivity
%                   (IEDExplorationNetwork only, so run13 only)
% spec.label        short column title
% spec.how          one sentence saying what was done
%
% THE ABLATIONS
%   intact             nothing changed.
%   sfa_removed        tau_a empty on every type; depression intact.
%   std_removed        the std field removed from every route; SFA intact.
%   both_removed       both of the above. (Then 'original' and 'chain_mts'
%                      are the same network.)
%   mean_matched       both removed, and each adapting variable replaced by
%                      a constant equal to its time average in the intact
%                      run:
%                      depression: each depressed route gets the weight
%                        scale <D>, the intact run's depression factor
%                        (std_route_factor) averaged over the route's
%                        presynaptic neurons and over time. One number per
%                        route, because a route's scale is a scalar.
%                      SFA: neuron i gets the constant input
%                        u_i = -<(c/K) sum_k a_ik>, its own time average.
%                        The rate is phi(x - (c/K) sum a - S_c); with
%                        y = x - abar the equation for y is the original
%                        with the adaptation term frozen at abar and the
%                        input u - abar, so a constant input of -abar is
%                        EXACTLY a per-neuron rise of the setpoint by abar
%                        (p3_check_offset tests this numerically).
%   mean_matched_rw    as mean_matched, but the route scale is the
%                      rate-weighted average <r D> / <r>: the constant that
%                      transmits the same time-averaged synaptic output
%                      r D as the intact run did. (During an event D is low
%                      while r is high, so this is below the plain <D>.)
%   mean_matched_pn    as mean_matched, but depression is matched per
%                      presynaptic neuron: column j of W times <D_j>.
%                      Needs replace_connectivity, so run13 only.
%   slow_removed       the multiple timescales removed, strength kept: each
%                      type keeps only its fastest SFA rung (c unchanged:
%                      c is the total budget, so the one rung carries all
%                      of it), each depressed route keeps only its first
%                      pool (a chain's releasable pool with the reserve
%                      held full; a product's first pool).
%   slow_removed_share as slow_removed, but the fastest SFA rung keeps only
%                      the share it had (c divided by the number of rungs),
%                      i.e. the slow rungs are deleted, not folded in.
%
% See also: p3_run_one, run_p3_ablation

names = {'intact', 'sfa_removed', 'std_removed', 'both_removed', 'mean_matched', ...
    'mean_matched_rw', 'mean_matched_pn', 'slow_removed', 'slow_removed_share'};
if nargin == 1 && strcmp(ablation, 'names')
    spec = names;
    return;
end
assert(ismember(ablation, names), 'Unknown ablation "%s".', ablation);
spec = struct('model_args', {{}}, 'drive', [], 'col_scale', [], 'label', '', 'how', '');
if strcmp(ablation, 'intact')
    spec.label = 'intact';
    spec.how = 'Nothing changed.';
    return;
end

C = numel(ref.names);
no_sfa = repmat({zeros(1, 0)}, 1, C);
switch ablation
    case 'sfa_removed'
        spec.model_args = {'tau_a', no_sfa};
        spec.label = 'SFA removed';
        spec.how = 'tau_a empty on every type; depression intact.';
    case 'std_removed'
        spec.model_args = {'synapse_config', strip_std(ref.synapse_config, ref.names, [])};
        spec.label = 'depression removed';
        spec.how = 'std removed from every route; SFA intact.';
    case 'both_removed'
        spec.model_args = {'tau_a', no_sfa, ...
            'synapse_config', strip_std(ref.synapse_config, ref.names, [])};
        spec.label = 'both removed';
        spec.how = 'tau_a empty and std removed.';
    case 'mean_matched'
        spec.model_args = {'tau_a', no_sfa, ...
            'synapse_config', strip_std(ref.synapse_config, ref.names, ref.D_route)};
        spec.drive = -ref.abar;
        spec.label = 'mean-matched';
        spec.how = ['Both removed; route scale = time- and population-averaged depression factor, ' ...
            'constant input = minus each neuron''s time-averaged adaptation term.'];
    case 'mean_matched_rw'
        spec.model_args = {'tau_a', no_sfa, ...
            'synapse_config', strip_std(ref.synapse_config, ref.names, ref.D_route_rw)};
        spec.drive = -ref.abar;
        spec.label = 'mean-matched, rate-weighted';
        spec.how = ['Both removed; route scale = <r D> / <r> (same mean synaptic output), ' ...
            'constant input = minus each neuron''s time-averaged adaptation term.'];
    case 'mean_matched_pn'
        spec.model_args = {'tau_a', no_sfa, ...
            'synapse_config', strip_std(ref.synapse_config, ref.names, [])};
        spec.drive = -ref.abar;
        spec.col_scale = neuron_scale(ref);
        spec.label = 'mean-matched, per neuron';
        spec.how = ['Both removed; column j of W times neuron j''s time-averaged depression factor, ' ...
            'constant input = minus each neuron''s time-averaged adaptation term.'];
    case {'slow_removed', 'slow_removed_share'}
        tau_a = ref.tau_a;
        for q = 1:C
            if ~isempty(tau_a{q}), tau_a{q} = tau_a{q}(1); end
        end
        spec.model_args = {'tau_a', tau_a, ...
            'synapse_config', first_pool(ref.synapse_config)};
        spec.label = 'slow parts removed';
        spec.how = 'Fastest SFA rung only (c unchanged, the total budget); first depression pool only.';
        if strcmp(ablation, 'slow_removed_share')
            spec.model_args = [spec.model_args, {'c', ref.c ./ max(1, ref.n_a)}];
            spec.label = 'slow parts removed, c/K';
            spec.how = 'Fastest SFA rung only, with the share of c it had (c/K); first depression pool only.';
        end
end
end

function sc = strip_std(sc, names, D_scale)
% Remove std from every route. With D_scale (C x C, pre by post, NaN where
% the route was not depressed) the route's scale is multiplied by it.
for i = 1:numel(names)
    if ~isfield(sc, names{i}), continue; end
    for j = 1:numel(names)
        if ~isfield(sc.(names{i}), names{j}), continue; end
        route = sc.(names{i}).(names{j});
        if ~isstruct(route) || ~isfield(route, 'std'), continue; end
        route = rmfield(route, 'std');
        if ~isempty(D_scale)
            s = 1;
            if isfield(route, 'scale') && ~isempty(route.scale), s = route.scale; end
            assert(isfinite(D_scale(i, j)) && D_scale(i, j) > 0, 'No depression average for route %s->%s.', names{i}, names{j});
            route.scale = s * D_scale(i, j);
        end
        sc.(names{i}).(names{j}) = route;
    end
end
end

function sc = first_pool(sc)
% Keep only the first pool of every depressed route, as a one-pool product
% route (for a chain: the releasable pool with the reserve held full).
pres = fieldnames(sc);
for i = 1:numel(pres)
    posts = fieldnames(sc.(pres{i}));
    for j = 1:numel(posts)
        route = sc.(pres{i}).(posts{j});
        if ~isstruct(route) || ~isfield(route, 'std'), continue; end
        route.std = struct('tau_rec', route.std.tau_rec(1), 'tau_rel', route.std.tau_rel(1));
        sc.(pres{i}).(posts{j}) = route;
    end
end
end

function d = neuron_scale(ref)
% One factor per presynaptic neuron. Only defined when all of a type's
% routes are depressed alike (then the factor is the same on each) or none.
C = numel(ref.names);
d = ones(ref.n, 1);
for pre = 1:C
    has = ~cellfun(@isempty, ref.D_neuron(pre, :));
    if ~any(has), continue; end
    assert(all(has), 'Type %s has depressed and undepressed routes; a per-neuron factor is not defined.', ref.names{pre});
    first = ref.D_neuron{pre, 1};
    for post = 2:C
        assert(max(abs(ref.D_neuron{pre, post} - first)) < 1e-12, ...
            'The routes of type %s differ; a per-neuron factor is not defined.', ref.names{pre});
    end
    d(ref.type_indices{pre}) = first;
end
end
