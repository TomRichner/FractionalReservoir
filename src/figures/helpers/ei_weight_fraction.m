function v = ei_weight_fraction(psa, res, verbose)
% EI_WEIGHT_FRACTION The realised E:I weight balance of the network a param-space result ran on.
%
%   v = EI_WEIGHT_FRACTION(psa, res)
%   v = EI_WEIGHT_FRACTION(psa, res, verbose)
%
%   v = |sum W(:, E)| / (|sum W(:, E)| + |sum W(:, I)|)
%
% summed over presynaptic columns: total excitatory drive against total
% inhibitory drive, in (0, 1), read as an E:I weight ratio. The network is
% REBUILT (psa.rebuild_model(res), then build()), since a stored result
% holds scalars, never W.
%
% WHY ONE HELPER WITH A MEMORY. fig_EI_weights_param_space and
% fig_lle_vs_rate each carried a copy of this function with a cache per
% call, keyed by config_idx. A rebuild takes up to about 2 s at n = 1000,
% and the joint sample of the fast run holds 128 networks, so main Figure 4,
% which draws both, paid for every network twice, and every redraw paid
% again (about 25 minutes for one call of fig_EI_weights_param_space on
% 2026-11-22, with the machine loaded). Here the values are remembered for
% the MATLAB session, keyed by what fixes W: the model class, the sweep's
% model defaults, the vector-parameter lookup, the grid point's
% configuration and its network seed. Not the condition, which sets the
% adaptation and not W (the per-call caches it replaces relied on the same).
% The key is a digest of those values, not a path, so a copy of a run
% directory (test_canvas_standard draws from one) finds the values its
% original left. `clear ei_weight_fraction` forgets them.
%
% See also: fig_EI_weights_param_space, fig_lle_vs_rate,
%           ParamSpaceAnalysis2/rebuild_model

if nargin < 3; verbose = 'near-none'; end
persistent memo
if isempty(memo); memo = containers.Map('KeyType', 'char', 'ValueType', 'double'); end
cfg_point = struct();
if isfield(res, 'config'); cfg_point = res.config; end
seed = [];   % none: rebuild_model below says why the network is out of reach
if isfield(res, 'network_seed'); seed = res.network_seed; end
key = digest({psa.model_class, psa.model_defaults, psa.vector_param_lookup, cfg_point, seed});
if isKey(memo, key)
    v = memo(key);
    return
end
m = psa.rebuild_model(res);
m.verbose = verbose;   % the model prints its build report only at 'verbose'
m.build();
ti = m.type_indices;
assert(numel(ti) >= 2, 'ei_weight_fraction:CellTypes', ...
    'ei_weight_fraction needs at least two cell types; got %d.', numel(ti));
S_E = full(sum(sum(m.W(:, ti{1}))));
S_I = full(sum(sum(m.W(:, ti{2}))));
denom = abs(S_E) + abs(S_I);
if denom == 0
    v = 0.5;   % a zero W has no balance; centre it rather than divide by zero
else
    v = abs(S_E) / denom;
end
memo(key) = v;
end

function key = digest(x)
md = java.security.MessageDigest.getInstance('SHA-1');
md.update(getByteStreamFromArray(x));
key = sprintf('%02x', typecast(md.digest(), 'uint8'));
end
