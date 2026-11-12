function model = route_pool()
%ROUTE_POOL Route model: vesicle pools in series (chained depression).
%
%   model = ROUTE_POOL()
%
% The vesicle-pool family behind the interface of route_model (read that
% header first). Registered as 'pool'; a route gets it by asking for
% synapse_config.<pre>.<post>.std.model = 'chain'.
%
% THE EQUATIONS, per presynaptic neuron with rate r. The state is
% s = [b1, b2]: b1 the RELEASABLE pool, b2 the RESERVE that refills it,
% both as fractions of their own capacity, both 1 at rest.
%
%   db1/dt = b2 (1 - b1)/tau_rec1 - b1 r/tau_rel                      (C1)
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2) b2 (1 - b1)/tau_rec1          (C2)
%   F      = b1                                                       (C3)
%
% Term by term: the releasable pool empties by release, b1 r/tau_rel, and
% is refilled from the reserve at a rate proportional to what the reserve
% holds (b2) and to the releasable pool's empty sites (1 - b1). That refill
% flux is what the reserve loses, divided by n2 because the reserve's
% capacity is n2 times the releasable pool's (the same vesicles are a
% smaller fraction of a larger pool). The reserve itself refills toward 1
% with tau_rec2. Only the releasable pool releases, so only b1 is the
% output factor; the reserve acts on the output through b1 alone.
%
% With the reserve held full (b2 = 1) (C1) is the single pool of
% route_product, db/dt = (1 - b)/tau_rec - b r/tau_rel; n2 -> infinity is
% that limit. A route that wants the single pool asks route_product for one
% component and carries no reserve variable.
%
% RANGE. b1 and b2 stay in [0, 1]: at b1 = 0 its derivative is
% b2/tau_rec1 >= 0 and at b1 = 1 it is -r/tau_rel <= 0; at b2 = 0 its
% derivative is 1/tau_rec2 > 0 and at b2 = 1 it is -(1 - b1)/(n2 tau_rec1)
% <= 0.
%
% THE LINEARISATION, with u = 1/tau_rec1, v = 1/tau_rec2, w = 1/tau_rel,
% q = 1/n2 (A(i, m, k) = d(ds_m)/d(s_k), B = d(ds)/dr, Fs = dF/ds):
%
%              d/db1            d/db2                  d/dr
%   db1/dt   -b2 u - r w      (1 - b1) u             -b1 w
%   db2/dt    q b2 u          -v - q (1 - b1) u        0
%   F          1                0
%
% The 2 x 2 is full and NOT symmetric, so P.coupling is true(2).
%
% THE STEADY STATE at a constant rate r. With k1 = r tau_rec1/tau_rel,
% k2 = r tau_rec2/(n2 tau_rel) and m = 1 + k1 - k2, (C2) = 0 gives
% b2 = 1 - k2 b1 and (C1) = 0 then k2 b1^2 - (1 + k1 + k2) b1 + 1 = 0, whose
% root in [0, 1] is
%
%   b1 = 2 / (1 + k1 + k2 + d),   d = sqrt(m^2 + 4 k1 k2)             (C4)
%   b2 = 1 - k2 b1
%
% (d^2 is (1 + k1 + k2)^2 - 4 k2 written as a sum of two nonnegative
% terms.) The quadratic formula's other form, (s - d)/(2 k2), cancels at
% small k2; this one does not, and it reduces to the single pool's
% 1/(1 + k1) at k2 = 0. b2 is evaluated as (m + d)/(1 + k1 + k2 + d), or as
% 4 k1 k2/((d - m)(1 + k1 + k2 + d)) where m < 0, which is the same number
% without the cancellation of 1 - k2 b1 when the reserve is nearly empty.
% The sustained output r b1 rises with the rate to the ceiling
% 1/(tau_rec1/tau_rel + tau_rec2/(n2 tau_rel)) (the product's falls).
%
% CONFIGURATION (synapse_config.<pre>.<post>.std):
%   model     'chain'
%   tau_rec   1 x 2 positive: [releasable-pool refill, reserve refill]
%   tau_rel   positive scalar: ONE release time, the releasable pool's
%   n2        positive finite scalar: reserve capacity over releasable-pool
%             capacity. No default.
% Refused, each with its own identifier (<prefix> is opts.err_id_prefix):
%   <prefix>:ChainTauRecCount     tau_rec has not exactly two values
%   <prefix>:ChainTauRelNotScalar tau_rel is not a scalar
%   <prefix>:ChainMissingN2       n2 is absent
%   <prefix>:ChainInvalidN2       n2 is not a positive finite scalar
%   <prefix>:ChainCombine         std.combine (the product's option)
%   <prefix>:ChainZeroFloor       opts.zero_floor (the class's std_zero_floor)
%   <prefix>:ChainStfNotBuilt     facilitation on the route (see below)
%   <prefix>:InvalidSynapseConfig anything else malformed
%
% P, beyond the fields every model has (route_model), PER NEURON, npre x 1:
%   tau_rec1, tau_rec2, tau_rel, n2
% and nominal: tau_rec (1 x 2), tau_rel, n2 as configured, std_model
% ('chain'), and tau_dec, tau_fac, G empty.
%
% ONE HANDLE BEYOND THE INTERFACE:
%   D = depression(b, P)   the depression factor from the route's b-block
%                          alone (npre x n_b, or npre x n_b x nt; D has
%                          the shape of prod(b, 2)): b1. It is what
%                          SRNNCellTypePairs.std_route_factor returns for
%                          display code that holds a b array and no g.
%
% WHAT COMES NEXT (FD-018; Plan A, section 4.6): release-coupled
% facilitation extends this file. A release-probability variable g
% multiplies the release term of (C1), b1 r/tau_rel -> g b1 r/tau_rel, and
% the output, F = g b1; optionally a calcium trace h drives g in place of
% the rate. They are g-block columns (n_g = 1 or 2), after the b-block in
% s, so the columns of b1 and b2 do not move; A gains the entry
% d(db1)/dg and the rows of g and h, and B and Fs their entries. A route
% may then also carry one pool or none, with the variable removed and not
% clamped. None of that is here: today the model has exactly the two
% pools, and facilitation on a chain route is refused.
%
% See also: route_model, route_product, SRNNCellTypePairs,
%           test_route_models, test_std_chain_model

model = struct('name', 'pool', ...
    'compile', @compile, 'initial', @initial, 'rhs', @rhs, 'output', @output, ...
    'linearize', @linearize, 'steady_state', @steady_state, ...
    'depression', @depression);
end

%% The interface

function P = compile(std_cfg, stf_cfg, opts)
if nargin < 3, opts = struct(); end
opts = fill_defaults(opts, struct('npre', 1, 'zero_floor', false, 'label', 'route', ...
    'err_id_prefix', 'route_pool'));
prefix = opts.err_id_prefix;
label = opts.label;
err_id = [prefix ':InvalidSynapseConfig'];

if ~is_absent(stf_cfg)
    error([prefix ':ChainStfNotBuilt'], ...
        ['Route %s has chained depression (std.model = ''chain'') and facilitation (stf). ' ...
        'Facilitation on a chain route is not built yet: the release-coupled form is ' ...
        'FD-018, and the output-only form is not offered on a chain route.'], label);
end
if opts.zero_floor
    error([prefix ':ChainZeroFloor'], ...
        ['std_zero_floor cannot be used with chained depression (route %s has ' ...
        'std.model = ''chain''). The floor maps a product route''s factor at r = 1 to ' ...
        'zero; the chain''s factor is the releasable pool itself and has no such map.'], label);
end
if ~isstruct(std_cfg) || ~isscalar(std_cfg)
    error(err_id, 'STD configuration on route %s must be a scalar struct.', label);
end
fields = fieldnames(std_cfg);
if any(strcmp('combine', fields))
    error([prefix ':ChainCombine'], ...
        ['std.combine on route %s: combine belongs to product depression and has no ' ...
        'meaning for std.model = ''chain'', whose factor is the releasable pool b1.'], label);
end
if ~any(strcmp('n2', fields))
    error([prefix ':ChainMissingN2'], ...
        ['Chained depression on route %s needs std.n2, the reserve''s capacity relative ' ...
        'to the releasable pool''s. It has no default.'], label);
end
for k = {'tau_rec', 'tau_rel'}
    if ~any(strcmp(k{1}, fields))
        error(err_id, 'Missing STD field on route %s: %s.', label, k{1});
    end
end
unknown = setdiff(fields, {'model', 'tau_rec', 'tau_rel', 'n2'});
if ~isempty(unknown)
    error(err_id, 'Unknown STD field on chain route %s: %s.', label, unknown{1});
end
if isfield(std_cfg, 'model') && ~strcmp(char(string(std_cfg.model)), 'chain')
    error(err_id, 'route_pool compiles std.model = ''chain''; route %s states ''%s''.', ...
        label, char(string(std_cfg.model)));
end

tau_rec = std_cfg.tau_rec;
tau_rel = std_cfg.tau_rel;
n2 = std_cfg.n2;
if ~is_positive(tau_rec)
    error(err_id, 'tau_rec on route %s must contain positive finite values.', label);
end
if numel(tau_rec) ~= 2
    error([prefix ':ChainTauRecCount'], ...
        ['Chained depression on route %s needs exactly two tau_rec values, ' ...
        '[releasable-pool refill, reserve refill]; got %d.'], label, numel(tau_rec));
end
if ~is_positive(tau_rel)
    error(err_id, 'tau_rel on route %s must contain positive finite values.', label);
end
if ~isscalar(tau_rel)
    error([prefix ':ChainTauRelNotScalar'], ...
        ['Chained depression on route %s needs a scalar tau_rel: only the releasable ' ...
        'pool releases. Got %d values.'], label, numel(tau_rel));
end
if ~(is_positive(n2) && isscalar(n2))
    error([prefix ':ChainInvalidN2'], ...
        'n2 on route %s must be a positive finite scalar.', label);
end

nominal = struct('tau_rec', reshape(double(tau_rec), 1, 2), 'tau_rel', double(tau_rel), ...
    'tau_dec', [], 'tau_fac', [], 'G', [], 'n2', double(n2), 'std_model', 'chain');

P = struct('model', 'pool');
P.n_b = 2;                    % b1 (column 1) and b2 (column 2)
P.n_g = 0;
P.ns = 2;
P.npre = opts.npre;
P.nominal = nominal;
P.signature = nominal;
P.signature.model = P.model;
P.coupling = true(2);         % each pool's derivative depends on both pools

% Per neuron, the nominal value replicated.
P.tau_rec1 = repmat(nominal.tau_rec(1), opts.npre, 1);
P.tau_rec2 = repmat(nominal.tau_rec(2), opts.npre, 1);
P.tau_rel  = repmat(nominal.tau_rel, opts.npre, 1);
P.n2       = repmat(nominal.n2, opts.npre, 1);
end

function s0 = initial(P) %#ok<INUSD>
% Both pools start full, their value at rest.
s0 = [1 1];
end

function [ds, F] = rhs(s, r, P)
b1 = s(:, 1);
b2 = s(:, 2);
refill = b2 .* (1 - b1) ./ P.tau_rec1;        % reserve -> releasable pool, in units of the releasable pool
ds = [refill - b1 .* r ./ P.tau_rel, (1 - b2) ./ P.tau_rec2 - refill ./ P.n2];
F = b1;
end

function F = output(s, P) %#ok<INUSD>
F = s(:, 1, :);
end

function [A, B, F, Fs] = linearize(s, r, P)
npre = size(s, 1);
b1 = s(:, 1);
b2 = s(:, 2);
u = 1 ./ P.tau_rec1;
q = 1 ./ P.n2;
empty_u = (1 - b1) .* u;                      % d(refill)/d(b2)
full_u = b2 .* u;                             % -d(refill)/d(b1)
% A(i, m, k) = d(ds_m)/d(s_k); the columns below are (1,1), (2,1), (1,2), (2,2).
A = reshape([-full_u - r ./ P.tau_rel, q .* full_u, empty_u, -1 ./ P.tau_rec2 - q .* empty_u], ...
    npre, 2, 2);
B = [-b1 ./ P.tau_rel, zeros(npre, 1)];       % the rate enters the release term only
F = b1;
Fs = [ones(npre, 1), zeros(npre, 1)];
end

function ss = steady_state(r, P)
% (C4); see the header for the derivation and why it is written this way.
k1 = r .* P.tau_rec1 ./ P.tau_rel;
k2 = r .* P.tau_rec2 ./ (P.n2 .* P.tau_rel);
m = 1 + k1 - k2;
d = sqrt(m .^ 2 + 4 .* k1 .* k2);
denom = 1 + k1 + k2 + d;
b1 = 2 ./ denom;
b2 = (m + d) ./ denom;
low = m < 0;                                  % reserve nearly empty: m + d cancels
if any(low)
    b2(low) = 4 .* k1(low) .* k2(low) ./ ((d(low) - m(low)) .* denom(low));
end
ss.s = [b1, b2];
ss.F = b1;
end

%% Beyond the interface

function D = depression(b, P) %#ok<INUSD>
%DEPRESSION The depression factor from the b-block alone: the releasable pool.
D = b(:, 1, :);
end

%% Helpers

function tf = is_absent(cfg)
tf = isempty(cfg) || (isstruct(cfg) && isempty(fieldnames(cfg)));
end

function tf = is_positive(value)
tf = isnumeric(value) && isvector(value) && all(isfinite(value)) && all(value > 0);
end

function opts = fill_defaults(opts, defaults)
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = defaults.(names{k});
    end
end
end
