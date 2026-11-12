function model = route_pool()
%ROUTE_POOL Route model: vesicle pools in series, with release-coupled facilitation.
%
%   model = ROUTE_POOL()
%
% The vesicle-pool family behind the interface of route_model (read that
% header first). Registered as 'pool'; a route gets it by asking for
% synapse_config.<pre>.<post>.std.model = 'chain' (chained depression) or
% .stf.model = 'chain' (release-coupled facilitation), or both.
%
% THE EQUATIONS, per presynaptic neuron with rate r. The full state is
% s = [b1, b2, g, h]:
%   b1  the RELEASABLE pool, as a fraction of its capacity      1 at rest
%   b2  the RESERVE that refills it, likewise                   1 at rest
%   g   the release probability, relative to its resting value  1 at rest
%   h   a calcium trace, a low-pass copy of the rate            0 at rest
%
%   db1/dt = b2 (1 - b1)/tau_rec1 - g b1 r/tau_rel                    (C1)
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2) b2 (1 - b1)/tau_rec1          (C2)
%   dg/dt  = (1 - g)/tau_dec + (G - g) d/tau_fac                      (C5)
%   dh/dt  = (r - h)/tau_h                                            (C6)
%   F      = g b1                                                     (C3)
% with the drive d = h.
%
% Term by term: the releasable pool empties by release, g b1 r/tau_rel,
% and is refilled from the reserve at a rate proportional to what the
% reserve holds (b2) and to the releasable pool's empty sites (1 - b1).
% That refill flux is what the reserve loses, divided by n2 because the
% reserve's capacity is n2 times the releasable pool's (the same vesicles
% are a smaller fraction of a larger pool). The reserve itself refills
% toward 1 with tau_rec2. Only the releasable pool releases, so only b1 is
% a depression factor; the reserve acts on the output through b1 alone.
% The release probability g relaxes to 1 with tau_dec and is driven toward
% its ceiling G by the drive d; it multiplies BOTH the output and the rate
% at which the releasable pool is used (that is the coupling: more release
% per spike empties the pool faster). The calcium trace follows the rate
% with tau_h and stands between the rate and g.
%
% A ROUTE CARRIES ONLY THE VARIABLES IT HAS. NONE IS EVER CLAMPED: an
% absent variable is removed from the state and from the equations.
%   no reserve       b2 is dropped and (C1)'s refill reads (1 - b1)/tau_rec1
%                    (the reserve held full; n2 -> infinity is that limit)
%   no calcium trace h is dropped and the rate drives g: d = r
%   no facilitation  g and h are dropped; the release term is b1 r/tau_rel
%                    and F = b1
%   no depression    b1 and b2 are dropped and F = g
% That gives eight variants, the b-block columns [b1, b2] crossed with the
% g-block columns [g, h]:
%
%   depression \ facilitation   none          g             g, h
%   none                        (no state:    [g]           [g, h]
%                                route_product)
%   single pool                 [b1]          [b1, g]       [b1, g, h]
%   chain                       [b1, b2]      [b1, b2, g]   [b1, b2, g, h]
%
% The single pool with no facilitation is route_product's one-component
% route (the same equation) and SRNNCellTypePairs gives such a route to
% route_product; compile serves it here too, for whoever calls it.
%
% RANGE. b1, b2 in [0, 1], g in [1, G], h in [0, 1] for r in [0, 1]: at
% b1 = 0 its derivative is b2/tau_rec1 >= 0 and at b1 = 1 it is
% -g r/tau_rel <= 0; at b2 = 0 its derivative is 1/tau_rec2 > 0 and at
% b2 = 1 it is -(1 - b1)/(n2 tau_rec1) <= 0; at g = 1 its derivative is
% (G - 1) d/tau_fac >= 0 and at g = G it is (1 - G)/tau_dec <= 0; at h = 0
% its derivative is r/tau_h >= 0 and at h = 1 it is (r - 1)/tau_h <= 0.
%
% THE LINEARISATION, with u = 1/tau_rec1, v = 1/tau_rec2, w = 1/tau_rel,
% q = 1/n2 (A(i, m, k) = d(ds_m)/d(s_k), B = d(ds)/dr, Fs = dF/ds), for the
% route with all four variables:
%
%            d/db1          d/db2             d/dg               d/dh           d/dr
%   db1/dt  -b2 u - g r w   (1 - b1) u       -b1 r w              0            -g b1 w
%   db2/dt   q b2 u         -v - q (1 - b1) u  0                  0              0
%   dg/dt    0               0               -1/tau_dec - d/tau_fac  (G - g)/tau_fac  0
%   dh/dt    0               0                0                  -1/tau_h       1/tau_h
%   F        g               0                b1                  0
%
% Without the reserve: delete the b2 row and column and put b2 = 1.
% Without the calcium trace: delete the h row and column, d = r, and the
% rate column of dg/dt is (G - g)/tau_fac. Without facilitation: delete
% the g and h rows and columns and put g = 1 (dF/db1 = 1). Without
% depression: delete both b rows and columns, F = g, dF/dg = 1.
% The block is NOT symmetric and not diagonal (except the lone g), and
% P.coupling states which entries can be nonzero for the variant. The one
% entry that is new in kind is (b1, g): a facilitation variable inside a
% depression equation, which is why a route is one system.
%
% THE STEADY STATE at a constant rate r. h = r, and (C5) = 0 gives, with
% or without the calcium trace,
%
%   g = (1 + G kappa r)/(1 + kappa r),   kappa = tau_dec/tau_fac
%
% The pools then see the EFFECTIVE RATE rho = g r (rho = r without
% facilitation). With k1 = rho tau_rec1/tau_rel,
% k2 = rho tau_rec2/(n2 tau_rel) and m = 1 + k1 - k2, (C2) = 0 gives
% b2 = 1 - k2 b1 and (C1) = 0 then k2 b1^2 - (1 + k1 + k2) b1 + 1 = 0, whose
% root in [0, 1] is
%
%   b1 = 2 / (1 + k1 + k2 + d),   d = sqrt(m^2 + 4 k1 k2)             (C4)
%   b2 = 1 - k2 b1
%
% (d^2 is (1 + k1 + k2)^2 - 4 k2 written as a sum of two nonnegative
% terms.) The quadratic formula's other form, (s - d)/(2 k2), cancels at
% small k2; this one does not, and it reduces to the single pool's
% 1/(1 + k1) at k2 = 0, which is what a route without a reserve returns.
% b2 is evaluated as (m + d)/(1 + k1 + k2 + d), or as
% 4 k1 k2/((d - m)(1 + k1 + k2 + d)) where m < 0, which is the same number
% without the cancellation of 1 - k2 b1 when the reserve is nearly empty.
% The sustained output r g b1 = rho b1(rho) rises with the rate to the
% ceiling 1/(tau_rec1/tau_rel + tau_rec2/(n2 tau_rel)) (the product's
% falls). Coupled facilitation therefore stretches the rate axis and
% nothing else at steady state: the output reaches the ceiling at lower
% rates, and the ceiling does not move.
%
% CONFIGURATION (synapse_config.<pre>.<post>).
%   std, the chain:
%     model     'chain'
%     tau_rec   1 x 2 positive: [releasable-pool refill, reserve refill]
%     tau_rel   positive scalar: ONE release time, the releasable pool's
%     n2        positive finite scalar: reserve capacity over
%               releasable-pool capacity. No default.
%   std, the single pool (only beside stf.model = 'chain', when the route
%   comes from SRNNCellTypePairs): tau_rec and tau_rel, positive scalars,
%   as a one-component product route gives them (no `model`, or 'product').
%   std absent: no depression.
%   stf, release-coupled:
%     model     'chain'
%     tau_dec   positive scalar: decay of g to 1
%     tau_fac   positive scalar: build-up of g
%     G         scalar >= 1: the ceiling of g
%     tau_h     positive finite scalar, OPTIONAL: the calcium trace's time
%               constant. Absent: no calcium trace, the rate drives g.
%   stf absent: no facilitation.
% Refused, each with its own identifier (<prefix> is opts.err_id_prefix):
%   <prefix>:ChainTauRecCount     chain: tau_rec has not exactly two values
%   <prefix>:ChainTauRelNotScalar chain: tau_rel is not a scalar
%   <prefix>:ChainMissingN2       chain: n2 is absent
%   <prefix>:ChainInvalidN2       chain: n2 is not a positive finite scalar
%   <prefix>:ChainCombine         std.combine (the product's option)
%   <prefix>:ChainZeroFloor       opts.zero_floor (the class's
%                                 std_zero_floor) on a route with depression
%   <prefix>:ChainStfOutputOnly   facilitation without stf.model = 'chain'
%                                 (today's output-only form) on a pool route
%   <prefix>:ChainStfProductStd   depression that is a product of two or
%                                 more components (neither the single pool
%                                 nor the chain)
%   <prefix>:ChainStfNotScalar    tau_dec, tau_fac or G is not a scalar
%   <prefix>:ChainStfInvalidG     G < 1
%   <prefix>:ChainStfInvalidTauH  tau_h is not a positive finite scalar
%   <prefix>:InvalidSynapseConfig anything else malformed
% The zero floor on a route WITHOUT depression is not in force and not
% refused, as on a product route: it has nothing to map.
%
% P, beyond the fields every model has (route_model):
%   has_g, has_h   logical: the route carries g; it carries h
%   s0             1 x ns, the starting state
%   PER NEURON, npre x 1, only those the variant has:
%     tau_rec1, tau_rel        depression
%     tau_rec2, n2             the reserve
%     tau_dec, tau_fac, G      facilitation
%     tau_h                    the calcium trace
%   nominal: tau_rec (1 x 2, scalar or []), tau_rel, n2, tau_dec, tau_fac,
%     G, tau_h as configured ([] where absent), std_model ('chain',
%     'single' or 'none') and stf_model ('chain' or '').
%   signature: the nominal values and the model's name. A route without
%     facilitation states no facilitation fields at all (stf_model, tau_h).
%
% TWO HANDLES FOR DISPLAY CODE, which holds a route's b array or g array
% and not its whole state:
%   D = depression(b, P)     (route_model) the depression factor: b1.
%   Fg = facilitation(g, P)  the factor facilitation puts on the output,
%                            from the route's g-block (npre x n_g, or
%                            npre x n_g x nt; Fg has the shape of
%                            prod(g, 2)): g, the first column. NOT g h:
%                            the calcium trace acts only through g.
%                            SRNNCellTypePairs.stf_route_factor calls it.
%
% THE PER-NEURON SPREAD (route_model: "THE SPREAD"; FD-015). With
% opts.std_sigma or opts.stf_sigma above zero the per-neuron columns hold
% drawn values and not the nominal: value = nominal * exp(sigma * z), z the
% value's column of the presynaptic neuron's standard normals. The columns
% this model reads:
%   depression (opts.std_z)     tau_rec1 1, tau_rel 2, tau_rec2 3, n2 5
%                               (column 4 is a product route's tau_rel_2)
%   facilitation (opts.stf_z)   tau_dec 1, tau_fac 2, tau_h 3, and 4 for
%                               the ceiling, jittered on its excess:
%                               G = 1 + (G_nominal - 1) * exp(sigma * z),
%                               so G >= 1 for every neuron
% The single pool reads the chain's columns 1 and 2, so a neuron has the
% same tau_rec1 and tau_rel whether its route is a single pool or a chain
% (and the same as on a product route's first component): equation (C1)
% is the same between them neuron by neuron.
% THE ORDERING RULE, chain only, where the nominal reserve is the slower
% pool (tau_rec2 > tau_rec1): a neuron whose log(tau_rec2 / tau_rec1)
% falls below a quarter of the nominal log-span has ITS RESERVE REFILL
% tau_rec2 raised until the gap equals that floor. tau_rec1 is never
% moved: it must stay the value the single pool has. P.clipped (npre x 1
% logical) says which neurons were moved. At sigma = 0 every column is
% the nominal replicated, exactly, and nothing is clipped.
%
% See also: route_model, route_product, SRNNCellTypePairs,
%           test_route_models, test_std_chain_model, test_stf_chain_model

model = struct('name', 'pool', ...
    'compile', @compile, 'initial', @initial, 'rhs', @rhs, 'output', @output, ...
    'linearize', @linearize, 'steady_state', @steady_state, ...
    'depression', @depression, 'facilitation', @facilitation);
end

%% The interface

function P = compile(std_cfg, stf_cfg, opts)
if nargin < 3, opts = struct(); end
opts = fill_defaults(opts, struct('npre', 1, 'zero_floor', false, 'label', 'route', ...
    'err_id_prefix', 'route_pool', 'std_sigma', 0, 'stf_sigma', 0));
prefix = opts.err_id_prefix;
label = opts.label;
err_id = [prefix ':InvalidSynapseConfig'];
has_std = ~is_absent(std_cfg);
has_stf = ~is_absent(stf_cfg);

% Which facilitation. Only the release-coupled form lives here.
if has_stf
    if ~isstruct(stf_cfg) || ~isscalar(stf_cfg)
        error(err_id, 'STF configuration on route %s must be a scalar struct.', label);
    end
    if ~(isfield(stf_cfg, 'model') && is_name(stf_cfg.model, 'chain'))
        error([prefix ':ChainStfOutputOnly'], ...
            ['Route %s is a vesicle-pool route (std.model = ''chain'') with output-only ' ...
            'facilitation (stf without a model, or stf.model = ''product''). That form is ' ...
            'not offered on a pool route: release probability also sets how fast the ' ...
            'releasable pool is used. Ask for the release-coupled form with ' ...
            'stf.model = ''chain''.'], label);
    end
end
% Which depression: the chain, the single pool, or none.
std_model = 'none';
if has_std
    if opts.zero_floor
        error([prefix ':ChainZeroFloor'], ...
            ['std_zero_floor cannot be used with a vesicle-pool route (route %s has ' ...
            'std.model = ''chain'' or stf.model = ''chain''). The floor maps a product ' ...
            'route''s factor at r = 1 to zero; a pool route''s factor is the releasable ' ...
            'pool itself and has no such map.'], label);
    end
    if ~isstruct(std_cfg) || ~isscalar(std_cfg)
        error(err_id, 'STD configuration on route %s must be a scalar struct.', label);
    end
    fields = fieldnames(std_cfg);
    if any(strcmp('combine', fields))
        error([prefix ':ChainCombine'], ...
            ['std.combine on route %s: combine belongs to product depression and has no ' ...
            'meaning on a vesicle-pool route (std.model = ''chain'' or stf.model = ' ...
            '''chain''), whose depression factor is the releasable pool b1.'], label);
    end
    std_model = 'single';
    if any(strcmp('model', fields))
        if is_name(std_cfg.model, 'chain')
            std_model = 'chain';
        elseif ~is_name(std_cfg.model, 'product')
            error(err_id, ['route_pool compiles std.model = ''chain'', or one depression ' ...
                'component without a model; route %s states ''%s''.'], label, ...
                char(string(std_cfg.model)));
        end
    end
elseif ~has_stf
    error(err_id, ['Route %s has neither depression nor facilitation: it has no state, ' ...
        'and route_product serves it.'], label);
end

nominal = struct('tau_rec', [], 'tau_rel', [], 'tau_dec', [], 'tau_fac', [], 'G', [], ...
    'n2', [], 'std_model', std_model, 'stf_model', '', 'tau_h', []);

% Depression.
if has_std
    if strcmp(std_model, 'chain') && ~any(strcmp('n2', fields))
        error([prefix ':ChainMissingN2'], ...
            ['Chained depression on route %s needs std.n2, the reserve''s capacity relative ' ...
            'to the releasable pool''s. It has no default.'], label);
    end
    for k = {'tau_rec', 'tau_rel'}
        if ~any(strcmp(k{1}, fields))
            error(err_id, 'Missing STD field on route %s: %s.', label, k{1});
        end
    end
    known = {'model', 'tau_rec', 'tau_rel'};
    if strcmp(std_model, 'chain'), known = [known, {'n2'}]; end
    unknown = setdiff(fields, known);
    if ~isempty(unknown)
        error(err_id, 'Unknown STD field on pool route %s: %s.', label, unknown{1});
    end
    tau_rec = std_cfg.tau_rec;
    tau_rel = std_cfg.tau_rel;
    if ~is_positive(tau_rec)
        error(err_id, 'tau_rec on route %s must contain positive finite values.', label);
    end
    if ~is_positive(tau_rel)
        error(err_id, 'tau_rel on route %s must contain positive finite values.', label);
    end
    if strcmp(std_model, 'chain')
        if numel(tau_rec) ~= 2
            error([prefix ':ChainTauRecCount'], ...
                ['Chained depression on route %s needs exactly two tau_rec values, ' ...
                '[releasable-pool refill, reserve refill]; got %d.'], label, numel(tau_rec));
        end
        if ~isscalar(tau_rel)
            error([prefix ':ChainTauRelNotScalar'], ...
                ['Chained depression on route %s needs a scalar tau_rel: only the releasable ' ...
                'pool releases. Got %d values.'], label, numel(tau_rel));
        end
        n2 = std_cfg.n2;
        if ~(is_positive(n2) && isscalar(n2))
            error([prefix ':ChainInvalidN2'], ...
                'n2 on route %s must be a positive finite scalar.', label);
        end
        nominal.n2 = double(n2);
    else
        if numel(tau_rec) ~= 1
            error([prefix ':ChainStfProductStd'], ...
                ['Route %s has release-coupled facilitation (stf.model = ''chain'') beside a ' ...
                'product of %d depression components. Release probability multiplies the ' ...
                'use of ONE releasable pool: give the route one depression component, the ' ...
                'chain (std.model = ''chain''), or none.'], label, numel(tau_rec));
        end
        if ~isscalar(tau_rel)
            error(err_id, 'tau_rel on route %s must be a scalar: the route has one pool.', label);
        end
    end
    nominal.tau_rec = reshape(double(tau_rec), 1, []);
    nominal.tau_rel = double(tau_rel);
end

% Facilitation.
if has_stf
    stf_fields = fieldnames(stf_cfg);
    for k = {'tau_dec', 'tau_fac', 'G'}
        if ~any(strcmp(k{1}, stf_fields))
            error(err_id, 'Missing STF field on route %s: %s.', label, k{1});
        end
    end
    unknown = setdiff(stf_fields, {'model', 'tau_dec', 'tau_fac', 'G', 'tau_h'});
    if ~isempty(unknown)
        error(err_id, 'Unknown STF field on pool route %s: %s.', label, unknown{1});
    end
    for k = {'tau_dec', 'tau_fac', 'G'}
        value = stf_cfg.(k{1});
        if ~is_positive(value)
            error(err_id, '%s on route %s must contain positive finite values.', k{1}, label);
        end
        if ~isscalar(value)
            error([prefix ':ChainStfNotScalar'], ...
                ['Release-coupled facilitation on route %s has ONE release-probability ' ...
                'variable: %s must be a scalar. Got %d values.'], label, k{1}, numel(value));
        end
    end
    if stf_cfg.G < 1
        error([prefix ':ChainStfInvalidG'], 'G on route %s must be at least one.', label);
    end
    nominal.tau_dec = double(stf_cfg.tau_dec);
    nominal.tau_fac = double(stf_cfg.tau_fac);
    nominal.G = double(stf_cfg.G);
    nominal.stf_model = 'chain';
    if any(strcmp('tau_h', stf_fields))
        tau_h = stf_cfg.tau_h;
        if ~(is_positive(tau_h) && isscalar(tau_h))
            error([prefix ':ChainStfInvalidTauH'], ...
                ['tau_h on route %s must be a positive finite scalar. Leave the field out ' ...
                'for a route without a calcium trace.'], label);
        end
        nominal.tau_h = double(tau_h);
    end
end

P = struct('model', 'pool');
P.n_b = numel(nominal.tau_rec);           % b1 (column 1) and, on a chain, b2 (column 2)
P.has_g = has_stf;
P.has_h = ~isempty(nominal.tau_h);
P.n_g = P.has_g + P.has_h;                % g (column n_b + 1) and h (column n_b + 2)
P.ns = P.n_b + P.n_g;
P.npre = opts.npre;
P.nominal = nominal;
% The signature: what the route's dynamics are. A route without
% facilitation states no facilitation fields.
P.signature = rmfield(nominal, {'stf_model', 'tau_h'});
if has_stf
    P.signature.stf_model = nominal.stf_model;
    P.signature.tau_h = nominal.tau_h;
end
P.signature.model = P.model;

% Which variable's derivative can depend on which (rows: derivative).
coupling = false(4);                      % of [b1, b2, g, h]
coupling(1, [1 2 3]) = true;              % db1/dt: b1, b2 (refill), g (release)
coupling(2, [1 2]) = true;                % db2/dt: b1 (refill), b2
coupling(3, [3 4]) = true;                % dg/dt: g, h (the drive)
coupling(4, 4) = true;                    % dh/dt: h
present = [P.n_b >= 1, P.n_b == 2, P.has_g, P.has_h];
P.coupling = coupling(present, present);

% Pools start full and release probability at rest, their values at
% r = 0; the calcium trace starts at 0, likewise.
P.s0 = [ones(1, P.n_b + P.has_g), zeros(1, P.has_h)];

% Per neuron, the nominal value replicated; only what the variant has.
npre = opts.npre;
if P.n_b >= 1
    P.tau_rec1 = repmat(nominal.tau_rec(1), npre, 1);
    P.tau_rel  = repmat(nominal.tau_rel, npre, 1);
end
if P.n_b == 2
    P.tau_rec2 = repmat(nominal.tau_rec(2), npre, 1);
    P.n2       = repmat(nominal.n2, npre, 1);
end
if P.has_g
    P.tau_dec = repmat(nominal.tau_dec, npre, 1);
    P.tau_fac = repmat(nominal.tau_fac, npre, 1);
    P.G       = repmat(nominal.G, npre, 1);
end
if P.has_h
    P.tau_h = repmat(nominal.tau_h, npre, 1);
end

% The per-neuron spread: each value times exp(sigma z) from its own column
% of the presynaptic neuron's draws (the header says which); of the
% ceiling, the excess over 1. Without a spread nothing here runs, and the
% values stay the nominal exactly.
z_std = spread_draws(opts, 'std', P.n_b >= 1, 2 + 3 * (P.n_b == 2));
z_stf = spread_draws(opts, 'stf', P.has_g, 4);
P.jittered = ~isempty(z_std) || ~isempty(z_stf);
P.clipped = false(npre, 1);
if ~isempty(z_std)
    e = exp(opts.std_sigma .* z_std);
    P.tau_rec1 = P.tau_rec1 .* e(:, 1);
    P.tau_rel  = P.tau_rel .* e(:, 2);
    if P.n_b == 2
        P.tau_rec2 = P.tau_rec2 .* e(:, 3);
        P.n2       = P.n2 .* e(:, 5);
        % The ordering rule (header): the reserve stays the slower pool by
        % at least a quarter of the nominal log-span. Only tau_rec2 moves.
        span = log(nominal.tau_rec(2)) - log(nominal.tau_rec(1));
        if span > 0
            floor_gap = span / 4;
            gap = span + opts.std_sigma .* (z_std(:, 3) - z_std(:, 1));
            short = gap < floor_gap;
            P.tau_rec2(short) = P.tau_rec1(short) .* exp(floor_gap);
            P.clipped = short;
        end
    end
end
if ~isempty(z_stf)
    e = exp(opts.stf_sigma .* z_stf);
    P.tau_dec = P.tau_dec .* e(:, 1);
    P.tau_fac = P.tau_fac .* e(:, 2);
    P.G       = 1 + (P.G - 1) .* e(:, 4);
    if P.has_h
        P.tau_h = P.tau_h .* e(:, 3);
    end
end
end

function s0 = initial(P)
s0 = P.s0;
end

function [ds, F] = rhs(s, r, P)
n_b = P.n_b;
if P.has_g
    g = s(:, n_b + 1);
    if P.has_h
        h = s(:, n_b + 2);
        dg = [(1 - g) ./ P.tau_dec + (P.G - g) .* h ./ P.tau_fac, (r - h) ./ P.tau_h];
    else
        dg = (1 - g) ./ P.tau_dec + (P.G - g) .* r ./ P.tau_fac;
    end
    if n_b == 0                               % no depression: a facilitating synapse
        ds = dg;
        F = g;
        return;
    end
end
b1 = s(:, 1);
if P.has_g
    release = g .* b1 .* r ./ P.tau_rel;      % release probability speeds the pool's use
else
    release = b1 .* r ./ P.tau_rel;
end
if n_b == 2
    b2 = s(:, 2);
    refill = b2 .* (1 - b1) ./ P.tau_rec1;    % reserve -> releasable pool, in units of the releasable pool
    db = [refill - release, (1 - b2) ./ P.tau_rec2 - refill ./ P.n2];
else
    db = (1 - b1) ./ P.tau_rec1 - release;    % no reserve: the refill as if b2 = 1
end
if P.has_g
    ds = [db, dg];
    F = g .* b1;
else
    ds = db;
    F = b1;
end
end

function F = output(s, P)
if P.n_b == 0
    F = s(:, 1, :);                           % g
elseif P.has_g
    F = s(:, P.n_b + 1, :) .* s(:, 1, :);     % g b1
else
    F = s(:, 1, :);                           % b1
end
end

function [A, B, F, Fs] = linearize(s, r, P)
% A(i, m, k) = d(ds_m)/d(s_k) is written at column (k - 1) ns + m of an
% npre x ns^2 array and reshaped at the end. Entries not written are zero.
npre = size(s, 1);
ns = P.ns;
n_b = P.n_b;
A = zeros(npre, ns * ns);
B = zeros(npre, ns);
Fs = zeros(npre, ns);
if P.has_g
    cg = n_b + 1;                             % g's column of s
    g = s(:, cg);
    at_gg = (cg - 1) * ns + cg;
    if P.has_h
        ch = cg + 1;                          % h's column of s
        h = s(:, ch);
        A(:, at_gg) = -1 ./ P.tau_dec - h ./ P.tau_fac;          % (g, g)
        A(:, (ch - 1) * ns + cg) = (P.G - g) ./ P.tau_fac;       % (g, h)
        A(:, (ch - 1) * ns + ch) = -1 ./ P.tau_h;                % (h, h)
        B(:, ch) = 1 ./ P.tau_h;                                 % the rate drives h, and g through it
    else
        A(:, at_gg) = -1 ./ P.tau_dec - r ./ P.tau_fac;          % (g, g)
        B(:, cg) = (P.G - g) ./ P.tau_fac;                       % the rate drives g
    end
    if n_b == 0
        F = g;
        Fs(:, cg) = 1;
        A = reshape(A, npre, ns, ns);
        return;
    end
end
b1 = s(:, 1);
if P.has_g
    use = g .* r ./ P.tau_rel;                % -d(release)/d(b1)
    A(:, (cg - 1) * ns + 1) = -b1 .* r ./ P.tau_rel;             % (b1, g)
    B(:, 1) = -g .* b1 ./ P.tau_rel;          % the rate enters the release term
    F = g .* b1;
    Fs(:, 1) = g;
    Fs(:, cg) = b1;
else
    use = r ./ P.tau_rel;
    B(:, 1) = -b1 ./ P.tau_rel;
    F = b1;
    Fs(:, 1) = 1;
end
if n_b == 2
    b2 = s(:, 2);
    u = 1 ./ P.tau_rec1;
    q = 1 ./ P.n2;
    empty_u = (1 - b1) .* u;                  % d(refill)/d(b2)
    full_u = b2 .* u;                         % -d(refill)/d(b1)
    A(:, 1) = -full_u - use;                                     % (b1, b1)
    A(:, 2) = q .* full_u;                                       % (b2, b1)
    A(:, ns + 1) = empty_u;                                      % (b1, b2)
    A(:, ns + 2) = -1 ./ P.tau_rec2 - q .* empty_u;              % (b2, b2)
else
    A(:, 1) = -1 ./ P.tau_rec1 - use;                            % (b1, b1)
end
A = reshape(A, npre, ns, ns);
end

function ss = steady_state(r, P)
% See the header: h = r, g from (C5) = 0, the pools from (C4) at the
% effective rate rho = g r.
rho = r;
if P.has_g
    kappa = P.tau_dec ./ P.tau_fac;
    g = (1 + P.G .* kappa .* r) ./ (1 + kappa .* r);
    rho = g .* r;
    g_block = g;
    if P.has_h
        g_block = [g, r .* ones(size(g))];
    end
end
if P.n_b == 2
    k1 = rho .* P.tau_rec1 ./ P.tau_rel;
    k2 = rho .* P.tau_rec2 ./ (P.n2 .* P.tau_rel);
    m = 1 + k1 - k2;
    d = sqrt(m .^ 2 + 4 .* k1 .* k2);
    denom = 1 + k1 + k2 + d;
    b1 = 2 ./ denom;
    b2 = (m + d) ./ denom;
    low = m < 0;                              % reserve nearly empty: m + d cancels
    if any(low)
        b2(low) = 4 .* k1(low) .* k2(low) ./ ((d(low) - m(low)) .* denom(low));
    end
    b_block = [b1, b2];
elseif P.n_b == 1
    b1 = 1 ./ (1 + rho .* P.tau_rec1 ./ P.tau_rel);
    b_block = b1;
end
if P.n_b == 0
    ss.s = g_block;
    ss.F = g;
elseif P.has_g
    ss.s = [b_block, g_block];
    ss.F = g .* b1;
else
    ss.s = b_block;
    ss.F = b1;
end
end

%% For display code

function D = depression(b, P) %#ok<INUSD>
%DEPRESSION The depression factor from the b-block alone: the releasable pool.
D = b(:, 1, :);
end

function Fg = facilitation(g, P) %#ok<INUSD>
%FACILITATION The factor facilitation puts on the output, from the g-block
% alone: the release probability g, its first column. Not g h.
Fg = g(:, 1, :);
end

%% Helpers

function tf = is_absent(cfg)
tf = isempty(cfg) || (isstruct(cfg) && isempty(fieldnames(cfg)));
end

function tf = is_name(value, name)
tf = ((ischar(value) && isrow(value)) || (isstring(value) && isscalar(value))) && ...
    strcmp(char(value), name);
end

function z = spread_draws(opts, mechanism, present, n_cols)
% The route's rows of its presynaptic type's standard normals for one
% mechanism ('std' or 'stf'), or [] when that mechanism has no spread on
% this route (sigma zero, or the route has no such variables).
z = [];
sigma = opts.([mechanism '_sigma']);
if ~present || sigma == 0
    return;
end
name = [mechanism '_z'];
if ~isfield(opts, name) || size(opts.(name), 1) ~= opts.npre || size(opts.(name), 2) < n_cols
    error([opts.err_id_prefix ':SpreadDraws'], ...
        ['Route %s has a per-neuron spread (opts.%s_sigma = %g) and needs opts.%s: one row of ' ...
        'standard normals per presynaptic neuron (%d rows) with at least %d columns.'], ...
        opts.label, mechanism, sigma, name, opts.npre, n_cols);
end
z = opts.(name);
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
