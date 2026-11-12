function model = route_product()
%ROUTE_PRODUCT Route model: product depression times output-only facilitation.
%
%   model = ROUTE_PRODUCT()
%
% The route model SRNNCellTypePairs has always had, behind the interface of
% route_model (read that header first). It serves every route without a
% `model` field: M depression components and L facilitation components,
% either count possibly zero, so also the route with no synapse dynamics
% at all (M = L = 0: no state, F = 1).
%
% THE EQUATIONS, per presynaptic neuron with rate r. The state is
% s = [b_1 .. b_M, g_1 .. g_L]:
%
%   db_m/dt = (1 - b_m)/tau_rec_m - b_m r/tau_rel_m           depression
%   dg_l/dt = (1 - g_l)/tau_dec_l + (G_l - g_l) r/tau_fac_l   facilitation
%
% Both are the same equation: a variable that relaxes to 1 at rest, with
% time constant tau_rest, and is driven by the rate toward a target, with
% time constant tau_drive / r. The target is 0 for depression and the
% ceiling G for facilitation:
%
%   ds_k/dt = (1 - s_k)/tau_rest_k + (target_k - s_k) r/tau_drive_k
%
% and that one line is what rhs evaluates. Nothing couples the components,
% so the linearisation A = d(ds)/ds is diagonal.
%
% THE OUTPUT FACTOR is F = D * prod_l g_l, with the depression factor D
% made from P_b = prod_m b_m in two optional steps:
%   combine     'product' (default): D = P_b. 'geomean': D = P_b^(1/M),
%               which at one ratio tau_rel/tau_rec for every component has
%               the one-component steady state 1/(1 + r tau_rec/tau_rel)
%               at every rate (the product raises it to the power M).
%   zero floor  (opts.zero_floor, the model class's std_zero_floor): D is
%               mapped linearly so that p_min -> 0 and 1 -> 1, where p_min
%               is D at the steady state of r = 1, the top of the logistic.
% Facilitation acts on the output only: g appears in no depression
% equation. (Release-coupled facilitation is another model.)
%
% CONFIGURATION (synapse_config.<pre>.<post>):
%   std.tau_rec   1 x M positive   recovery time constants
%   std.tau_rel   scalar or 1 x M  release time constants
%   std.combine   'product' | 'geomean'   (optional)
%   stf.tau_dec   1 x L positive   decay time constants
%   stf.tau_fac   scalar or 1 x L  build-up time constants
%   stf.G         scalar or 1 x L, >= 1   ceilings
%
% P, beyond the fields every model has (route_model):
%   tau_rest   npre x ns   [tau_rec, tau_dec]
%   tau_drive  npre x ns   [tau_rel, tau_fac]
%   target     npre x ns   [zeros,   G]
%   geomean    logical: the geometric-mean combine is in force (M > 1)
%   zero_floor logical: the zero floor is in force (M > 0)
%   p_min      npre x 1: the floor's p_min ([] when M = 0)
%   nominal    tau_rec, tau_rel, tau_dec, tau_fac, G as configured rows
%              ([] where absent), combine, zero_floor
%
% TWO HANDLES BEYOND THE INTERFACE, `combine` and `zero_floor`, are the
% two steps above. They exist so that SRNNCellTypePairs.combine_std_pair,
% apply_std_zero_floor_pair and std_zero_floor_pmin, which predate this
% file, are delegations to it and not second copies of the arithmetic.
% The model class's four hot functions do not call them.
%
% See also: route_model, SRNNCellTypePairs, test_route_models

model = struct('name', 'product', ...
    'compile', @compile, 'initial', @initial, 'rhs', @rhs, 'output', @output, ...
    'linearize', @linearize, 'steady_state', @steady_state, ...
    'combine', @combine, 'zero_floor', @zero_floor);
end

%% The interface

function P = compile(std_cfg, stf_cfg, opts)
if nargin < 3, opts = struct(); end
opts = fill_defaults(opts, struct('npre', 1, 'zero_floor', false, 'label', 'route', ...
    'err_id_prefix', 'route_product'));
err_id = [opts.err_id_prefix ':InvalidSynapseConfig'];
nominal = struct('tau_rec', [], 'tau_rel', [], 'tau_dec', [], 'tau_fac', [], 'G', [], ...
    'combine', 'product', 'zero_floor', logical(opts.zero_floor));

if ~is_absent(std_cfg)
    check_fields(std_cfg, {'tau_rec', 'tau_rel'}, {'combine'}, 'STD', opts.label, err_id);
    nominal.tau_rec = positive_row(std_cfg.tau_rec, [], 'tau_rec', 'STD', opts.label, err_id);
    nominal.tau_rel = positive_row(std_cfg.tau_rel, numel(nominal.tau_rec), 'tau_rel', 'STD', ...
        opts.label, err_id);
    if isfield(std_cfg, 'combine') && ~isempty(std_cfg.combine)
        mode = std_cfg.combine;
        is_text = (ischar(mode) && isrow(mode)) || (isstring(mode) && isscalar(mode));
        if ~is_text || ~any(strcmp(char(mode), {'product', 'geomean'}))
            error(err_id, 'STD combine on route %s must be ''product'' or ''geomean''.', opts.label);
        end
        nominal.combine = char(mode);
    end
end
if ~is_absent(stf_cfg)
    check_fields(stf_cfg, {'tau_dec', 'tau_fac', 'G'}, {}, 'STF', opts.label, err_id);
    nominal.tau_dec = positive_row(stf_cfg.tau_dec, [], 'tau_dec', 'STF', opts.label, err_id);
    L = numel(nominal.tau_dec);
    nominal.tau_fac = positive_row(stf_cfg.tau_fac, L, 'tau_fac', 'STF', opts.label, err_id);
    nominal.G = positive_row(stf_cfg.G, L, 'G', 'STF', opts.label, err_id);
    if any(nominal.G < 1)
        error(err_id, 'G on route %s must be at least one.', opts.label);
    end
end

P = struct('model', 'product');
P.n_b = numel(nominal.tau_rec);
P.n_g = numel(nominal.tau_dec);
P.ns = P.n_b + P.n_g;
P.npre = opts.npre;
P.nominal = nominal;
P.signature = nominal;
P.signature.model = P.model;
P.signature.zero_floor = nominal.zero_floor && P.n_b > 0;    % in force, not merely asked for

% Per neuron, the nominal row replicated.
P.tau_rest  = per_neuron([nominal.tau_rec, nominal.tau_dec], opts.npre);
P.tau_drive = per_neuron([nominal.tau_rel, nominal.tau_fac], opts.npre);
P.target    = per_neuron([zeros(1, P.n_b), nominal.G], opts.npre);
P.geomean = strcmp(nominal.combine, 'geomean') && P.n_b > 1;
P.coupling = logical(eye(P.ns));      % no component's derivative depends on another

% The floor's p_min: the depression factor, before the floor, at the
% steady state of r = 1.
P.zero_floor = false;
P.p_min = [];
if P.n_b > 0
    at_ceiling = steady_state(1, P);
    P.p_min = combine(prod(at_ceiling.s(:, 1:P.n_b), 2), P);
    P.zero_floor = nominal.zero_floor;
end
end

function s0 = initial(P)
% Every component starts at 1, its value at rest.
s0 = ones(1, P.ns);
end

function [ds, F] = rhs(s, r, P)
ds = (1 - s) ./ P.tau_rest + ((P.target - s) .* r) ./ P.tau_drive;
F = output(s, P);
end

function F = output(s, P)
if P.n_g == 0                 % every column is a b
    F = depression(prod(s, 2), P);
else
    F = depression(prod(s(:, 1:P.n_b, :), 2), P) .* prod(s(:, P.n_b + 1:end, :), 2);
end
end

function [A, B, F, Fs] = linearize(s, r, P)
npre = size(s, 1);
ns = P.ns;
A = zeros(npre, ns * ns);
A(:, 1:ns + 1:end) = -1 ./ P.tau_rest - r ./ P.tau_drive;    % the diagonal (m, m)
A = reshape(A, npre, ns, ns);
B = (P.target - s) ./ P.tau_drive;

% F = D(P_b) prod(g), so
%   dF/db_m = (dD/dP_b) (the product of the other b) prod(g)
%   dF/dg_l = D (the product of the other g)
if P.n_g == 0                 % every column is a b: the same lines without g
    [F, dD_dP] = depression(prod(s, 2), P);
    Fs = dD_dP .* product_of_others(s);
else
    b = s(:, 1:P.n_b);
    g = s(:, P.n_b + 1:end);
    [D, dD_dP] = depression(prod(b, 2), P);
    facilitation = prod(g, 2);
    F = D .* facilitation;
    Fs = [(dD_dP .* facilitation) .* product_of_others(b), D .* product_of_others(g)];
end
end

function ss = steady_state(r, P)
% ds/dt = 0 at constant r: s = (1 + target k r) / (1 + k r), k = tau_rest/tau_drive.
% Depression (target 0): 1/(1 + r tau_rec/tau_rel). Facilitation (target
% G): (1 + G k r)/(1 + k r).
k = P.tau_rest ./ P.tau_drive;
ss.s = (1 + P.target .* k .* r) ./ (1 + k .* r);
ss.F = output(ss.s, P);
end

%% The depression factor, and its two steps

function [D, dD_dP] = depression(P_b, P)
% D from P_b = prod_m b_m (the combine, then the floor), and dD/dP_b.
D = P_b;
dD_dP = 1;
if P.geomean
    [D, dD_dP] = combine(P_b, P);
end
if P.zero_floor
    [D, floor_gain] = zero_floor(D, P);
    dD_dP = floor_gain .* dD_dP;
end
end

function [D, dD_dP] = combine(P_b, P)
%COMBINE The depression factor D from P_b = prod_m b_m, and dD/dP_b.
% b_m > 0 always (db/dt = 1/tau_rec at b = 0), so P_b > 0.
if P.geomean
    D = P_b .^ (1 / P.n_b);
    dD_dP = D ./ (P.n_b .* P_b);
else
    D = P_b;
    dD_dP = 1;
end
end

function [D, gain] = zero_floor(D, P)
%ZERO_FLOOR The floor's linear map of D, and its slope.
if P.zero_floor
    D = (D - P.p_min) ./ (1 - P.p_min);
    gain = 1 ./ (1 - P.p_min);
else
    gain = 1;
end
end

%% Helpers

function Q = product_of_others(v)
% Q(:, m) = the product of every column of v except m (ones where v has
% one column; npre x 0 where it has none).
% One and two columns, the usual counts, are written out: this is inside
% every jacobian_times call.
k = size(v, 2);
if k == 1
    Q = ones(size(v));
elseif k == 2
    Q = v(:, [2 1]);
else
    Q = zeros(size(v));
    for m = 1:k
        Q(:, m) = prod(v(:, [1:m - 1, m + 1:k]), 2);
    end
end
end

function tf = is_absent(cfg)
tf = isempty(cfg) || (isstruct(cfg) && isempty(fieldnames(cfg)));
end

function opts = fill_defaults(opts, defaults)
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts, names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = defaults.(names{k});
    end
end
end

function check_fields(cfg, required, optional, mechanism, label, err_id)
if ~isstruct(cfg) || ~isscalar(cfg)
    error(err_id, '%s configuration on route %s must be a scalar struct.', mechanism, label);
end
% (Plain loops: this runs for every route on every get_params, and setdiff
% on two short cell arrays costs more than the rest of compile.)
fields = fieldnames(cfg);
known = [required, optional];
for k = 1:numel(required)
    if ~any(strcmp(required{k}, fields))
        missing = setdiff(required, fields);
        error(err_id, 'Missing %s field on route %s: %s.', mechanism, label, missing{1});
    end
end
for k = 1:numel(fields)
    if ~any(strcmp(fields{k}, known))
        unknown = setdiff(fields, known);
        error(err_id, 'Unknown %s field on route %s: %s.', mechanism, label, unknown{1});
    end
end
end

function M = per_neuron(row, npre)
% The nominal row, one copy per presynaptic neuron (npre x 0 when empty).
M = repmat(reshape(row, 1, []), npre, 1);
end

function value = positive_row(value, count, field_name, mechanism, label, err_id)
% A row of positive finite values. COUNT empty: any nonempty length (this
% field sets the component count). COUNT given: a scalar is repeated COUNT
% times, anything else must have COUNT values.
if ~isnumeric(value) || ~isvector(value) || any(~isfinite(value)) || any(value <= 0)
    error(err_id, '%s on route %s must contain positive finite values.', field_name, label);
end
value = reshape(value, 1, []);
if isempty(count)
    if isempty(value)
        error(err_id, '%s %s must be nonempty on route %s.', mechanism, field_name, label);
    end
elseif isscalar(value)
    value = repmat(value, 1, count);
elseif numel(value) ~= count
    error(err_id, '%s on route %s must be scalar or have %d values.', field_name, label, count);
end
end
