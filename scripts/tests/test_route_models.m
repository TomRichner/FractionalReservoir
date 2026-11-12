% test_route_models.m - every route model, alone, with no network.
%
% A route model (src/model/routes/, see route_model) is the local system a
% route's synapses carry: a state per presynaptic neuron, driven by that
% neuron's rate, and an output factor F with synaptic output r .* F.
% SRNNCellTypePairs calls its handles and knows nothing else about it, so
% a model can be checked here by itself, against its own equations.
%
% TO TEST A NEW MODEL, add entries to model_cases() at the bottom: a
% configuration, how to draw a state in the model's range, and a
% `reference` that evaluates the model's equations for ONE neuron in plain
% scalar loops, written from the equations and not from the vectorised
% file. Every check of section 2 then runs on it.
%
% Checks:
%   1. The registry: route_model() lists the names; each resolves to a
%      struct with the interface's handles and its own name; an unknown
%      name errors route_model:UnknownModel.
%   2. For every case:
%      a. shapes: P.ns = P.n_b + P.n_g, P.npre, the per-neuron parameter
%         matrices; initial 1 x ns; rhs, linearize, steady_state; output on
%         an npre x ns x nt array equals output page by page
%      b. initial is the stated starting state and a rest state: the
%         derivative there is zero at r = 0
%      c. rhs equals the per-neuron reference loop (derivative and F)
%      d. linearize equals central differences of the model's OWN rhs and
%         output: A = d(ds)/ds, B = d(ds)/dr, Fs = dF/ds; L.F equals output
%      e. steady_state: the derivative vanishes there, and a long
%         integration from `initial` at constant rate arrives there
%      f. parameters are per neuron: replacing one neuron's rows of P by
%         another configuration's changes that neuron's derivative and
%         output to the other configuration's, and no other neuron's
%   3. route_product in particular: the geometric-mean combine has the
%      one-component steady state at every rate and the product its
%      square; the zero floor maps the steady state of r = 1 to F = 0 and
%      rest to F = 1; the route with no synapse dynamics has no state and
%      F = 1; a scalar tau_rel / tau_fac / G is repeated; the signature
%      tells routes apart; bad configurations error
%      <prefix>:InvalidSynapseConfig.
%
% Prints PASS/FAIL per check and a final banner. Assumes setup_paths.
%
% See also: route_model, route_product, test_jacobian_times,
%           test_SRNNCellTypePairs

fprintf('=== Testing the route models ===\n\n');
all_passed = true;
handles = {'compile', 'initial', 'rhs', 'output', 'linearize', 'steady_state'};

%% 1. the registry
names = route_model();
ok = iscellstr(names) && ~isempty(names); %#ok<ISCLSTR>
for k = 1:numel(names)
    model = route_model(names{k});
    ok = ok && isstruct(model) && strcmp(model.name, names{k}) && all(isfield(model, handles)) && ...
        all(cellfun(@(h) isa(model.(h), 'function_handle'), handles));
end
all_passed = check(sprintf('the registry resolves its %d name(s) to the interface', numel(names)), ok) && all_passed;
all_passed = check('an unknown model name errors route_model:UnknownModel', ...
    throws_id(@() route_model('no_such_model'), 'route_model:UnknownModel')) && all_passed;

%% 2. every case, through the interface only
cases = model_cases();
rates = [0; 0.02; 0.1; 0.25; 0.5; 1];        % the rate's range is [0, 1]
npre = numel(rates);
tol_rhs = 1e-13;                             % reference loop: same arithmetic, other order
tol_fd = 1e-7;                               % central differences, h = 1e-6
tol_ss = 1e-8;                               % integration to the steady state
rng(20261111, 'twister');
for c = 1:numel(cases)
    cs = cases(c);
    model = route_model(cs.model);
    opts = cs.opts;
    opts.npre = npre;
    P = model.compile(cs.std, cs.stf, opts);
    ns = P.ns;
    fprintf('  %s (ns = %d)\n', cs.label, ns);
    s = cs.draw(npre, ns);
    r = rates;

    % a. shapes
    s0 = model.initial(P);
    [ds, F] = model.rhs(s, r, P);
    L = model.linearize(s, r, P);
    ss = model.steady_state(r, P);
    nt = 3;
    s3 = cat(3, s, cs.draw(npre, ns), cs.draw(npre, ns));
    F3 = model.output(s3, P);
    by_page = true;
    for k = 1:nt
        by_page = by_page && isequal(size(F3), [npre 1 nt]) && isequal(F3(:, 1, k), model.output(s3(:, :, k), P));
    end
    ok = ns == P.n_b + P.n_g && P.npre == npre && strcmp(P.model, cs.model) && ...
        isequal(size(s0), [1 ns]) && isequal(size(ds), [npre ns]) && isequal(size(F), [npre 1]) && ...
        isequal(size(L.A, 1:3), [npre ns ns]) && isequal(size(L.B), [npre ns]) && ...
        isequal(size(L.F), [npre 1]) && isequal(size(L.Fs), [npre ns]) && ...
        isequal(size(ss.s), [npre ns]) && isequal(size(ss.F), [npre 1]) && ...
        isstruct(P.nominal) && isstruct(P.signature);
    all_passed = check('    shapes of P, initial, rhs, linearize, steady_state', ok) && all_passed;
    all_passed = check('    output on npre x ns x nt equals output page by page', by_page) && all_passed;

    % b. initial
    at_rest = model.rhs(repmat(s0, npre, 1), zeros(npre, 1), P);
    all_passed = check('    initial is the stated state, and at rest when r = 0', ...
        isequal(s0, cs.initial) && all(at_rest(:) == 0)) && all_passed;

    % c. rhs against the per-neuron loop
    ds_ref = zeros(npre, ns);
    F_ref = zeros(npre, 1);
    for i = 1:npre
        [ds_ref(i, :), F_ref(i)] = cs.reference(s(i, :), r(i), cs);
    end
    e = max([abs(ds(:) - ds_ref(:)); abs(F - F_ref); abs(model.output(s, P) - F_ref)]);
    all_passed = check(sprintf('    rhs and output equal the per-neuron reference loop (max %.1e)', e), ...
        e <= tol_rhs) && all_passed;

    % d. linearize against central differences of the model's own rhs and output
    h = 1e-6;
    A_fd = zeros(npre, ns, ns);
    Fs_fd = zeros(npre, ns);
    for k = 1:ns
        dk = zeros(1, ns);
        dk(k) = h;
        A_fd(:, :, k) = (model.rhs(s + dk, r, P) - model.rhs(s - dk, r, P)) / (2 * h);
        Fs_fd(:, k) = (model.output(s + dk, P) - model.output(s - dk, P)) / (2 * h);
    end
    B_fd = (model.rhs(s, r + h, P) - model.rhs(s, r - h, P)) / (2 * h);
    eA = rel_err(L.A, A_fd);
    eB = rel_err(L.B, B_fd);
    eF = rel_err(L.Fs, Fs_fd);
    all_passed = check(sprintf('    linearize equals central differences (A %.1e, B %.1e, Fs %.1e)', eA, eB, eF), ...
        max([eA, eB, eF]) <= tol_fd && isequal(L.F, F)) && all_passed;

    % e. steady state
    d_ss = model.rhs(ss.s, r, P);
    s_end = zeros(npre, ns);
    if ns > 0
        ode_opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-13);
        flow = @(~, v) reshape(model.rhs(reshape(v, npre, ns), r, P), [], 1);
        [~, v] = ode45(flow, [0 cs.T_settle], repmat(s0, npre, 1), ode_opts);
        s_end = reshape(v(end, :), npre, ns);
    end
    e_int = max(abs(s_end(:) - ss.s(:)));
    if isempty(e_int), e_int = 0; end
    all_passed = check(sprintf('    steady_state: derivative zero (max %.1e), reached by integration (max %.1e)', ...
        max(abs([d_ss(:); 0])), e_int), ...
        max(abs([d_ss(:); 0])) <= 1e-13 && e_int <= tol_ss && isequal(ss.F, model.output(ss.s, P))) && all_passed;

    % f. parameters are carried per neuron
    if ns > 0
        opts_alt = cs.opts;
        opts_alt.npre = npre;
        P_alt = model.compile(cs.std_alt, cs.stf_alt, opts_alt);
        Q = P;
        i_alt = 2;
        fields = fieldnames(P);
        n_moved = 0;
        for k = 1:numel(fields)
            v = P.(fields{k});
            if isnumeric(v) && size(v, 1) == npre && ~isscalar(v)
                v(i_alt, :) = P_alt.(fields{k})(i_alt, :);
                Q.(fields{k}) = v;
                n_moved = n_moved + 1;
            end
        end
        [ds_q, F_q] = model.rhs(s, r, Q);
        [ds_a, F_a] = model.rhs(s, r, P_alt);
        L_q = model.linearize(s, r, Q);
        L_a = model.linearize(s, r, P_alt);
        others = setdiff(1:npre, i_alt);
        ok = n_moved > 0 && isequal(ds_q(i_alt, :), ds_a(i_alt, :)) && isequal(F_q(i_alt), F_a(i_alt)) && ...
            isequal(ds_q(others, :), ds(others, :)) && isequal(F_q(others), F(others)) && ...
            ~isequal(ds_q(i_alt, :), ds(i_alt, :)) && ...
            isequal(L_q.A(i_alt, :, :), L_a.A(i_alt, :, :)) && isequal(L_q.B(i_alt, :), L_a.B(i_alt, :)) && ...
            isequal(L_q.Fs(i_alt, :), L_a.Fs(i_alt, :)) && isequal(L_q.A(others, :, :), L.A(others, :, :));
        all_passed = check(sprintf('    one neuron given another configuration''s rows (%d matrices): only it changes', ...
            n_moved), ok) && all_passed;
    end
end

%% 3. route_product in particular
fprintf('  route_product\n');
product = route_model('product');
r = rates;
% One ratio rho = tau_rel / tau_rec = 0.125 on both components.
rho = 0.125;
one  = product.compile(struct('tau_rec', 2, 'tau_rel', 0.25), [], struct('npre', npre));
two  = product.compile(struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5]), [], struct('npre', npre));
geo  = product.compile(struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5], 'combine', 'geomean'), [], struct('npre', npre));
F_one = product.steady_state(r, one).F;
e_geo = max(abs(product.steady_state(r, geo).F - F_one));
e_sq  = max(abs(product.steady_state(r, two).F - F_one .^ 2));
e_cf  = max(abs(F_one - 1 ./ (1 + r ./ rho)));
all_passed = check(sprintf(['    steady state: one component 1/(1 + r/rho) (%.1e); geomean of two equals it (%.1e); ' ...
    'product of two is its square (%.1e)'], e_cf, e_geo, e_sq), max([e_cf, e_geo, e_sq]) <= 1e-15) && all_passed;
geo1 = product.compile(struct('tau_rec', 2, 'tau_rel', 0.25, 'combine', 'geomean'), [], struct('npre', npre));
all_passed = check('    geomean on one component is the plain factor', ...
    ~geo1.geomean && isequal(product.steady_state(r, geo1).F, F_one)) && all_passed;

for combine = {'product', 'geomean'}
    fl = product.compile(struct('tau_rec', [2 5], 'tau_rel', [0.25 0.4], 'combine', combine{1}), ...
        struct('tau_dec', 1.5, 'tau_fac', 0.25, 'G', 3), struct('npre', 1, 'zero_floor', true));
    nf = product.compile(struct('tau_rec', [2 5], 'tau_rel', [0.25 0.4], 'combine', combine{1}), ...
        struct('tau_dec', 1.5, 'tau_fac', 0.25, 'G', 3), struct('npre', 1));
    top = product.steady_state(1, fl);
    D_top = product.output([top.s(1:2), 1], fl);            % facilitation taken out
    D_rest = product.output([1 1 1], fl);
    mid = [0.6 0.7 1];
    D_mid = product.output(mid, nf);
    ok = fl.zero_floor && ~nf.zero_floor && abs(D_top) <= 1e-15 && D_rest == 1 && ...
        abs(product.output(mid, fl) - (D_mid - fl.p_min) / (1 - fl.p_min)) <= 1e-15 && ...
        abs(fl.p_min - product.output([top.s(1:2), 1], nf)) <= 1e-15;
    all_passed = check(sprintf('    zero floor (%s): F = 0 at the steady state of r = 1, F = 1 at rest, linear between', ...
        combine{1}), ok) && all_passed;
end
st = product.compile([], struct('tau_dec', 1, 'tau_fac', 0.5, 'G', 2), struct('zero_floor', true));
all_passed = check('    zero floor on a route without depression is not in force', ...
    ~st.zero_floor && product.output(1.5, st) == 1.5) && all_passed;

empty = product.compile([], [], struct('npre', npre));
empty2 = product.compile(struct(), struct(), struct('npre', npre));
[ds_e, F_e] = product.rhs(zeros(npre, 0), r, empty);
L_e = product.linearize(zeros(npre, 0), r, empty);
ok = empty.ns == 0 && isequal(empty, empty2) && isequal(size(ds_e), [npre 0]) && isequal(F_e, ones(npre, 1)) && ...
    isequal(size(L_e.A, 1:3), [npre 0 0]) && isequal(size(L_e.B), [npre 0]) && isequal(L_e.F, ones(npre, 1)) && ...
    isequal(size(L_e.Fs), [npre 0]) && isequal(size(product.initial(empty)), [1 0]) && ...
    isequal(product.output(zeros(npre, 0, 4), empty), ones(npre, 1, 4));
all_passed = check('    the route with no synapse dynamics: no state, F = 1', ok) && all_passed;

sc = product.compile(struct('tau_rec', [2 4 8], 'tau_rel', 0.25), struct('tau_dec', [1 2], 'tau_fac', 0.5, 'G', 2));
ok = isequal(sc.nominal.tau_rel, [0.25 0.25 0.25]) && isequal(sc.nominal.tau_fac, [0.5 0.5]) && ...
    isequal(sc.nominal.G, [2 2]) && sc.n_b == 3 && sc.n_g == 2 && sc.npre == 1 && ...
    isequal(sc.tau_rest, [2 4 8 1 2]) && isequal(sc.tau_drive, [0.25 0.25 0.25 0.5 0.5]) && ...
    isequal(sc.target, [0 0 0 2 2]);
all_passed = check('    a scalar tau_rel, tau_fac or G is repeated; the stacked rows are [depression, facilitation]', ok) && all_passed;

sig = @(std_cfg, stf_cfg, varargin) product.compile(std_cfg, stf_cfg, varargin{:}).signature;
base = struct('tau_rec', [2 4], 'tau_rel', 0.25);
ok = isequal(sig(base, []), sig(base, struct())) && ...
    isequal(sig(base, [], struct('npre', 7)), sig(base, [])) && ...
    ~isequal(sig(base, []), sig(struct('tau_rec', [2 4], 'tau_rel', 0.3), [])) && ...
    ~isequal(sig(base, []), sig(struct('tau_rec', [2 4], 'tau_rel', 0.25, 'combine', 'geomean'), [])) && ...
    ~isequal(sig(base, []), sig(base, struct('tau_dec', 1, 'tau_fac', 1, 'G', 2))) && ...
    ~isequal(sig(base, []), sig(base, [], struct('zero_floor', true))) && ...
    ~isequal(sig(base, []), sig([], []));
all_passed = check('    the signature tells routes apart and ignores npre', ok) && all_passed;

bad = { ...
    {struct('tau_rec', 2), []}, ...                                      % tau_rel missing
    {struct('tau_rec', 2, 'tau_rel', 1, 'tau', 3), []}, ...              % unknown field
    {struct('tau_rec', [2 4], 'tau_rel', [1 2 3]), []}, ...              % wrong count
    {struct('tau_rec', -2, 'tau_rel', 1), []}, ...                       % not positive
    {struct('tau_rec', zeros(1, 0), 'tau_rel', 1), []}, ...              % no component
    {struct('tau_rec', 2, 'tau_rel', 1, 'combine', 'mean'), []}, ...     % bad combine
    {struct('tau_rec', 2, 'tau_rel', 1, 'model', 'product'), []}, ...    % not a field of this model
    {7, []}, ...                                                         % not a struct
    {[], struct('tau_dec', 1, 'tau_fac', 1)}, ...                        % G missing
    {[], struct('tau_dec', 1, 'tau_fac', 1, 'G', 0.5)}, ...              % G below one
    {[], struct('tau_dec', [1 2], 'tau_fac', 1, 'G', [2 3 4])}};         % wrong count
ok = true;
for k = 1:numel(bad)
    ok = ok && throws_id(@() product.compile(bad{k}{1}, bad{k}{2}), 'route_product:InvalidSynapseConfig') && ...
        throws_id(@() product.compile(bad{k}{1}, bad{k}{2}, struct('err_id_prefix', 'SRNNCellTypePairs')), ...
        'SRNNCellTypePairs:InvalidSynapseConfig');
end
all_passed = check(sprintf('    %d bad configurations error <prefix>:InvalidSynapseConfig', numel(bad)), ok) && all_passed;

fprintf('\n');
if all_passed
    fprintf('=== ALL route model TESTS PASSED ===\n');
else
    fprintf('=== SOME route model TESTS FAILED ===\n');
end

%% The cases. One entry per model and configuration.
function cases = model_cases()
% Fields:
%   label, model      a name for the printout; the registry name
%   std, stf, opts    the route's configuration and compile options
%   std_alt, stf_alt  another configuration of the same shape (check 2f)
%   initial           the starting state the model must report
%   draw(npre, ns)    a random state in the model's range
%   reference(s, r, cs)  [ds, F] for ONE neuron (s 1 x ns, r scalar),
%                     from the equations, in scalar loops
%   T_settle          how long to integrate to reach the steady state (s)
none = struct([]);
cases = struct('label', {}, 'model', {}, 'std', {}, 'stf', {}, 'opts', {}, 'std_alt', {}, 'stf_alt', {}, ...
    'initial', {}, 'draw', {}, 'reference', {}, 'T_settle', {});

std1 = struct('tau_rec', 2, 'tau_rel', 0.25);
std1_alt = struct('tau_rec', 3, 'tau_rel', 0.4);
std2 = struct('tau_rec', [2 5], 'tau_rel', [0.25 0.4]);
std2_alt = struct('tau_rec', [1.5 7], 'tau_rel', [0.3 0.2]);
std2_geo = struct('tau_rec', [2 5], 'tau_rel', [0.25 0.4], 'combine', 'geomean');
std2_geo_alt = struct('tau_rec', [1.5 7], 'tau_rel', [0.3 0.2], 'combine', 'geomean');
stf1 = struct('tau_dec', 1.5, 'tau_fac', 0.25, 'G', 3);
stf1_alt = struct('tau_dec', 1, 'tau_fac', 0.5, 'G', 2);
stf2 = struct('tau_dec', [1.5 4], 'tau_fac', [0.25 1], 'G', [3 1.5]);
stf2_alt = struct('tau_dec', [1 3], 'tau_fac', [0.5 2], 'G', [2 4]);
floor_on = struct('zero_floor', true);

% The product model's range: b in (0, 1], g in [1, G]. G is at most 4 here.
draw_product = @(n_b, n_g) @(npre, ns) [0.15 + 0.85 * rand(npre, n_b), 1 + 3 * rand(npre, n_g)];
add = @(cases, label, std_cfg, stf_cfg, opts, std_alt, stf_alt, n_b, n_g, T) [cases, struct( ...
    'label', ['product: ' label], 'model', 'product', 'std', std_cfg, 'stf', stf_cfg, 'opts', opts, ...
    'std_alt', std_alt, 'stf_alt', stf_alt, 'initial', ones(1, n_b + n_g), ...
    'draw', draw_product(n_b, n_g), 'reference', @product_reference, 'T_settle', T)];

cases = add(cases, 'no synapse dynamics', none, none, struct(), none, none, 0, 0, 1);
cases = add(cases, 'one depression component', std1, none, struct(), std1_alt, none, 1, 0, 80);
cases = add(cases, 'two depression components', std2, none, struct(), std2_alt, none, 2, 0, 200);
cases = add(cases, 'facilitation only', none, stf1, struct(), none, stf1_alt, 0, 1, 60);
cases = add(cases, 'one depression, one facilitation', std1, stf1, struct(), std1_alt, stf1_alt, 1, 1, 80);
cases = add(cases, 'two depression, two facilitation', std2, stf2, struct(), std2_alt, stf2_alt, 2, 2, 200);
cases = add(cases, 'geometric mean of two', std2_geo, none, struct(), std2_geo_alt, none, 2, 0, 200);
cases = add(cases, 'two depression, zero floor', std2, none, floor_on, std2_alt, none, 2, 0, 200);
cases = add(cases, 'geometric mean, zero floor, two facilitation', std2_geo, stf2, floor_on, std2_geo_alt, stf2_alt, 2, 2, 200);
end

function [ds, F] = product_reference(s, r, cs)
% Product depression times output-only facilitation for ONE neuron, from
% the equations (docs/EquationsParametersDocs/Equations_stability_paper.md):
%   db_m/dt = (1 - b_m)/tau_rec_m - b_m r/tau_rel_m
%   dg_l/dt = (1 - g_l)/tau_dec_l + (G_l - g_l) r/tau_fac_l
%   F = floor(combine(prod b)) * prod g
M = 0;
L = 0;
if ~isempty(cs.std), M = numel(cs.std.tau_rec); end
if ~isempty(cs.stf), L = numel(cs.stf.tau_dec); end
ds = zeros(1, M + L);
D = 1;
p_min = 1;
for m = 1:M
    tau_rec = cs.std.tau_rec(m);
    tau_rel = cs.std.tau_rel(min(m, numel(cs.std.tau_rel)));
    b = s(m);
    ds(m) = (1 - b) / tau_rec - b * r / tau_rel;
    D = D * b;
    p_min = p_min * tau_rel / (tau_rec + tau_rel);       % b_m at the steady state of r = 1
end
if M > 1 && isfield(cs.std, 'combine') && strcmp(cs.std.combine, 'geomean')
    D = nthroot(D, M);
    p_min = nthroot(p_min, M);
end
if M > 0 && isfield(cs.opts, 'zero_floor') && cs.opts.zero_floor
    D = (D - p_min) / (1 - p_min);
end
F = D;
for l = 1:L
    tau_dec = cs.stf.tau_dec(l);
    tau_fac = cs.stf.tau_fac(min(l, numel(cs.stf.tau_fac)));
    G = cs.stf.G(min(l, numel(cs.stf.G)));
    g = s(M + l);
    ds(M + l) = (1 - g) / tau_dec + (G - g) * r / tau_fac;
    F = F * g;
end
end

function e = rel_err(A, B)
% The largest difference, relative to the largest magnitude in B (absolute
% where B is empty or zero).
scale = max(abs(B(:)));
if isempty(scale) || scale == 0, scale = 1; end
e = max(abs(A(:) - B(:))) / scale;
if isempty(e), e = 0; end
end

function tf = throws_id(fn, id)
tf = false;
try
    fn();
catch err
    tf = strcmp(err.identifier, id);
end
end

function ok = check(label, ok)
if ok
    fprintf('  PASS  %s\n', label);
else
    fprintf('  FAIL  %s\n', label);
end
end
