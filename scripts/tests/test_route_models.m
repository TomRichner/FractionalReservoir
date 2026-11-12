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
%         output: A = d(ds)/ds, B = d(ds)/dr, Fs = dF/ds; its F equals output
%      e. steady_state: the derivative vanishes there, and a long
%         integration from `initial` at constant rate arrives there
%      f. parameters are per neuron: replacing one neuron's rows of P by
%         another configuration's changes that neuron's derivative and
%         output to the other configuration's, and no other neuron's
%      g. the per-neuron spread (opts.std_sigma, std_z, stf_sigma, stf_z):
%         with parameters that differ from neuron to neuron, rhs, output
%         and steady_state are, neuron by neuron, those of a route
%         configured with that neuron's values (jittered_config, written
%         from route_model's rule and column table), and linearize equals
%         central differences; a zero spread compiles to exactly the route
%         without one
%   3. route_product in particular: the geometric-mean combine has the
%      one-component steady state at every rate and the product its
%      square; the zero floor maps the steady state of r = 1 to F = 0 and
%      rest to F = 1; the route with no synapse dynamics has no state and
%      F = 1; a scalar tau_rel / tau_fac / G is repeated; the signature
%      tells routes apart; bad configurations error
%      <prefix>:InvalidSynapseConfig.
%   3b. route_pool in particular: the steady state against the roots of
%      its quadratic and, at a very large reserve, against the single
%      pool; it is accurate where the reserve is nearly empty; the
%      sustained output rises to its ceiling; the pools' derivatives point
%      inward on the boundary of [0, 1]^2; the 2 x 2 is full and not
%      symmetric; the signature tells n2 apart and a chain from a product
%      with the same time constants; each refusal has its own identifier.
%   3c. route_pool with release-coupled facilitation: every variant
%      carries exactly its variables and states the coupling it has; the
%      steady state is g = (1 + G k r)/(1 + k r), h = r and the
%      facilitation-free pools at the effective rate g r, so the sustained
%      output is that of the route without facilitation at g r and its
%      ceiling does not move (Plan A section 4.6's numbers); G = 1 is the
%      route without facilitation, exactly; g and h cannot leave [1, G]
%      and [0, 1]; the single pool without facilitation is route_product's
%      one component; the display factors are b1 and g (not g h); the
%      signature tells every variant apart; each refusal has its own
%      identifier.
%   3d. the per-neuron spread across the models: tau_rec_1 and tau_rel_1
%      are the same per neuron on a one-component route, a product, a
%      single pool and a chain, and tau_rec_2 comes from one standard
%      normal on the product and the chain; every quantity reads its
%      documented column; G = 1 + (G - 1) exp(sigma z) >= 1; values are
%      positive with the nominal as median and sigma as log-spread; the
%      chain's ordering rule raises tau_rec_2 alone, to exactly the floor,
%      and acts nowhere else; the zero floor's p_min and the two display
%      handles (both models have `facilitation`) go by row; a spread
%      without usable draws errors <prefix>:SpreadDraws.
%   4.SRNNCellTypePairs.jacobian_times applies a route's A in two ways (its
%      diagonal where the model states A is diagonal, the full per-neuron
%      product otherwise): on a network with every kind of route, forcing
%      the full way changes nothing and still equals J * Y. The one check
%      here that builds a network. The product model never takes the full
%      way; the pool model always does, and test_std_chain_model checks
%      that on a network of its own.
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
    L = linearized(model, s, r, P);
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
        isstruct(P.nominal) && isstruct(P.signature) && ...
        islogical(P.coupling) && isequal(size(P.coupling), [ns ns]);
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
    % P.coupling may overstate what depends on what, never understate: where
    % it says false, the derivative does not move at all (the differences
    % are exactly zero) and A is exactly zero.
    all_passed = check('    P.coupling does not understate: A is zero wherever it is false', ...
        all(A_fd(:, ~P.coupling) == 0, 'all') && all(L.A(:, ~P.coupling) == 0, 'all')) && all_passed;

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
        L_q = linearized(model, s, r, Q);
        L_a = linearized(model, s, r, P_alt);
        others = setdiff(1:npre, i_alt);
        ok = n_moved > 0 && isequal(ds_q(i_alt, :), ds_a(i_alt, :)) && isequal(F_q(i_alt), F_a(i_alt)) && ...
            isequal(ds_q(others, :), ds(others, :)) && isequal(F_q(others), F(others)) && ...
            ~isequal(ds_q(i_alt, :), ds(i_alt, :)) && ...
            isequal(L_q.A(i_alt, :, :), L_a.A(i_alt, :, :)) && isequal(L_q.B(i_alt, :), L_a.B(i_alt, :)) && ...
            isequal(L_q.Fs(i_alt, :), L_a.Fs(i_alt, :)) && isequal(L_q.A(others, :, :), L.A(others, :, :));
        all_passed = check(sprintf('    one neuron given another configuration''s rows (%d matrices): only it changes', ...
            n_moved), ok) && all_passed;
    end

    % g. the per-neuron spread: every neuron its own parameters. The route
    % is compiled with a spread on both mechanisms and fixed standard
    % normals; neuron i must then behave as a route configured with the
    % values the documented rule gives it (jittered_config below, written
    % from route_model's column table and not from the models), through
    % the same per-neuron reference loop as check c.
    if ns > 0
        sig_std = 0.3;
        sig_stf = 0.2;
        z_std = randn(npre, 9);
        z_stf = randn(npre, 10);
        opts_j = cs.opts;
        opts_j.npre = npre;
        opts_j.std_sigma = sig_std; opts_j.std_z = z_std;
        opts_j.stf_sigma = sig_stf; opts_j.stf_z = z_stf;
        Pj = model.compile(cs.std, cs.stf, opts_j);
        [ds_j, F_j] = model.rhs(s, r, Pj);
        ds_ref = zeros(npre, ns);
        F_ref = zeros(npre, 1);
        ss_ref = zeros(npre, ns);
        for i = 1:npre
            cs_i = jittered_config(cs, sig_std, z_std(i, :), sig_stf, z_stf(i, :));
            [ds_ref(i, :), F_ref(i)] = cs.reference(s(i, :), r(i), cs_i);
            % the same neuron as a route of its own, at its nominal values
            opts_i = cs.opts;
            P_i = model.compile(cs_i.std, cs_i.stf, opts_i);
            ss_i = model.steady_state(r(i), P_i);
            ss_ref(i, :) = ss_i.s;
        end
        e = max([abs(ds_j(:) - ds_ref(:)); abs(F_j - F_ref); abs(model.output(s, Pj) - F_ref)]);
        L_j = linearized(model, s, r, Pj);
        A_fd = zeros(npre, ns, ns);
        Fs_fd = zeros(npre, ns);
        for k = 1:ns
            dk = zeros(1, ns);
            dk(k) = h;
            A_fd(:, :, k) = (model.rhs(s + dk, r, Pj) - model.rhs(s - dk, r, Pj)) / (2 * h);
            Fs_fd(:, k) = (model.output(s + dk, Pj) - model.output(s - dk, Pj)) / (2 * h);
        end
        B_fd = (model.rhs(s, r + h, Pj) - model.rhs(s, r - h, Pj)) / (2 * h);
        e_lin = max([rel_err(L_j.A, A_fd), rel_err(L_j.B, B_fd), rel_err(L_j.Fs, Fs_fd)]);
        ss_j = model.steady_state(r, Pj);
        d_ss = model.rhs(ss_j.s, r, Pj);
        e_ss = max(abs(ss_j.s(:) - ss_ref(:)));
        differ = ~isequal(ds_j(2:end, :), ds(2:end, :));       % (row 1 has r = 0)
        all_passed = check(sprintf(['    with a spread: rhs and output equal the reference loop on each neuron''s own ' ...
            'values (max %.1e); linearize equals central differences (%.1e)'], e, e_lin), ...
            Pj.jittered && ~P.jittered && differ && e <= tol_rhs && e_lin <= tol_fd && isequal(L_j.F, F_j)) && all_passed;
        all_passed = check(sprintf(['    with a spread: steady_state is each neuron''s own (max %.1e) and the ' ...
            'derivative vanishes there (max %.1e)'], e_ss, max(abs(d_ss(:)))), ...
            e_ss <= 1e-13 && max(abs(d_ss(:))) <= 1e-13 && isequal(ss_j.F, model.output(ss_j.s, Pj))) && all_passed;
        % A zero spread is the route without one, to the last bit, with
        % or without draws handed in; and the two mechanisms are separate.
        opts_0 = opts_j;
        opts_0.std_sigma = 0; opts_0.stf_sigma = 0;
        P0 = model.compile(cs.std, cs.stf, opts_0);
        opts_s = opts_j; opts_s.stf_sigma = 0;
        opts_f = opts_j; opts_f.std_sigma = 0;
        Ps = model.compile(cs.std, cs.stf, opts_s);
        Pf = model.compile(cs.std, cs.stf, opts_f);
        only_std = Ps.jittered == (P.n_b > 0) && Pf.jittered == (P.n_g > 0);
        all_passed = check('    a zero spread compiles to the route without one, exactly; each spread is in force only where the route has the mechanism', ...
            isequaln(P0, P) && only_std) && all_passed;
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
L_e = linearized(product, zeros(npre, 0), r, empty);
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

%% 3b. route_pool in particular
fprintf('  route_pool\n');
pool = route_model('pool');
chain_cfg = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
r = [rates; 3; 30];                          % beyond the rate's range too: the formula is general
P = pool.compile(chain_cfg, [], struct('npre', numel(r)));
ss = pool.steady_state(r, P);
% (C4) from the quadratic k2 b1^2 - (1 + k1 + k2) b1 + 1 = 0, by roots().
e_root = 0;
for i = 1:numel(r)
    k1 = r(i) * 2 / 0.25;
    k2 = r(i) * 5 / (2 * 0.25);
    if k2 == 0
        b1 = 1 / (1 + k1);
    else
        b1 = min(roots([k2, -(1 + k1 + k2), 1]));        % the root in [0, 1]
    end
    e_root = max([e_root, abs(ss.s(i, 1) - b1), abs(ss.s(i, 2) - (1 - k2 * b1))]);
end
all_passed = check(sprintf('    steady state equals the root of the quadratic in [0, 1], and b2 = 1 - k2 b1 (max %.1e)', e_root), ...
    e_root <= 1e-14 && all(ss.s(:) >= 0 & ss.s(:) <= 1)) && all_passed;
% The values Plan A section 1 quotes for these parameters.
quoted = pool.steady_state([0.1; 0.25; 0.5; 1], pool.compile(chain_cfg, [], struct('npre', 4))).s(:, 1);
all_passed = check('    b1 at r = 0.1, 0.25, 0.5, 1 is 0.420, 0.200, 0.106, 0.054', ...
    isequal(round(quoted, 3), [0.420; 0.200; 0.106; 0.054])) && all_passed;

% A very large reserve is the single pool.
big = chain_cfg;
big.n2 = 1e12;
ss_big = pool.steady_state(r, pool.compile(big, [], struct('npre', numel(r))));
single = 1 ./ (1 + r * 2 / 0.25);
e_single = max(abs(ss_big.s(:, 1) - single) ./ single);
all_passed = check(sprintf('    n2 = 1e12: b1 is the single pool''s 1/(1 + k1) (rel %.1e), b2 within 1e-9 of 1', e_single), ...
    e_single <= 1e-10 && all(abs(ss_big.s(:, 2) - 1) <= 1e-9)) && all_passed;
% A reserve that is nearly empty (k2 large against 1 + k1): b2 is small,
% and still satisfies (C1) = 0 to round-off relative to its own size.
slow = struct('model', 'chain', 'tau_rec', [0.05 400], 'tau_rel', 0.25, 'n2', 0.5);
r_slow = [0.25; 1; 30];
ss_slow = pool.steady_state(r_slow, pool.compile(slow, [], struct('npre', 3)));
k1 = r_slow * 0.05 / 0.25;
k2 = r_slow * 400 / (0.5 * 0.25);
e_b2 = max(abs(ss_slow.s(:, 2) .* (1 - ss_slow.s(:, 1)) - k1 .* ss_slow.s(:, 1)) ./ (k1 .* ss_slow.s(:, 1)));
all_passed = check(sprintf('    reserve nearly empty (b2 down to %.1e): b2 (1 - b1) = k1 b1 holds (rel %.1e)', ...
    min(ss_slow.s(:, 2)), e_b2), all(k2 > 1 + k1) && e_b2 <= 1e-13 && all(ss_slow.s(:, 2) > 0)) && all_passed;

% The sustained output r b1 rises with the rate to
% 1/(tau_rec1/tau_rel + tau_rec2/(n2 tau_rel)) = 1/18.
r_up = [0.01 0.05 0.1 0.25 0.5 1 10 1e3 1e6]';
out = r_up .* pool.steady_state(r_up, pool.compile(chain_cfg, [], struct('npre', numel(r_up)))).F;
all_passed = check(sprintf('    sustained output rises with the rate to its ceiling 1/18 (at r = 1e6: 1/%.4f)', 1 / out(end)), ...
    all(diff(out) > 0) && out(end) < 1 / 18 && abs(out(end) - 1 / 18) <= 1e-6) && all_passed;

% The derivative points inward on the boundary of [0, 1]^2.
edge = rand(50, 1);
rr = rand(50, 1);
P50 = pool.compile(chain_cfg, [], struct('npre', 50));
d_b1_0 = pool.rhs([zeros(50, 1), edge], rr, P50);
d_b1_1 = pool.rhs([ones(50, 1), edge], rr, P50);
d_b2_0 = pool.rhs([edge, zeros(50, 1)], rr, P50);
d_b2_1 = pool.rhs([edge, ones(50, 1)], rr, P50);
all_passed = check('    the pools cannot leave [0, 1]: the derivative points inward on each edge', ...
    all(d_b1_0(:, 1) >= 0) && all(d_b1_1(:, 1) <= 0) && all(d_b2_0(:, 2) > 0) && all(d_b2_1(:, 2) <= 0)) && all_passed;

% The 2 x 2 is full and not symmetric, and the model says so.
[A, B, ~, Fs] = pool.linearize([0.4 0.7], 0.3, pool.compile(chain_cfg, [], struct()));
A = reshape(A, 2, 2);
all_passed = check('    A is a full 2 x 2, not symmetric; P.coupling is full; the rate and the output touch b1 only', ...
    all(A(:) ~= 0) && A(1, 2) ~= A(2, 1) && isequal(P.coupling, true(2)) && B(2) == 0 && B(1) ~= 0 && ...
    isequal(Fs, [1 0])) && all_passed;
all_passed = check('    depression(b, P) is the releasable pool, with the shape of prod(b, 2)', ...
    isequal(pool.depression(cat(3, [0.4 0.7; 0.2 0.9], [0.1 0.3; 0.5 0.6]), P), ...
    cat(3, [0.4; 0.2], [0.1; 0.5])) && isequal(pool.initial(P), [1 1]) && P.n_b == 2 && P.n_g == 0) && all_passed;

sig = @(cfg, varargin) pool.compile(cfg, [], varargin{:}).signature;
other_n2 = chain_cfg;
other_n2.n2 = 3;
other_tau = chain_cfg;
other_tau.tau_rec = [2 6];
as_product = product.compile(struct('tau_rec', [2 5], 'tau_rel', 0.25), []).signature;
ok = isequal(sig(chain_cfg), sig(chain_cfg, struct('npre', 7))) && ~isequal(sig(chain_cfg), sig(other_n2)) && ...
    ~isequal(sig(chain_cfg), sig(other_tau)) && ~isequal(sig(chain_cfg), as_product);
all_passed = check('    the signature tells n2 and the time constants apart, a chain from a product, and ignores npre', ok) && all_passed;

without = @(name) rmfield(chain_cfg, name);
with = @(name, value) setfield(chain_cfg, name, value);
stf_cfg = struct('tau_dec', 1.5, 'tau_fac', 0.25, 'G', 3);
refusals = { ...
    {with('tau_rec', 2), [], struct(), 'ChainTauRecCount'}, ...
    {with('tau_rec', [2 5 9]), [], struct(), 'ChainTauRecCount'}, ...
    {with('tau_rel', [0.25 0.4]), [], struct(), 'ChainTauRelNotScalar'}, ...
    {without('n2'), [], struct(), 'ChainMissingN2'}, ...
    {with('n2', 0), [], struct(), 'ChainInvalidN2'}, ...
    {with('n2', Inf), [], struct(), 'ChainInvalidN2'}, ...
    {with('n2', [2 3]), [], struct(), 'ChainInvalidN2'}, ...
    {with('combine', 'geomean'), [], struct(), 'ChainCombine'}, ...
    {with('combine', 'product'), [], struct(), 'ChainCombine'}, ...
    {chain_cfg, [], struct('zero_floor', true), 'ChainZeroFloor'}, ...
    {chain_cfg, stf_cfg, struct(), 'ChainStfOutputOnly'}, ...
    {with('tau_rec', [2 -5]), [], struct(), 'InvalidSynapseConfig'}, ...
    {with('tau_rel', 0), [], struct(), 'InvalidSynapseConfig'}, ...
    {without('tau_rel'), [], struct(), 'InvalidSynapseConfig'}, ...
    {with('tau', 3), [], struct(), 'InvalidSynapseConfig'}, ...
    {with('model', 'product'), [], struct(), 'InvalidSynapseConfig'}, ...
    {[], [], struct(), 'InvalidSynapseConfig'}};
ok = true;
for k = 1:numel(refusals)
    bad_k = refusals{k};
    as_class = bad_k{3};
    as_class.err_id_prefix = 'SRNNCellTypePairs';
    ok_k = throws_id(@() pool.compile(bad_k{1}, bad_k{2}, bad_k{3}), ['route_pool:' bad_k{4}]) && ...
        throws_id(@() pool.compile(bad_k{1}, bad_k{2}, as_class), ['SRNNCellTypePairs:' bad_k{4}]);
    if ~ok_k, fprintf('      refusal %d (%s) did not raise its identifier\n', k, bad_k{4}); end
    ok = ok && ok_k;
end
all_passed = check(sprintf('    %d bad chain configurations each raise their own <prefix>:<identifier>', numel(refusals)), ok) && all_passed;

%% 3c. route_pool with release-coupled facilitation
fprintf('  route_pool, release-coupled facilitation\n');
% Tom's values (Plan A, section 4.6): G = 5, tau_fac = tau_rel = 0.25,
% tau_dec = 1.5, tau_h = 0.5.
fac = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5);
fac_h = fac;
fac_h.tau_h = 0.5;
one_pool = struct('tau_rec', 2, 'tau_rel', 0.25);
variants = { ...          std        stf     n_b n_g  coupling over the variables present
    'single pool',        one_pool,  [],     1, 0, true; ...
    'chain',              chain_cfg, [],     2, 0, true(2); ...
    'g',                  [],        fac,    0, 1, true; ...
    'g, h',               [],        fac_h,  0, 2, [true true; false true]; ...
    'single pool, g',     one_pool,  fac,    1, 1, [true true; false true]; ...
    'single pool, g, h',  one_pool,  fac_h,  1, 2, [true true false; false true true; false false true]; ...
    'chain, g',           chain_cfg, fac,    2, 1, [true true true; true true false; false false true]; ...
    'chain, g, h',        chain_cfg, fac_h,  2, 2, [true true true false; true true false false; ...
                                                    false false true true; false false false true]};
ok = true;
sigs = cell(1, size(variants, 1));
for k = 1:size(variants, 1)
    Pk = pool.compile(variants{k, 2}, variants{k, 3}, struct('npre', 3));
    sigs{k} = Pk.signature;
    s0 = [ones(1, variants{k, 4}), ones(1, min(variants{k, 5}, 1)), zeros(1, variants{k, 5} == 2)];
    ok_k = Pk.n_b == variants{k, 4} && Pk.n_g == variants{k, 5} && isequal(Pk.coupling, variants{k, 6}) && ...
        isequal(pool.initial(Pk), s0) && Pk.has_g == (variants{k, 5} >= 1) && Pk.has_h == (variants{k, 5} == 2) && ...
        isfield(Pk, 'tau_rec1') == (variants{k, 4} >= 1) && isfield(Pk, 'tau_rec2') == (variants{k, 4} == 2) && ...
        isfield(Pk, 'n2') == (variants{k, 4} == 2) && isfield(Pk, 'tau_dec') == (variants{k, 5} >= 1) && ...
        isfield(Pk, 'tau_h') == (variants{k, 5} == 2);
    if ~ok_k, fprintf('      variant ''%s'' is not as stated\n', variants{k, 1}); end
    ok = ok && ok_k;
end
all_passed = check(sprintf(['    %d variants: each carries exactly its variables (columns [b1, b2 | g, h]), ' ...
    'its parameters and its coupling; pools and g start at 1, h at 0'], size(variants, 1)), ok) && all_passed;
distinct = true;
for k = 1:numel(sigs)
    for j = k + 1:numel(sigs)
        distinct = distinct && ~isequal(sigs{k}, sigs{j});
    end
end
sig2 = @(std_cfg, stf_cfg) pool.compile(std_cfg, stf_cfg, struct()).signature;
other_h = fac_h;
other_h.tau_h = 0.6;
other_G = fac;
other_G.G = 4;
output_only = product.compile(one_pool, struct('tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5)).signature;
ok = distinct && ~isequal(sig2(one_pool, fac_h), sig2(one_pool, other_h)) && ...
    ~isequal(sig2(one_pool, fac), sig2(one_pool, other_G)) && ...
    ~isequal(sig2(one_pool, fac), output_only) && ...
    isequal(sig2(chain_cfg, fac_h), pool.compile(chain_cfg, fac_h, struct('npre', 9)).signature) && ...
    ~isfield(sig2(chain_cfg, []), 'tau_h') && ~isfield(sig2(chain_cfg, []), 'stf_model');
all_passed = check(['    the signature tells the variants apart, tau_h and G apart, and a coupled route from an ' ...
    'output-only one with the same numbers; a route without facilitation states no facilitation field'], ok) && all_passed;

% The steady state. kappa = tau_dec/tau_fac = 6.
r4 = [0.02; 0.1; 0.25; 1];
g_cf = (1 + 5 * 6 * r4) ./ (1 + 6 * r4);
e_g = 0;
e_out = 0;
e_h = 0;
rel_out = zeros(4, 2);
for k = 5:8                                  % the four variants with depression and facilitation
    Pk = pool.compile(variants{k, 2}, variants{k, 3}, struct('npre', 4));
    P_plain = pool.compile(variants{k, 2}, [], struct('npre', 4));
    ss_k = pool.steady_state(r4, Pk);
    n_b = variants{k, 4};
    e_g = max(e_g, max(abs(ss_k.s(:, n_b + 1) - g_cf)));
    if variants{k, 5} == 2, e_h = max(e_h, max(abs(ss_k.s(:, n_b + 2) - r4))); end
    % r g b1 is the facilitation-free route's sustained output at the rate g r.
    at_gr = pool.steady_state(g_cf .* r4, P_plain);
    e_out = max(e_out, max(abs(r4 .* ss_k.F - g_cf .* r4 .* at_gr.F)));
    e_out = max(e_out, max(abs(ss_k.s(:, 1:n_b) - at_gr.s), [], 'all'));
    rel_out(:, n_b) = ss_k.F ./ pool.steady_state(r4, P_plain).F;
end
for k = 3:4                                  % facilitation alone: F = g
    ss_k = pool.steady_state(r4, pool.compile([], variants{k, 3}, struct('npre', 4)));
    e_g = max([e_g; abs(ss_k.s(:, 1) - g_cf); abs(ss_k.F - g_cf)]);
end
all_passed = check(sprintf(['    steady state: g = (1 + G k r)/(1 + k r) (max %.1e), h = r (max %.1e), pools and output ' ...
    'those of the route without facilitation at the rate g r (max %.1e)'], e_g, e_h, e_out), ...
    e_g <= 1e-15 && e_h == 0 && e_out <= 1e-15) && all_passed;
all_passed = check(['    Plan A 4.6''s numbers at r = 0.02, 0.1, 0.25, 1: g = 1.43, 2.50, 3.40, 4.43; sustained output ' ...
    'over no facilitation 1.35, 1.50, 1.31, 1.09 (single pool) and 1.32, 1.19, 1.08, 1.02 (chain)'], ...
    isequal(round(g_cf, 2), [1.43; 2.50; 3.40; 4.43]) && isequal(round(rel_out(:, 1), 2), [1.35; 1.50; 1.31; 1.09]) && ...
    isequal(round(rel_out(:, 2), 2), [1.32; 1.19; 1.08; 1.02])) && all_passed;
% Coupled facilitation does not raise the chain's ceiling 1/18.
out_f = r_up .* pool.steady_state(r_up, pool.compile(chain_cfg, fac_h, struct('npre', numel(r_up)))).F;
all_passed = check(sprintf(['    the sustained output with facilitation rises to the same ceiling 1/18 (at r = 1e6: ' ...
    '1/%.4f) and is above the output without it at every rate'], 1 / out_f(end)), ...
    all(diff(out_f) > 0) && out_f(end) < 1 / 18 && abs(out_f(end) - 1 / 18) <= 1e-6 && all(out_f > out)) && all_passed;

% G = 1: g has nowhere to go, and at g = 1 the pools' equations are the
% facilitation-free route's, to the last bit.
flat = fac_h;
flat.G = 1;
P_flat = pool.compile(chain_cfg, flat, struct('npre', 50));
s_flat = [rand(50, 2), ones(50, 1), rand(50, 1)];
[d_flat, F_flat] = pool.rhs(s_flat, rr, P_flat);
[d_plain, F_plain] = pool.rhs(s_flat(:, 1:2), rr, P50);
all_passed = check('    G = 1: dg/dt = 0 at g = 1, and the pools and the output are the chain''s without facilitation, exactly', ...
    all(d_flat(:, 3) == 0) && isequal(d_flat(:, 1:2), d_plain) && isequal(F_flat, F_plain)) && all_passed;

% g stays in [1, G] and h in [0, 1] (for a rate in [0, 1]); b1 stays in
% [0, 1] with g in its release term.
P_edge = pool.compile(chain_cfg, fac_h, struct('npre', 50));
P_edge1 = pool.compile(one_pool, fac, struct('npre', 50));
at = @(g, h) pool.rhs([edge, edge, g, h], rr, P_edge);
d_g1 = at(ones(50, 1), edge);
d_gG = at(5 * ones(50, 1), edge);
d_h0 = at(1 + 4 * edge, zeros(50, 1));
d_h1 = at(1 + 4 * edge, ones(50, 1));
d_g1r = pool.rhs([edge, ones(50, 1)], rr, P_edge1);
d_gGr = pool.rhs([edge, 5 * ones(50, 1)], rr, P_edge1);
d_b10 = pool.rhs([zeros(50, 1), edge, 1 + 4 * edge, edge], rr, P_edge);
d_b11 = pool.rhs([ones(50, 1), edge, 1 + 4 * edge, edge], rr, P_edge);
all_passed = check('    g cannot leave [1, G], h cannot leave [0, 1], b1 cannot leave [0, 1] with g in the release term', ...
    all(d_g1(:, 3) >= 0) && all(d_gG(:, 3) <= 0) && all(d_h0(:, 4) >= 0) && all(d_h1(:, 4) <= 0) && ...
    all(d_g1r(:, 2) >= 0) && all(d_gGr(:, 2) <= 0) && all(d_b10(:, 1) >= 0) && all(d_b11(:, 1) <= 0)) && all_passed;

% The single pool without facilitation is route_product's one component.
P_single = pool.compile(one_pool, [], struct('npre', npre));
P_prod = product.compile(one_pool, [], struct('npre', npre));
s1 = 0.05 + 0.9 * rand(npre, 1);
[d_s, F_s] = pool.rhs(s1, rates, P_single);
[d_p, F_p] = product.rhs(s1, rates, P_prod);
L_s = linearized(pool, s1, rates, P_single);
L_p = linearized(product, s1, rates, P_prod);
e_same = max([abs(d_s - d_p); abs(F_s - F_p); abs(L_s.A(:) - L_p.A(:)); abs(L_s.B - L_p.B); abs(L_s.Fs - L_p.Fs); ...
    abs(pool.steady_state(rates, P_single).s - product.steady_state(rates, P_prod).s)]);
all_passed = check(sprintf('    the single pool without facilitation is the product''s one component (max %.1e)', e_same), ...
    e_same <= 1e-15) && all_passed;

% The display factors: b1, and g (not g h).
P_all = pool.compile(chain_cfg, fac_h, struct('npre', 2));
g_arr = cat(3, [2 0.3; 4 0.1], [1.5 0.9; 3 0.5]);
all_passed = check('    facilitation(g, P) is g, the first g-block column, with the shape of prod(g, 2); never g h', ...
    isequal(pool.facilitation(g_arr, P_all), g_arr(:, 1, :)) && ...
    isequal(pool.depression(cat(3, [0.4 0.7; 0.2 0.9], [0.1 0.3; 0.5 0.6]), P_all), cat(3, [0.4; 0.2], [0.1; 0.5])) && ...
    isequal(pool.output([0.4 0.7 2 0.3; 0.2 0.9 4 0.1], P_all), [0.8; 0.8])) && all_passed;

two_pools = struct('tau_rec', [2 5], 'tau_rel', 0.25);
set_fac = @(name, value) setfield(fac, name, value);
set_std = @(cfg, name, value) setfield(cfg, name, value);
refusals = { ...
    {chain_cfg, rmfield(fac, 'model'), struct(), 'ChainStfOutputOnly'}, ...
    {chain_cfg, set_fac('model', 'product'), struct(), 'ChainStfOutputOnly'}, ...
    {[], rmfield(fac, 'model'), struct(), 'ChainStfOutputOnly'}, ...
    {two_pools, fac, struct(), 'ChainStfProductStd'}, ...
    {struct('tau_rec', [2 5 9], 'tau_rel', 0.25), fac_h, struct(), 'ChainStfProductStd'}, ...
    {set_std(one_pool, 'combine', 'geomean'), fac, struct(), 'ChainCombine'}, ...
    {set_std(chain_cfg, 'combine', 'product'), fac, struct(), 'ChainCombine'}, ...
    {one_pool, fac, struct('zero_floor', true), 'ChainZeroFloor'}, ...
    {chain_cfg, fac_h, struct('zero_floor', true), 'ChainZeroFloor'}, ...
    {one_pool, set_fac('tau_dec', [1.5 3]), struct(), 'ChainStfNotScalar'}, ...
    {one_pool, set_fac('tau_fac', [0.25 1]), struct(), 'ChainStfNotScalar'}, ...
    {chain_cfg, set_fac('G', [5 2]), struct(), 'ChainStfNotScalar'}, ...
    {one_pool, set_fac('G', 0.5), struct(), 'ChainStfInvalidG'}, ...
    {[], set_fac('tau_h', 0), struct(), 'ChainStfInvalidTauH'}, ...
    {[], set_fac('tau_h', -1), struct(), 'ChainStfInvalidTauH'}, ...
    {one_pool, set_fac('tau_h', Inf), struct(), 'ChainStfInvalidTauH'}, ...
    {one_pool, set_fac('tau_h', [0.5 1]), struct(), 'ChainStfInvalidTauH'}, ...
    {one_pool, set_fac('tau_h', []), struct(), 'ChainStfInvalidTauH'}, ...
    {one_pool, rmfield(fac, 'G'), struct(), 'InvalidSynapseConfig'}, ...
    {one_pool, set_fac('tau_dec', -1), struct(), 'InvalidSynapseConfig'}, ...
    {one_pool, set_fac('tau', 3), struct(), 'InvalidSynapseConfig'}, ...
    {set_std(one_pool, 'n2', 2), fac, struct(), 'InvalidSynapseConfig'}, ...
    {set_std(one_pool, 'tau_rel', [0.25 0.4]), fac, struct(), 'InvalidSynapseConfig'}, ...
    {one_pool, 7, struct(), 'InvalidSynapseConfig'}};
ok = true;
for k = 1:numel(refusals)
    bad_k = refusals{k};
    as_class = bad_k{3};
    as_class.err_id_prefix = 'SRNNCellTypePairs';
    ok_k = throws_id(@() pool.compile(bad_k{1}, bad_k{2}, bad_k{3}), ['route_pool:' bad_k{4}]) && ...
        throws_id(@() pool.compile(bad_k{1}, bad_k{2}, as_class), ['SRNNCellTypePairs:' bad_k{4}]);
    if ~ok_k, fprintf('      refusal %d (%s) did not raise its identifier\n', k, bad_k{4}); end
    ok = ok && ok_k;
end
floor_free = pool.compile([], fac_h, struct('zero_floor', true));
all_passed = check(sprintf(['    %d bad configurations with facilitation each raise their own <prefix>:<identifier>; ' ...
    'the zero floor on a route without depression is not in force and not refused'], numel(refusals)), ...
    ok && floor_free.n_b == 0 && floor_free.n_g == 2) && all_passed;

%% 3d. the per-neuron spread across the models
% What section 2g cannot see case by case: that a quantity has the SAME
% value whichever model the route carries; the chain's ordering rule; the
% ceiling; what a model refuses.
fprintf('  the per-neuron spread\n');
product = route_model('product');
pool = route_model('pool');
nn = 4000;
rng(77, 'twister');
zs = randn(nn, 7);
zf = randn(nn, 7);
sg = 0.25;
o = struct('npre', nn, 'std_sigma', sg, 'std_z', zs, 'stf_sigma', sg, 'stf_z', zf);
one = struct('tau_rec', 2, 'tau_rel', 0.25);
two = struct('tau_rec', [2 5], 'tau_rel', [0.25 0.4]);
three = struct('tau_rec', [2 5 9], 'tau_rel', [0.25 0.4 0.3]);
chain_cfg = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
fac_cfg = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5, 'tau_h', 0.5);
fac_out = struct('tau_dec', [1.5 4], 'tau_fac', [0.25 1], 'G', [5 1.5]);
P1 = product.compile(one, [], o);
P2 = product.compile(two, [], o);
P3 = product.compile(three, [], o);
Pc = pool.compile(chain_cfg, [], o);
Ps1 = pool.compile(one, rmfield(fac_cfg, 'tau_h'), o);
Pcf = pool.compile(chain_cfg, fac_cfg, o);
Pfo = product.compile([], fac_out, o);
% the same releasable-pool refill and release time under every model
same_first = isequal(P1.tau_rest, P2.tau_rest(:, 1), P3.tau_rest(:, 1), Pc.tau_rec1, Ps1.tau_rec1, Pcf.tau_rec1) && ...
    isequal(P1.tau_drive, P2.tau_drive(:, 1), P3.tau_drive(:, 1), Pc.tau_rel, Ps1.tau_rel, Pcf.tau_rel) && ...
    isequal(P1.tau_rest, 2 * exp(sg * zs(:, 1))) && isequal(P1.tau_drive, 0.25 * exp(sg * zs(:, 2)));
all_passed = check(['    tau_rec_1 and tau_rel_1 are the same per neuron on a one-component route, a product of ' ...
    'two and of three, a single pool and a chain'], same_first) && all_passed;
% the second refill time: one standard normal for the product and the chain
free = ~Pc.clipped;
same_second = isequal(P2.tau_rest(:, 2), 5 * exp(sg * zs(:, 3))) && isequal(P3.tau_rest(:, 2), P2.tau_rest(:, 2)) && ...
    isequal(Pc.tau_rec2(free), P2.tau_rest(free, 2)) && isequal(Pcf.tau_rec2, Pc.tau_rec2);
cols_ok = isequal(P2.tau_drive(:, 2), 0.4 * exp(sg * zs(:, 4))) && isequal(Pc.n2, 2 * exp(sg * zs(:, 5))) && ...
    isequal(P3.tau_rest(:, 3), 9 * exp(sg * zs(:, 6))) && isequal(P3.tau_drive(:, 3), 0.3 * exp(sg * zs(:, 7)));
all_passed = check(['    tau_rec_2 comes from one standard normal on the product and on the chain; tau_rel_2, n2 ' ...
    'and a third component read columns 4, 5 and 6, 7'], same_second && cols_ok) && all_passed;
% facilitation's columns, and the ceiling
fcols = isequal(Pcf.tau_dec, 1.5 * exp(sg * zf(:, 1))) && isequal(Pcf.tau_fac, 0.25 * exp(sg * zf(:, 2))) && ...
    isequal(Pcf.tau_h, 0.5 * exp(sg * zf(:, 3))) && isequal(Pcf.G, 1 + 4 * exp(sg * zf(:, 4))) && ...
    isequal(Pfo.tau_rest, [Pcf.tau_dec, 4 * exp(sg * zf(:, 5))]) && ...
    isequal(Pfo.tau_drive, [Pcf.tau_fac, 1 * exp(sg * zf(:, 6))]) && ...
    isequal(Pfo.target, [Pcf.G, 1 + 0.5 * exp(sg * zf(:, 7))]) && isequal(Ps1.G, Pcf.G);
G_one = pool.compile([], struct('model', 'chain', 'tau_dec', 1, 'tau_fac', 1, 'G', 1), o);
all_passed = check(['    facilitation: tau_dec, tau_fac, tau_h and G - 1 read columns 1 to 4, a second component ' ...
    '5 to 7; G >= 1 for every neuron and G = 1 stays 1'], ...
    fcols && all(Pcf.G >= 1) && all(Pfo.target(:) >= 1) && all(G_one.G == 1)) && all_passed;
% every drawn value positive, with the nominal as its median and sigma as its log-spread
lg = @(v, nominal) [mean(log(v / nominal)), std(log(v / nominal))];
stats = [lg(Pc.tau_rec1, 2); lg(Pc.tau_rel, 0.25); lg(Pc.n2, 2); lg(P2.tau_rest(:, 2), 5); ...
    lg(Pcf.tau_dec, 1.5); lg(Pcf.tau_h, 0.5); lg(Pcf.G - 1, 4)];
se = sg / sqrt(nn);
all_passed = check(sprintf(['    drawn values are positive; log-mean within 4 standard errors of the nominal ' ...
    '(max %.2g of %.2g) and log-spread within 5%% of sigma'], max(abs(stats(:, 1))), 4 * se), ...
    all([Pc.tau_rec1; Pc.tau_rel; Pc.tau_rec2; Pc.n2; Pcf.tau_dec; Pcf.tau_fac; Pcf.tau_h] > 0) && ...
    all(abs(stats(:, 1)) <= 4 * se) && all(abs(stats(:, 2) / sg - 1) <= 0.05)) && all_passed;
% the ordering rule: the reserve refill alone is raised to the floor
span = log(5 / 2);
gap_drawn = span + sg * (zs(:, 3) - zs(:, 1));
should = gap_drawn < span / 4;
gap = log(Pc.tau_rec2 ./ Pc.tau_rec1);
ok = isequal(Pc.clipped, should) && nnz(should) > 20 && nnz(should) < 0.06 * nn && ...
    all(gap >= span / 4 - 1e-12) && max(abs(gap(should) - span / 4)) <= 1e-12 && ...
    isequal(Pc.tau_rec1, 2 * exp(sg * zs(:, 1))) && all(Pc.tau_rec2(should) > 5 * exp(sg * zs(should, 3)));
all_passed = check(sprintf(['    the chain''s ordering rule: %d of %d neurons below a quarter of the nominal log-span ' ...
    'have tau_rec_2 raised to it exactly; tau_rec_1 is as drawn for all'], nnz(should), nn), ok) && all_passed;
% no ordering where the nominal reserve is the faster pool, on a product route or with no spread
fast_reserve = pool.compile(struct('model', 'chain', 'tau_rec', [3 0.8], 'tau_rel', 0.3, 'n2', 0.5), [], o);
o_zero = o; o_zero.std_sigma = 0;
unspread = pool.compile(chain_cfg, [], o_zero);
all_passed = check('    no neuron is moved where the nominal reserve is the faster pool, nor with no spread; a product route has no ordering', ...
    ~any(fast_reserve.clipped) && isequal(fast_reserve.tau_rec2, 0.8 * exp(sg * zs(:, 3))) && ...
    ~any(unspread.clipped) && all(unspread.tau_rec2 == 5) && ~isfield(P2, 'clipped') && ...
    any(P2.tau_rest(:, 2) < P2.tau_rest(:, 1))) && all_passed;
% the display factors take per-neuron parameters as they come
b = 0.2 + 0.7 * rand(nn, 2, 3);
o_floor = o; o_floor.zero_floor = true;
Pfl = product.compile(two, [], o_floor);
k = Pfl.tau_rest ./ Pfl.tau_drive;
p_min = prod(1 ./ (1 + k), 2);
D = product.depression(b, Pfl);
g = 1 + rand(nn, 2, 3);
ok = max(abs(Pfl.p_min - p_min)) <= 1e-15 && std(p_min) > 0.1 * mean(p_min) && ...
    max(abs(D - (prod(b, 2) - p_min) ./ (1 - p_min)), [], 'all') <= 1e-14 && ...
    isequal(product.facilitation(g, Pfo), prod(g, 2)) && isequal(pool.facilitation(g, Pcf), g(:, 1, :)) && ...
    isequal(pool.depression(b, Pc), b(:, 1, :));
all_passed = check('    the zero floor''s p_min is each neuron''s own, and depression and facilitation (both models have the handle) use it by row', ...
    ok) && all_passed;
% refusals
o_none = rmfield(o, 'std_z');
o_short = o; o_short.std_z = zs(1:10, :);
o_narrow = o; o_narrow.std_z = zs(:, 1:4);
ok = throws_id(@() product.compile(one, [], o_none), 'route_product:SpreadDraws') && ...
    throws_id(@() pool.compile(chain_cfg, [], o_short), 'route_pool:SpreadDraws') && ...
    throws_id(@() pool.compile(chain_cfg, [], o_narrow), 'route_pool:SpreadDraws') && ...
    throws_id(@() product.compile(three, [], o_narrow), 'route_product:SpreadDraws') && ...
    ~throws_id(@() product.compile(two, [], o_narrow), 'route_product:SpreadDraws') && ...
    ~throws_id(@() product.compile([], fac_out, o_none), 'route_product:SpreadDraws');
all_passed = check('    a spread without draws, with too few rows or too few columns errors <prefix>:SpreadDraws; a mechanism the route lacks needs none', ...
    ok) && all_passed;

%% 4. the class's two ways of applying A agree
% jacobian_times applies a route's A = d(ds)/ds in one of two ways: only
% its diagonal where the model states that A is diagonal (P.coupling), the
% full per-neuron product otherwise. The product model never takes the
% second way, so its routes are put through it here by telling the class
% that every route is coupled: the product must not change, and must still
% be J * Y. (The pool model always takes it: test_std_chain_model.)
fprintf('  SRNNCellTypePairs.jacobian_times, diagonal and full application of A\n');
rng(3, 'twister');
net = all_blocks_network();
net.build();
p = net.get_params();
p_full = p;
n_diagonal = 0;
for k = 1:numel(p_full.route)
    n_diagonal = n_diagonal + (p_full.route{k}.A_is_diagonal && p_full.route{k}.P.ns > 1);
    p_full.route{k}.A_is_diagonal = false;
end
N = numel(net.S0);
e_two = 0;
e_J = 0;
for trial = 1:3
    S = net.S0 + 0.2 * rand(N, 1) .* (net.S0 > 0);        % b, g off their rest values
    S(p.state_layout.x) = 0.5 * randn(p.n, 1);
    Y = randn(N, 6);
    Z_diag = SRNNCellTypePairs.jacobian_times(S, Y, p);
    Z_full = SRNNCellTypePairs.jacobian_times(S, Y, p_full);
    JY = SRNNCellTypePairs.compute_Jacobian_fast(S, p) * Y;
    e_two = max(e_two, rel_err(Z_full, Z_diag));
    e_J = max(e_J, rel_err(Z_full, JY));
end
all_passed = check(sprintf(['    %d multi-variable routes: full application equals the diagonal one (rel %.1e) ' ...
    'and J * Y (rel %.1e)'], n_diagonal, e_two, e_J), n_diagonal >= 2 && e_two <= 1e-14 && e_J <= 1e-13) && all_passed;

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
% Three of each, with a scalar tau_rel: no preset has three components, and
% the leave-one-out products are a loop only from three columns on.
std3 = struct('tau_rec', [2 4 8], 'tau_rel', 0.3, 'combine', 'geomean');
std3_alt = struct('tau_rec', [1 3 9], 'tau_rel', 0.2, 'combine', 'geomean');
stf3 = struct('tau_dec', [1 2 4], 'tau_fac', [0.25 0.5 1], 'G', [3 2 1.5]);
stf3_alt = struct('tau_dec', [1.5 3 5], 'tau_fac', [0.5 0.25 2], 'G', [2 4 3]);
cases = add(cases, 'three depression (geometric mean, zero floor), three facilitation', ...
    std3, stf3, floor_on, std3_alt, stf3_alt, 3, 3, 350);

% The pool model: the chain, a releasable pool and its reserve, both in
% [0, 1]. The paper's values, then a reserve smaller than the releasable
% pool that refills faster than it, so no ordering of the two is assumed.
draw_pool = @(npre, ns) 0.05 + 0.9 * rand(npre, ns);
add_pool = @(cases, label, std_cfg, std_alt, T) [cases, struct( ...
    'label', ['pool: ' label], 'model', 'pool', 'std', std_cfg, 'stf', none, 'opts', struct(), ...
    'std_alt', std_alt, 'stf_alt', none, 'initial', [1 1], ...
    'draw', draw_pool, 'reference', @chain_reference, 'T_settle', T)];
chain = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
chain_alt = struct('model', 'chain', 'tau_rec', [1.5 7], 'tau_rel', 0.4, 'n2', 3.5);
small = struct('model', 'chain', 'tau_rec', [3 0.8], 'tau_rel', 0.3, 'n2', 0.5);
small_alt = struct('model', 'chain', 'tau_rec', [2 1.2], 'tau_rel', 0.2, 'n2', 0.8);
cases = add_pool(cases, 'chain, reserve twice the releasable pool', chain, chain_alt, 300);
cases = add_pool(cases, 'chain, small fast reserve', small, small_alt, 200);

% The pool model's other variants: the single pool, and release-coupled
% facilitation (a release probability g in [1, G], and with tau_h a
% calcium trace h in [0, 1]) on no depression, the single pool and the
% chain. Tom's values (Plan A, section 4.6), and an alternative in which
% nothing equals anything else. G is at most 5 here.
fac = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5);
fac_alt = struct('model', 'chain', 'tau_dec', 0.9, 'tau_fac', 0.4, 'G', 2.86);
fac_h = struct('model', 'chain', 'tau_dec', 1.5, 'tau_fac', 0.25, 'G', 5, 'tau_h', 0.5);
fac_h_alt = struct('model', 'chain', 'tau_dec', 0.9, 'tau_fac', 0.4, 'G', 2.86, 'tau_h', 0.2);
draw_variant = @(n_b, has_g, has_h) @(npre, ns) [0.05 + 0.9 * rand(npre, n_b), ...
    1 + 4 * rand(npre, has_g), rand(npre, has_h)];
add_variant = @(cases, label, std_cfg, stf_cfg, std_alt, stf_alt, n_b, has_g, has_h, T) [cases, struct( ...
    'label', ['pool: ' label], 'model', 'pool', 'std', std_cfg, 'stf', stf_cfg, 'opts', struct(), ...
    'std_alt', std_alt, 'stf_alt', stf_alt, 'initial', [ones(1, n_b + has_g), zeros(1, has_h)], ...
    'draw', draw_variant(n_b, has_g, has_h), 'reference', @pool_reference, 'T_settle', T)];
cases = add_variant(cases, 'single pool', std1, none, std1_alt, none, 1, 0, 0, 80);
cases = add_variant(cases, 'g alone (rate-driven)', none, fac, none, fac_alt, 0, 1, 0, 60);
cases = add_variant(cases, 'g and calcium trace h', none, fac_h, none, fac_h_alt, 0, 1, 1, 60);
cases = add_variant(cases, 'single pool, g', std1, fac, std1_alt, fac_alt, 1, 1, 0, 100);
cases = add_variant(cases, 'single pool, g, h', std1, fac_h, std1_alt, fac_h_alt, 1, 1, 1, 100);
cases = add_variant(cases, 'chain, g', chain, fac, chain_alt, fac_alt, 2, 1, 0, 300);
cases = add_variant(cases, 'chain, g, h', chain, fac_h, chain_alt, fac_h_alt, 2, 1, 1, 300);
cases = add_variant(cases, 'small fast reserve, g, h', small, fac_h_alt, small_alt, fac_h, 2, 1, 1, 200);
end

function [ds, F] = pool_reference(s, r, cs)
% Every variant of the pool model for ONE neuron, from the equations of
% docs/Plan_A_chained_vessicles.md, section 4.6, written with an explicit
% list of which variables the route has (and NOT by setting an absent
% variable to a constant in the full equations):
%   db1/dt = b2 (1 - b1)/tau_rec1 - g b1 r/tau_rel       b2 absent: (1 - b1)/tau_rec1
%                                                        g absent:  b1 r/tau_rel
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2) b2 (1 - b1)/tau_rec1
%   dg/dt  = (1 - g)/tau_dec + (G - g) d/tau_fac         d = h, or r when h is absent
%   dh/dt  = (r - h)/tau_h
%   F = g b1; b1 when g is absent; g when b1 is absent
has_b1 = ~isempty(cs.std);
has_b2 = has_b1 && numel(cs.std.tau_rec) == 2;
has_g = ~isempty(cs.stf);
has_h = has_g && isfield(cs.stf, 'tau_h');
ds = zeros(1, has_b1 + has_b2 + has_g + has_h);
col = 0;
F = 1;
if has_g
    g = s(has_b1 + has_b2 + 1);
end
if has_b1
    b1 = s(1);
    col = col + 1;
    if has_g
        used = g * b1 * r / cs.std.tau_rel;
    else
        used = b1 * r / cs.std.tau_rel;
    end
    if has_b2
        b2 = s(2);
        ds(col) = b2 * (1 - b1) / cs.std.tau_rec(1) - used;
        col = col + 1;
        ds(col) = (1 - b2) / cs.std.tau_rec(2) - (1 / cs.std.n2) * b2 * (1 - b1) / cs.std.tau_rec(1);
    else
        ds(col) = (1 - b1) / cs.std.tau_rec(1) - used;
    end
    F = F * b1;
end
if has_g
    col = col + 1;
    if has_h
        h = s(col + 1);
        ds(col) = (1 - g) / cs.stf.tau_dec + (cs.stf.G - g) * h / cs.stf.tau_fac;
        ds(col + 1) = (r - h) / cs.stf.tau_h;
    else
        ds(col) = (1 - g) / cs.stf.tau_dec + (cs.stf.G - g) * r / cs.stf.tau_fac;
    end
    F = F * g;
end
end

function [ds, F] = chain_reference(s, r, cs)
% The chained pools for ONE neuron, from the equations (C1) to (C3) of
% docs/Plan_A_chained_vessicles.md, section 1:
%   db1/dt = b2 (1 - b1)/tau_rec1 - b1 r/tau_rel
%   db2/dt = (1 - b2)/tau_rec2 - (1/n2) b2 (1 - b1)/tau_rec1
%   F = b1
tau_rec1 = cs.std.tau_rec(1);
tau_rec2 = cs.std.tau_rec(2);
b1 = s(1);
b2 = s(2);
ds = zeros(1, 2);
ds(1) = b2 * (1 - b1) / tau_rec1 - b1 * r / cs.std.tau_rel;
ds(2) = (1 - b2) / tau_rec2 - (1 / cs.std.n2) * b2 * (1 - b1) / tau_rec1;
F = b1;
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

function cs = jittered_config(cs, sig_std, z_std, sig_stf, z_stf)
% The configuration ONE neuron has under a per-neuron spread, from its row
% of standard normals: written from the rule and the column table in
% route_model's header ("THE SPREAD"), not from the models.
%   value = nominal * exp(sigma * z(column));  G = 1 + (G - 1) * exp(sigma * z(column))
%   depression    1 tau_rec_1, 2 tau_rel_1, 3 tau_rec_2, 4 tau_rel_2, 5 n2,
%                 6, 7 tau_rec_3, tau_rel_3, ...
%   facilitation  1 tau_dec_1, 2 tau_fac_1, 3 tau_h, 4 G_1 - 1,
%                 5, 6, 7 tau_dec_2, tau_fac_2, G_2 - 1, ...
% and the chain's ordering rule: where log(tau_rec_2 / tau_rec_1) falls
% below a quarter of the nominal log-span (nominal reserve the slower),
% tau_rec_2 alone is raised to that floor.
is_chain = ~isempty(cs.std) && isfield(cs.std, 'model') && strcmp(cs.std.model, 'chain');
if ~isempty(cs.std)
    nominal = cs.std;
    M = numel(nominal.tau_rec);
    tau_rel = nominal.tau_rel .* ones(1, M);
    if is_chain
        tau_rel = nominal.tau_rel;                 % one release time
    end
    for m = 1:M
        if m <= 2, col = 2 * m - 1; else, col = 2 * m; end
        cs.std.tau_rec(m) = nominal.tau_rec(m) * exp(sig_std * z_std(col));
        if m <= numel(tau_rel)
            tau_rel(m) = tau_rel(m) * exp(sig_std * z_std(col + 1));
        end
    end
    cs.std.tau_rel = tau_rel;
    if is_chain
        cs.std.n2 = nominal.n2 * exp(sig_std * z_std(5));
        span = log(nominal.tau_rec(2) / nominal.tau_rec(1));
        if span > 0 && log(cs.std.tau_rec(2) / cs.std.tau_rec(1)) < span / 4
            cs.std.tau_rec(2) = cs.std.tau_rec(1) * exp(span / 4);
        end
    end
end
if ~isempty(cs.stf)
    nominal = cs.stf;
    L = numel(nominal.tau_dec);
    tau_fac = nominal.tau_fac .* ones(1, L);
    G = nominal.G .* ones(1, L);
    for l = 1:L
        if l == 1, cols = [1 2 4]; else, cols = 3 * l - 1 + (0:2); end
        cs.stf.tau_dec(l) = nominal.tau_dec(l) * exp(sig_stf * z_stf(cols(1)));
        tau_fac(l) = tau_fac(l) * exp(sig_stf * z_stf(cols(2)));
        G(l) = 1 + (G(l) - 1) * exp(sig_stf * z_stf(cols(3)));
    end
    cs.stf.tau_fac = tau_fac;
    cs.stf.G = G;
    if isfield(nominal, 'tau_h')
        cs.stf.tau_h = nominal.tau_h * exp(sig_stf * z_stf(3));
    end
end
end

function L = linearized(model, s, r, P)
% The four outputs of linearize, as one struct.
[L.A, L.B, L.F, L.Fs] = model.linearize(s, r, P);
end

function model = all_blocks_network()
% Three cell types with every kind of product route: two depression
% components, one, depression with two facilitation components,
% facilitation alone, none; the geometric mean on one; the zero floor on.
synapse_config = struct();
synapse_config.E.E.std = struct('tau_rec', [0.3 1], 'tau_rel', 0.2, 'combine', 'geomean');
synapse_config.E.PV.stf = struct('tau_dec', 0.5, 'tau_fac', 0.25, 'G', 2);
synapse_config.E.SST.std = struct('tau_rec', 0.8, 'tau_rel', 0.25);
synapse_config.E.SST.stf = struct('tau_dec', [0.4 1.2], 'tau_fac', [0.2 0.3], 'G', [1.5 2]);
synapse_config.PV.PV.std = struct('tau_rec', 0.6, 'tau_rel', 0.3);
synapse_config.SST.PV.stf = struct('tau_dec', 0.9, 'tau_fac', 0.4, 'G', 1.8);
model = SRNNCellTypePairs( ...
    'n', 12, 'indegree', 4, 'n_cellTypes', 3, 'cell_type_names', {'E', 'PV', 'SST'}, ...
    'f', [0.5 0.3 0.2], 'mu_tilde_relative', [0.1 -0.2 -0.08], ...
    'sigma_tilde_relative', [0.01 0.02 0.015], 'tau_a', {[0.25 1], [], 2}, ...
    'c', [0.05 0 0.03], 'synapse_config', synapse_config, 'std_zero_floor', true, ...
    'T_range', [0 0.1], 'fs', 200, 'lya_method', 'none', 'verbose', 'near-none');
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
