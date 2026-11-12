function paths = write_manuscript_tables(cfg)
% WRITE_MANUSCRIPT_TABLES Generate the manuscript's equation and parameter tables.
%
%   paths = WRITE_MANUSCRIPT_TABLES()
%   paths = WRITE_MANUSCRIPT_TABLES('preset_name', p)
%   paths = WRITE_MANUSCRIPT_TABLES('out_dir', d)
%
% Writes four files DIRECTLY into out_dir:
%   equation_table.md               equations + parameter table, annotated
%   adaptation_conditions.md        the regimes as they run, annotated
%   parameters_table.md             the parameter table ALONE
%   adaptation_conditions_table.md  the conditions table ALONE
% The bare pair carry no title, provenance banner or prose -- they are for
% pasting into a manuscript -- and share their row writers with the annotated
% pair, so the two versions of a table cannot disagree.
%
% THESE ARE PIPELINE OUTPUT, not documentation checked into the tree. They are
% rebuilt from a live model on every run, describe one preset, and are
% regenerable -- the same properties a figure has. So they go where the figures
% go: fig_doc_tables makes them a registry entry, and make_all_paper_figures
% hands them <fig_root>/doc_tables. They used to be written into a
% doc_equations_table/ folder derived from this file's own location, which meant
% the one output make_all_paper_figures produced OUTSIDE its output root.
%
% The PDF renders are not generated here. Run pandoc against these .md files if
% you want them.
%
% It used to write a third file, fig_equations/parameter_table.md, whose table
% was IDENTICAL row-for-row to equation_table.md's "Parameters as run" section.
% Dropped 2026-08-27 along with the fig_equations/ folder, which held no fig_*.m
% and so was never a figure in the registry.
%
% WHY THESE ARE GENERATED. They were hand-written and had gone stale in a way
% that mattered: equation_table.md documented a SINGLE STD timescale, the
% LOGISTIC activation and c_E = 0.5/3 on SRNNModel2, while the paper's preset is
% DUAL STD, piecewise, on SRNNCellTypePairs. adaptation_conditions.md listed
% n_a_E / n_a_I / n_b_E / n_b_I counts, which are not SRNNCellTypePairs
% properties at all -- that class names routes individually. A reader checking
% the methods against the code would have found neither matching.
%
% Every value below is read off a BUILT model, under the sfa_and_std condition,
% so the tables cannot disagree with what the sweeps ran. This is the same
% principle write_run_parameters_md applies to a run directory.
%
% WHAT IS NOT GENERATED. The prose explaining what each mechanism MEANS stays
% hand-written -- generation is for the numbers and the shapes, which drift.
% The equations themselves are emitted because their FORM depends on the preset
% (how many STD timescales, whether STF is present, whether there is noise).
%
% See also: srnn_param_preset, write_run_parameters_md, build_from_preset

arguments
    cfg.preset_name (1,:) char = 'celltype_pairs_Sc0p2_noise0p025_dualStd_7cond'
    % Where the two .md files land. Empty defaults to figs/doc_tables/ -- under
    % figs/ like every other output, and deliberately NOT next to this .m file,
    % which is how the pipeline's tables ended up outside its own output root.
    cfg.out_dir     (1,:) char = ''
    cfg.verbose     (1,:) char    = 'minimal'   % 'verbose' | 'minimal' | 'near-none' (see verbose_level)
    % Condition names the tables cover, in order; empty = every condition of
    % the preset (select_conditions). See "WHICH CONDITIONS" below.
    cfg.conditions                = {}
end
% WHICH CONDITIONS (FD-006). The equations and the conditions table cover the
% covered conditions: every one by default, the selection when cfg.conditions
% names one. The parameter table describes ONE built model, under the most
% adapted of the covered conditions (full_adaptation_condition: the chain wins
% a tie with the product). The depression the equations state is read off the
% compiled routes of the covered conditions (route_model: P.model,
% P.nominal), never off a condition's name:
%   * product routes (independent pools, route_product) -> the product
%     equation for b_im; the "product SQUARES the depression" prose appears
%     only when a covered condition has a product of TWO OR MORE pools;
%   * chain routes (route_pool) -> (C1)-(C3) of
%     docs/EquationsParametersDocs/Equations_stability_paper.md, the chain's
%     parameters (both refill times, the one release time, n2) and its steady
%     state (C4), evaluated with the model's own steady_state handle.

setup_paths();
if isempty(cfg.out_dir)
    project_root = fileparts(which('setup_paths'));
    out_dir = fullfile(project_root, 'figs', 'doc_tables');
else
    out_dir = cfg.out_dir;
end
out_dir = ensure_dir(out_dir);

[preset, model_class, all_conditions] = srnn_param_preset(cfg.preset_name);
all_names = cellfun(@(c) c.name, all_conditions, 'UniformOutput', false);
[~, keep] = select_conditions(all_names, cfg.conditions);
conditions = all_conditions(keep);
% The most-adapted COVERED regime, whatever it is called. The tables document
% the fullest configuration, and which condition that is varies by preset
% since names now state their timescale counts. A selection with no
% conforming name (no_adaptation alone) falls back to the preset's own.
try
    full_name = full_adaptation_condition(conditions);
catch
    full_name = full_adaptation_condition(all_conditions);
end
model = build_from_preset(cfg.preset_name, full_name);
M = model_facts(model, model_class, preset);
% The regime the model above was built under, for the provenance banner. This
% was the literal 'sfa_and_std' for a week after that name stopped existing.
M.condition = full_name;
M.n_all_conditions = numel(all_conditions);
M = depression_facts(M, conditions, model);

paths = {};
paths{end+1} = write_equation_table(out_dir, cfg.preset_name, M);
paths{end+1} = write_conditions_table(out_dir, cfg.preset_name, conditions, M);
% The same two tables with nothing around them, for pasting into a manuscript.
paths{end+1} = write_param_table_bare(out_dir, M);
paths{end+1} = write_conditions_table_bare(out_dir, conditions, M);

vprintf(cfg.verbose, 'verbose', '[manuscript tables] preset %s (%s)\n', cfg.preset_name, model_class);
vprintf(cfg.verbose, 'verbose', '  %s\n', paths{:});
end

%% ------------------------------------------------------------------------
function M = model_facts(model, model_class, preset)
% Everything the tables need, read off the BUILT model.
M = struct();
M.model_class = model_class;
M.n           = model.n;
M.indegree    = model.indegree;
M.tau_d       = model.tau_d;
M.activation  = model.activation;
M.S_a         = model.S_a;
M.S_c         = model.S_c;
% Per-neuron setpoint spread, one entry per cell type; zeros when unset. On
% Pairs build() fills sigma_S_c in as zeros(1, C); SRNNModel2 keeps E and I
% scalars. mu_S_c empty means the draw is centred on S_c.
if strcmp(model_class, 'SRNNCellTypePairs')
    M.sigma_S_c = get_or(model, 'sigma_S_c', []);
else
    M.sigma_S_c = [get_or(model, 'sigma_S_c_E', 0), get_or(model, 'sigma_S_c_I', 0)];
end
if isempty(M.sigma_S_c); M.sigma_S_c = 0; end
% Per-neuron SFA ladder spread (Pairs only; build() fills it as zeros(1, C)).
M.tau_a_spread = get_or(model, 'tau_a_spread', 0);
if isempty(M.tau_a_spread); M.tau_a_spread = 0; end
M.level_of_chaos = model.level_of_chaos;
M.fs          = model.fs;
M.sigma_u_noise = get_or(model, 'sigma_u_noise', 0);
M.is_pairs    = strcmp(model_class, 'SRNNCellTypePairs');

if M.is_pairs
    M.types   = model.cell_type_names;
    M.f       = model.f;
    M.n_a     = model.n_a;
    M.tau_a   = model.tau_a;
    M.c       = model.c;
    M.mu      = model.mu_tilde_relative;
    M.sigma   = model.sigma_tilde_relative;
    M.routes  = route_facts(model);
else
    M.types   = {'E', 'I'};
    M.f       = [model.f, 1 - model.f];
    M.n_a     = [model.n_a_E, model.n_a_I];
    M.tau_a   = {model.tau_a_E, model.tau_a_I};
    M.c       = [model.c_E, model.c_I];
    M.mu      = [model.mu_E_tilde_relative, model.mu_I_tilde_relative];
    M.sigma   = [model.sigma_E_tilde_relative, model.sigma_I_tilde_relative];
    M.routes  = struct('label', {'E (all outgoing)'}, ...
        'tau_rec', {model.tau_b_E_rec}, 'tau_rel', {model.tau_b_E_rel}, ...
        'has_stf', {false}, 'tau_dec', {[]}, 'tau_fac', {[]}, 'G', {[]});
end

M.n_b_max = 0;
M.n_g_max = 0;
M.any_stf = false;
for k = 1:numel(M.routes)
    M.n_b_max = max(M.n_b_max, numel(M.routes(k).tau_rec));
    M.n_g_max = max(M.n_g_max, numel(M.routes(k).tau_dec));
    M.any_stf = M.any_stf || M.routes(k).has_stf;
end
M.has_noise = any(M.sigma_u_noise(:) > 0);
M.preset_fields = fieldnames(preset);
end

function M = depression_facts(M, conditions, model)
% What depression the COVERED conditions carry, from their compiled routes.
%   M.n_b_max       the widest PRODUCT route (the M of the product equation),
%                   over the covered conditions; 0 with no product depression
%   M.any_product   some covered route has product depression (n_b >= 1)
%   M.any_chain     some covered route is a chain
%   M.chains        the distinct chain routes (one per P.signature): their
%                   nominal values and steady state, for the prose
%   M.cond_dep      per covered condition, the table column's text
%   M.any_stf       also true when a covered condition has facilitation
% For SRNNModel2 (no routes) the built model's facts stand as they are.
M.any_chain = false;
M.product_ratios = [];   % tau_rec./tau_rel of the first covered product route of n_b_max pools
M.chains = struct('nominal', {}, 'probe_r', {}, 'b1', {}, 'b2', {}, 'ceiling', {}, 'conditions', {});
M.cond_dep = repmat({'none'}, 1, numel(conditions));
if ~M.is_pairs
    M.any_product = M.n_b_max > 0;
    if M.n_b_max > 0
        M.product_ratios = M.routes(1).tau_rec ./ M.routes(1).tau_rel;
    end
    return
end
% The rep model's own routes decide M.n_b_max as before, widened by any
% covered condition's product routes.
M.n_b_max = 0;
for k = 1:numel(M.routes)
    if ~M.routes(k).is_chain
        M.n_b_max = max(M.n_b_max, numel(M.routes(k).tau_rec));
    end
    M.any_chain = M.any_chain || M.routes(k).is_chain;
end
probe_r = [0.1 0.25 0.5 1];
best_nb = 0;
zero_floor =isprop(model, 'std_zero_floor') && model.std_zero_floor;
for c = 1:numel(conditions)
    routes = condition_routes(conditions{c}, zero_floor);
    texts = cell(1, numel(routes));
    for k = 1:numel(routes)
        P = routes(k).P;
        if P.n_b > 0 && strcmp(P.model, 'product')
            M.n_b_max = max(M.n_b_max, P.n_b);
            if P.n_b > best_nb
                best_nb = P.n_b;
                M.product_ratios = P.nominal.tau_rec ./ P.nominal.tau_rel;
            end
        end
        M.any_stf = M.any_stf || P.n_g > 0;
        if strcmp(P.model, 'pool')
            M.any_chain = true;
            j = find(arrayfun(@(x) isequal(x.nominal, P.signature), M.chains), 1);
            if isempty(j)
                [b1, b2, ceiling] = chain_steady_state(routes(k), probe_r);
                M.chains(end+1) = struct('nominal', P.signature, 'probe_r', probe_r, ...
                    'b1', b1, 'b2', b2, 'ceiling', ceiling, 'conditions', {{}});
                j = numel(M.chains);
            end
            if ~any(strcmp(M.chains(j).conditions, conditions{c}.name))
                M.chains(j).conditions{end+1} = conditions{c}.name;
            end
        end
        texts{k} = depression_text(P);
    end
    if isempty(routes)
        continue
    end
    labels = {routes.label};
    if isscalar(unique(texts))
        M.cond_dep{c} = texts{1};          % every depressing route alike
    else
        M.cond_dep{c} = strjoin(strcat(labels, {': '}, texts), '; ');
    end
end
M.any_product = M.n_b_max > 0;
end

function routes = condition_routes(cd_, zero_floor)
% The condition's depressing routes, compiled by the model class's own route
% compiler (SRNNCellTypePairs.compile_route), labelled "EE", "EI", ...
routes = struct('label', {}, 'model', {}, 'P', {});
if ~isfield(cd_, 'synapse_config') || ~isstruct(cd_.synapse_config); return; end
sc = cd_.synapse_config;
pres = fieldnames(sc);
for a = 1:numel(pres)
    posts = fieldnames(sc.(pres{a}));
    for b = 1:numel(posts)
        e = sc.(pres{a}).(posts{b});
        if ~(isfield(e, 'std') && ~isempty(e.std)); continue; end
        stf = [];
        if isfield(e, 'stf'); stf = e.stf; end
        R = SRNNCellTypePairs.compile_route(struct('std', e.std, 'stf', stf), ...
            struct('npre', 1, 'zero_floor', zero_floor, ...
            'label', sprintf('%s->%s', pres{a}, posts{b}), 'err_id_prefix', 'write_manuscript_tables'));
        routes(end+1) = struct('label', [pres{a} posts{b}], 'model', R.model, 'P', R.P); %#ok<AGROW>
    end
end
end

function [b1, b2, ceiling] = chain_steady_state(route, probe_r)
% The chain's steady state at the probe rates and its sustained-output
% ceiling, from the route model's steady_state handle (route_pool's (C4)).
b1 = zeros(size(probe_r)); b2 = b1;
for k = 1:numel(probe_r)
    ss = route.model.steady_state(probe_r(k), route.P);
    b1(k) = route.model.depression(ss.s(:, 1:route.P.n_b), route.P);
    b2(k) = ss.s(1, 2);
end
r_inf = 1e9;                       % r b1 -> the ceiling as r -> infinity
ss = route.model.steady_state(r_inf, route.P);
ceiling = r_inf * ss.F;
end

function s = depression_text(P)
% One route's depression, as the conditions table states it (from P.nominal).
n = P.nominal;
if strcmp(P.model, 'pool')
    s = sprintf('chained pools: $\\tau_{rec}$ = %s s, $\\tau_{rel}$ = %s s, $n_2$ = %g', ...
        vec(n.tau_rec), vec(n.tau_rel), n.n2);
elseif P.n_b == 1
    s = sprintf('one pool: $\\tau_{rec}$ = %s s, $\\tau_{rel}$ = %s s', vec(n.tau_rec), vec(n.tau_rel));
else
    s = sprintf('product of %d pools: $\\tau_{rec}$ = %s s, $\\tau_{rel}$ = %s s', ...
        P.n_b, vec(n.tau_rec), vec(n.tau_rel));
end
if strcmp(P.model, 'product') && isfield(n, 'combine') && strcmp(n.combine, 'geomean')
    s = [s ', geometric-mean combine'];
end
if P.n_g > 0
    s = sprintf('%s; facilitation $\\tau_{dec}$ = %s s, $\\tau_{fac}$ = %s s, $G$ = %s', ...
        s, vec(n.tau_dec), vec(n.tau_fac), vec(n.G));
end
end

function s = vec(v)
if isscalar(v)
    s = sprintf('%g', v);
else
    s = ['[' strjoin(arrayfun(@(x) sprintf('%g', x), v, 'UniformOutput', false), ', ') ']'];
end
end

function R = route_facts(model)
% One entry per depressing/facilitating route of the BUILT model. The values
% are the configuration as written; whether the route is a chain is the
% compiled route's (get_params: params.route{pre, post}.P.model).
R = struct('label', {}, 'tau_rec', {}, 'tau_rel', {}, ...
           'has_stf', {}, 'tau_dec', {}, 'tau_fac', {}, 'G', {}, 'is_chain', {}, 'n2', {});
sc = model.synapse_config;
if ~isstruct(sc); return; end
params = model.get_params();
type_names = model.cell_type_names;
pres = fieldnames(sc);
for a = 1:numel(pres)
    posts = fieldnames(sc.(pres{a}));
    for b = 1:numel(posts)
        e = sc.(pres{a}).(posts{b});
        has_std = isfield(e, 'std') && ~isempty(e.std);
        has_stf = isfield(e, 'stf') && ~isempty(e.stf);
        if ~has_std && ~has_stf; continue; end
        r = struct('label', sprintf('%s to %s', pres{a}, posts{b}), ...
            'tau_rec', [], 'tau_rel', [], 'has_stf', has_stf, ...
            'tau_dec', [], 'tau_fac', [], 'G', [], 'is_chain', false, 'n2', []);
        if has_std
            r.tau_rec = e.std.tau_rec(:)';
            r.tau_rel = e.std.tau_rel(:)';
            P = params.route{strcmp(type_names, pres{a}), strcmp(type_names, posts{b})}.P;
            r.is_chain = strcmp(P.model, 'pool');
            if r.is_chain; r.n2 = P.nominal.n2; end
        end
        if has_stf
            r.tau_dec = e.stf.tau_dec(:)';
            r.tau_fac = e.stf.tau_fac(:)';
            r.G       = e.stf.G(:)';
        end
        R(end+1) = r; %#ok<AGROW>
    end
end
end

%% ------------------------------------------------------------------------
function p = write_equation_table(dir_out, preset_name, M)
p = fullfile(dir_out, 'equation_table.md');
fid = open_md(p, 'SRNN model equations', preset_name, M);
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '## System equations\n\n');

% ONE ALIGNED BLOCK, IN THE ORDER AND NOTATION OF THE CANONICAL FILE
% docs/EquationsParametersDocs/Equations_stability_paper.md: dx with theta_j on
% the right-hand side, then theta, r, a and b. Timescale counts are the SYMBOLS
% M and K, not this preset's numbers -- the numbers are in the parameter table
% as the M and K rows. (Until 2026-09-10 the product was spelled out inline in
% the dx equation with numeric limits, and theta was named only afterwards;
% TR asked for the six equations to match the paper's statement.)
%
% What still tracks the preset: the noise term is present only for a stochastic
% preset, the b equation only when there is depression, and facilitation adds
% its factor and equation only when a route carries STF.
%
% Chained pools (FD-006): with chain routes alone theta's factor is b_{i1},
% the releasable pool, as the canonical file writes it; with product routes
% as well (1TS is a product of one pool), theta carries a depression factor
% D_i, and a line states D_i for each form. Product-only presets are
% written exactly as before.
theta_factors = {};
mixed = M.any_chain && M.any_product;
if mixed
    theta_factors{end+1} = 'D_{i}';
elseif M.any_chain
    theta_factors{end+1} = 'b_{i1}';
elseif M.n_b_max > 0
    theta_factors{end+1} = '\prod_{m=1}^{M} b_{im}';
end
if M.any_stf;     theta_factors{end+1} = '\prod_{n=1}^{n_f} g_{in}'; end
% Spelled exactly as the canonical file spells it -- r_{i} and plain spaces --
% so the generated block can be diffed against it.
if isempty(theta_factors)
    theta_rhs = 'r_{i}';
else
    theta_rhs = ['r_{i} ' strjoin(theta_factors, ' ')];
end

if M.has_noise
    dx_line = ['dx_i &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, \theta_j}' ...
        '{\tau_d}\, dt \;+\; \frac{\sigma_u}{\tau_d}\, dW_i'];
else
    dx_line = ['\frac{dx_i}{dt} &= \frac{-x_i + u_i + \sum_{j=1}^{N} w_{ij}\, ' ...
        '\theta_j}{\tau_d}'];
end
eqs = { ...
    dx_line, ...
    ['\theta_i &= ' theta_rhs]};
if mixed
    eqs{end+1} = ['D_{i} &= \prod_{m=1}^{M} b_{im} \quad \text{(independent pools)}, ' ...
        '\qquad D_{i} = b_{i1} \quad \text{(chained pools)}'];
end
eqs = [eqs, { ...
    'r_i &= \phi\left( x_i - a_{0_i} - \frac{c}{K} \sum_{k=1}^{K} a_{ik} \right)', ...
    '\frac{da_{ik}}{dt} &= \frac{-a_{ik} + r_i}{\tau_{a_k}}, \qquad k = 1, \dots, K'}];
if M.n_b_max > 0
    eqs{end+1} = ['\frac{db_{im}}{dt} &= \frac{1-b_{im}}{\tau_{rec_m}} - ' ...
        '\frac{b_{im}\, r_i}{\tau_{rel_m}}, \qquad m = 1, \dots, M'];
end
if M.any_chain
    % (C1) and (C2), spelled as the canonical file spells them.
    eqs{end+1} = ['\frac{db_{i1}}{dt} &= \frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}} - ' ...
        '\frac{b_{i1}\, r_i}{\tau_{rel}}'];
    eqs{end+1} = ['\frac{db_{i2}}{dt} &= \frac{1 - b_{i2}}{\tau_{rec_2}} - ' ...
        '\frac{1}{n_2}\,\frac{b_{i2}\,(1 - b_{i1})}{\tau_{rec_1}}'];
end
if M.any_stf
    eqs{end+1} = ['\frac{dg_{in}}{dt} &= \frac{1-g_{in}}{\tau_{dec_n}} + ' ...
        '\frac{(G - g_{in})\, r_i}{\tau_{fac_n}}, \qquad n = 1, \dots, n_f'];
end
% Lines are DATA to %s, so single backslashes; the '\\[8pt]' row separator is
% in the format string and so is doubled. Note '\\[8pt]' after every line but
% the last, exactly as the canonical block.
fprintf(fid, '$$\n\\begin{aligned}\n%s\n\\end{aligned}\n$$\n\n', ...
    strjoin(eqs, ' \\\\[8pt]\n'));

if M.has_noise
    fprintf(fid, ['Additive Wiener noise enters **only** $x$, which keeps the ' ...
        'diffusion constant: Ito and Stratonovich coincide, the Milstein term ' ...
        'vanishes, the variational equation is untouched, and the noise ' ...
        'cancels in Benettin''s trajectory difference. $\\sigma_u$ is ' ...
        'input-referred, so it is comparable to the stimulus amplitude.\n\n']);
end

% a_{0_i} is the setpoint (the property S_c), subtracted here so phi is the
% zero-centred function. The code passes the centre into phi instead
% (piecewiseSigmoid(x, S_a, S_c)); the two agree because a translation commutes,
% phi_{S_c}(z) = phi_0(z - S_c).
if any(M.sigma_S_c(:) > 0)
    setpoint_note = sprintf(['drawn per neuron as $a_{0_i} = %s + \\sigma_{S_c}' ...
        '\\,\\xi_i$, $\\xi_i \\sim \\mathcal{N}(0,1)$, with $\\sigma_{S_c}$ = %s ' ...
        'per cell type'], fmt(M.S_c), fmt(M.sigma_S_c));
else
    setpoint_note = sprintf('%s for every neuron here', fmt(M.S_c));
end
if M.any_chain && M.any_product
    m_phrase = ['$M$ the number of independent depression pools of a product ' ...
        'route ($M = 1$ is the single pool)'];
elseif M.any_chain
    m_phrase = '$b_{i1}$, $b_{i2}$ the two chained depression pools';
else
    m_phrase = '$M$ the number of depression timescales';
end
fprintf(fid, ['$\\phi$ is centred on zero; $a_{0_i}$ is the **setpoint** ' ...
    '(`S_c`), %s. $K$ is the number of SFA timescales and %s; their values ' ...
    'for this preset are in the parameter table below.\n\n'], setpoint_note, m_phrase);
fprintf(fid, ['The rate $r_i$ is the **pre-depression** output of the ' ...
    'nonlinearity. Depression and facilitation enter presynaptically, as ' ...
    'factors multiplying $r_j$ inside the recurrent sum, so both mechanisms are ' ...
    'driven by the raw rate $r_i$ and not by the depressed output.\n\n']);

if M.any_stf
    after = 'after depression and facilitation';
else
    after = 'after depression';
end
fprintf(fid, ['$\\theta_i$ is the **synaptic output**: the rate ' ...
    '%s, and the quantity the recurrent sum transmits.\n\n'], after);
if M.n_b_max > 1
    % Quote the actual depth this preset reaches. The ratio tau_rec/tau_rel is
    % the ONLY combination that sets the steady state, so it is what the numbers
    % below are computed from.
    % From the first covered PRODUCT route of the widest count (depression_facts);
    % a chain's two tau_rec share one tau_rel and are not two factors of a
    % product.
    ratios = M.product_ratios;
    lead = 'The';
    if M.any_chain
        % Mixed: say which routes this paragraph is about.
        lead = 'On a route with **independent** pools (not a chain), the';
    end
    fprintf(fid, '%s', lead);
    fprintf(fid, [' depression timescales enter as a **product**, which is ' ...
        'the reason for using more than one: %d of them SQUARE the depression ' ...
        'at steady state, so the extra timescale buys DEPTH and not merely a ' ...
        'slower approach. There is no split of $\\tau$ that would make them ' ...
        'match a single timescale.\n\n'], M.n_b_max);
    if ~isempty(ratios) && all(abs(ratios - ratios(1)) < 1e-12)
        rr = ratios(1);
        one = 1 / (1 + rr);
        many = 1 / (1 + rr)^M.n_b_max;
        if M.n_b_max == 2
            power_phrase = 'exactly the square';
        else
            power_phrase = sprintf('exactly the %dth power', M.n_b_max);
        end
        fprintf(fid, ['With $\\tau_{rec}/\\tau_{rel} = %g$ on every timescale, ' ...
            'each factor relaxes to $b_k(r) = 1/(1 + %g\\,r)$, so at $r = 1$ a ' ...
            'single timescale gives a gain of %.4g (a %.0f-fold reduction ' ...
            'against an undepressed synapse) while %d give %.4g -- a ' ...
            '**%.0f-fold** reduction, %s.\n\n'], ...
            rr, rr, one, 1/one, M.n_b_max, many, 1/many, power_phrase);
    end
    fprintf(fid, ['Adaptation, by contrast, enters as a **sum**, and is ' ...
        'normalised by its timescale count: the rate subtracts ' ...
        '$(c/K)\\sum_{k=1}^{K} a_k$. Since every $a_k$ relaxes to the rate, ' ...
        '$\\sum_k a_k \\to K r$ regardless of the $\\tau_k$, so the steady state ' ...
        'is $c\\,r$ whatever $K$ is -- adding SFA timescales changes the route, ' ...
        'not the destination. $c$ is therefore the **total** adaptation budget. ' ...
        'No such normalisation applies to depression, which multiplies rather ' ...
        'than sums.\n\n']);
end

if M.any_chain
    write_chain_prose(fid, M);
end

fprintf(fid, ['The sparse weights are drawn element-wise and scaled by the ' ...
    'synaptic gain $g$:\n\n$$\nw_{ij} = g\\, S_{ij}\\left( \\tilde{\\sigma}_{ij}\\, ' ...
    '\\xi_{ij} + \\tilde{\\mu}_{ij} \\right), \\qquad \\xi_{ij} \\sim ' ...
    '\\mathcal{N}(0,1), \\qquad S_{ij} \\sim \\mathrm{Bernoulli}(\\alpha)\n$$\n\n']);
fprintf(fid, ['with $(\\tilde{\\mu}, \\tilde{\\sigma})$ indexed by ' ...
    '**(postsynaptic, presynaptic)** cell type, and given as multiples of ' ...
    '$F = 1/\\sqrt{n\\,\\alpha(2-\\alpha)}$.\n\n']);

fprintf(fid, '## Parameters as run\n\n');
write_param_rows(fid, M);
end

function write_chain_prose(fid, M)
% What the chained pools are, and their steady state, for each distinct chain
% route of the covered conditions. The values are the route model's
% (route_pool: steady_state), evaluated in depression_facts.
fprintf(fid, ['On a **chained** route the two depression variables are pools in ' ...
    'series, not factors of a product: $b_{i1}$ is the **releasable pool**, the ' ...
    'only one released from and the only one in $\\theta_i$, and $b_{i2}$ the ' ...
    '**reserve** that refills it, at a rate proportional to how full the reserve ' ...
    'is; the reserve refills toward full with its own, slower time constant and ' ...
    'loses what it hands on, divided by its relative capacity $n_2$. Both start ' ...
    'at 1 and stay in $[0, 1]$.\n\n']);
if M.any_product
    fprintf(fid, ['With the reserve held full ($b_{i2} \\equiv 1$) the releasable-pool ' ...
        'equation is the single-pool equation ($M = 1$), so single-pool and chained ' ...
        'depression share the releasable-pool equation and differ by the reserve. ' ...
        'The reserve starts full, so the two coincide at the onset of activity and ' ...
        'come apart as sustained activity draws the reserve down. They are not ' ...
        'matched in steady-state depth.\n\n']);
end
fprintf(fid, ['At a constant rate $r$, with $k_1 = r\\,\\tau_{rec_1}/\\tau_{rel}$, ' ...
    '$k_2 = r\\,\\tau_{rec_2}/(n_2\\,\\tau_{rel})$ and $s = 1 + k_1 + k_2$, the ' ...
    'steady state is\n\n$$\nb_1 = \\frac{2}{s + \\sqrt{s^2 - 4k_2}}, \\qquad ' ...
    'b_2 = 1 - k_2\\, b_1 ,\n$$\n\n']);
fprintf(fid, ['and the sustained output $r\\, b_1$ rises with the rate to the ceiling ' ...
    '$\\big(\\tau_{rec_1}/\\tau_{rel} + \\tau_{rec_2}/(n_2\\,\\tau_{rel})\\big)^{-1}$; ' ...
    'it has no peak.\n\n']);
for j = 1:numel(M.chains)
    c = M.chains(j);
    fprintf(fid, ['With $\\tau_{rec} = %s$ s, $\\tau_{rel} = %s$ s and $n_2 = %g$ ' ...
        '(%s): $b_1$ = %s and $b_2$ = %s at $r$ = %s; the ceiling is %.4g.\n\n'], ...
        vec(c.nominal.tau_rec), vec(c.nominal.tau_rel), c.nominal.n2, ...
        strjoin(cellfun(@(n) ['`' n '`'], c.conditions, 'UniformOutput', false), ', '), ...
        num_list(c.b1, '%.3f'), num_list(c.b2, '%.3f'), num_list(c.probe_r, '%g'), c.ceiling);
end
end

function s = num_list(v, f)
s = strjoin(arrayfun(@(x) sprintf(f, x), v, 'UniformOutput', false), ', ');
end

function tf = route_is_chain(r)
% A route entry of M.routes is a chain (SRNNModel2's entries have no field).
tf = isfield(r, 'is_chain') && r.is_chain;
end

%% ------------------------------------------------------------------------
% A bare parameter table existed before, as fig_equations/parameter_table.md,
% and was removed on 2026-08-27 as a second copy of the section
% equation_table.md already carries. It is back (TR, 2026-09-09) for a
% different reason: a manuscript wants the TABLE ALONE, with none of the
% title, provenance banner or prose around it. Both files draw their rows from
% write_param_rows, so there is still only one place the table is spelled out.
function p = write_param_table_bare(dir_out, M)
% The parameter table and nothing else -- no title, no provenance, no heading
% -- for dropping straight into a manuscript. Shares write_param_rows with
% equation_table.md's "Parameters as run" section, so the two cannot disagree.
p = fullfile(dir_out, 'parameters_table.md');
fid = fopen(p, 'w');
if fid < 0; error('write_manuscript_tables:CannotOpen', 'Could not open %s', p); end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
write_param_rows(fid, M);
end

function write_param_rows(fid, M)
% Header row and one row per parameter; the only place the table is spelled
% out. The "## Parameters as run" heading is the annotated caller's, so this
% can also be written bare.
fprintf(fid, '| Symbol | Name | Value | Units |\n|---|---|---|---|\n');

row(fid, '$n$',          'Network size',              M.n,        'neurons');
% NOTE the single backslashes below. These strings are passed as DATA to a %s
% conversion, not as a format string, so they are written verbatim -- a '\\'
% here would land in the markdown as a literal double backslash and break the
% LaTeX. The sprintf() calls further down are the opposite case and correctly
% use '\\', because sprintf collapses it.
row(fid, '$\alpha$',    'Connection probability',    M.indegree / M.n, '--');
row(fid, 'indegree',     'Expected in-degree',        M.indegree, 'synapses');
row(fid, '$f$',          'Fraction per cell type',    M.f,        '--');
row(fid, '$\tau_d$',    'Dendritic time constant',   M.tau_d,    's');
row(fid, '$g$',          'Synaptic gain',             M.level_of_chaos, '--');
row(fid, '$\tilde{\mu}$',    'Weight means (post, pre)',  M.mu,    'multiples of $F$');
row(fid, '$\tilde{\sigma}$', 'Weight std devs (post, pre)', M.sigma, 'multiples of $F$');
row(fid, '$\phi$',      'Nonlinearity',              M.activation, '--');
if ~strcmp(M.activation, 'tanh')
    row(fid, '$S_c$',    'Nonlinearity setpoint',     M.S_c,      '--');
    % Per-neuron spread of the setpoint, one entry per cell type; 0 = shared.
    row(fid, '$\sigma_{S_c}$', 'Setpoint std per cell type', M.sigma_S_c, '--');
end
if strcmp(M.activation, 'piecewise')
    row(fid, '$S_a$',    'Piecewise slope parameter', M.S_a,      '--');
end
row(fid, '$K$',          'SFA timescales per type',   M.n_a,      '--');
for t = 1:numel(M.types)
    if M.n_a(t) > 0
        row(fid, sprintf('$\\tau_a$ (%s)', M.types{t}), 'SFA time constants', ...
            M.tau_a{t}, 's');
    end
end
% Per-neuron ladders: log-normal spread at both ends of the ladder, geometric
% between (Equations_stability_paper.md). ALWAYS printed (TR, 2026-09-14): a
% zero row says every neuron carries the nominal ladder, which a reader of the
% table could not otherwise tell from a preset that has no such knob.
row(fid, '$\sigma_{\tau_a}$', ...
    'SFA ladder spread per cell type (log-normal, both ends; 0 = every neuron has the nominal ladder)', ...
    M.tau_a_spread, '--');
row(fid, '$c$',          'SFA coupling per type',     M.c,        '--');
% M as in the equations: depression timescales. Uniform across routes in every
% preset so far; if a preset ever differs by route, this is the widest one.
if M.any_chain
    % With chains in play M counts only the pools of a PRODUCT route.
    row(fid, '$M$',      'STD pools of a product route (covered conditions)', M.n_b_max, '--');
else
    row(fid, '$M$',      'STD timescales',            M.n_b_max,  '--');
end
for k = 1:numel(M.routes)
    r = M.routes(k);
    if route_is_chain(r)
        % (C1)-(C3): two refill times, ONE release time, the reserve's size.
        row(fid, sprintf('$\\tau_{rec_1}$ (%s)', r.label), 'Chained STD: releasable-pool refill', r.tau_rec(1), 's');
        row(fid, sprintf('$\\tau_{rec_2}$ (%s)', r.label), 'Chained STD: reserve refill', r.tau_rec(2), 's');
        row(fid, sprintf('$\\tau_{rel}$ (%s)', r.label),   'Chained STD: release (releasable pool only)', r.tau_rel, 's');
        row(fid, sprintf('$n_2$ (%s)', r.label), ...
            'Chained STD: reserve capacity / releasable-pool capacity', r.n2, '--');
    elseif ~isempty(r.tau_rec)
        row(fid, sprintf('$\\tau_{rec}$ (%s)', r.label), 'STD recovery', r.tau_rec, 's');
        row(fid, sprintf('$\\tau_{rel}$ (%s)', r.label), 'STD release',  r.tau_rel, 's');
    end
    if r.has_stf
        row(fid, sprintf('$\\tau_{dec}$ (%s)', r.label), 'STF decay',    r.tau_dec, 's');
        row(fid, sprintf('$\\tau_{fac}$ (%s)', r.label), 'STF rate',     r.tau_fac, 's');
        row(fid, sprintf('$G$ (%s)', r.label),           'STF ceiling',  r.G, '--');
    end
end
if M.has_noise
    row(fid, '$\sigma_u$', 'Input-referred noise',    M.sigma_u_noise, 'input units');
    row(fid, '$\sigma_x$', 'Stationary std of $x$',   M.sigma_u_noise/sqrt(2*M.tau_d), '--');
end
row(fid, '$f_s$',        'Sample rate',               M.fs,       'Hz');
fprintf(fid, '\n');
end

%% ------------------------------------------------------------------------
function p = write_conditions_table(dir_out, preset_name, conditions, M)
p = fullfile(dir_out, 'adaptation_conditions.md');
fid = open_md(p, 'Adaptation conditions', preset_name, M);
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, ['Every sweep runs each grid point under all %d regimes, on the ' ...
    '**same network** (the weight seed is shared), so the comparison is ' ...
    'paired.\n\n'], M.n_all_conditions);
if numel(conditions) < M.n_all_conditions
    fprintf(fid, 'This table covers the %d of them this figure set draws.\n\n', numel(conditions));
end
write_conditions_rows(fid, conditions, M);
end

function p = write_conditions_table_bare(dir_out, conditions, M)
% The same table and nothing else -- no title, no provenance, no prose -- for
% dropping straight into a manuscript. Shares write_conditions_rows with the
% annotated file, so the two cannot disagree.
p = fullfile(dir_out, 'adaptation_conditions_table.md');
fid = fopen(p, 'w');
if fid < 0; error('write_manuscript_tables:CannotOpen', 'Could not open %s', p); end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
write_conditions_rows(fid, conditions, M);
end

function write_conditions_rows(fid, conditions, M)
% Header row and one row per regime; the only place the table is spelled out.
if M.is_pairs
    % The Depression column (FD-006) names each condition's depression model
    % and its time constants, read off the compiled routes (depression_facts):
    % without it the chained and the product MTS print identical rows (same
    % routes, same tau_a).
    fprintf(fid, ['| Condition | SFA timescales $\\tau_a$ (s) | Depressing routes ' ...
        '| Depression | Facilitating routes |\n']);
    fprintf(fid, '|---|---|---|---|---|\n');
    for k = 1:numel(conditions)
        cd_ = conditions{k};
        R = cond_routes(cd_);
        % The TIMESCALES, not a count. Conditions carry tau_a and n_a is derived
        % from it on the model, so there is no n_a field here to read -- and the
        % values are what a reader actually wants: two regimes can share a count
        % and differ in which timescales they use.
        fprintf(fid, '| %s | %s | %s | %s | %s |\n', ...
            tidy(cd_.name), tau_a_str(cd_), R.std, M.cond_dep{k}, R.stf);
    end
    fprintf(fid, '\n');
else
    fprintf(fid, '| Condition | $n_{a_E}$ | $n_{b_E}$ |\n|---|---|---|\n');
    for k = 1:numel(conditions)
        cd_ = conditions{k};
        fprintf(fid, '| %s | %d | %d |\n', tidy(cd_.name), cd_.n_a_E, cd_.n_b_E);
    end
    fprintf(fid, '\n');
end
end

function s = tau_a_str(cd_)
% A condition's SFA timescales, as a readable list. Conditions carry tau_a as a
% 1 x C cell; SFA is on the first cell type, and an empty row means none.
if ~isfield(cd_, 'tau_a') || isempty(cd_.tau_a) || isempty(cd_.tau_a{1})
    s = 'none';
    return
end
t = cd_.tau_a{1};
s = strjoin(arrayfun(@(x) sprintf('%.4g', x), t, 'UniformOutput', false), ', ');
end

function R = cond_routes(cd_)
R = struct('std', 'none', 'stf', 'none');
if ~isfield(cd_, 'synapse_config') || isempty(fieldnames(cd_.synapse_config))
    return
end
sc = cd_.synapse_config;
s_list = {}; f_list = {};
pres = fieldnames(sc);
for a = 1:numel(pres)
    posts = fieldnames(sc.(pres{a}));
    for b = 1:numel(posts)
        e = sc.(pres{a}).(posts{b});
        % Routes read as "EE", "EI", ... here -- presynaptic type first, as in
        % the mu_EE_relative sweep aliases. The parameter table keeps the longer
        % "E to E" form, where a route label sits next to a tau symbol.
        lbl = sprintf('%s%s', pres{a}, posts{b});
        if isfield(e, 'std') && ~isempty(e.std); s_list{end+1} = lbl; end %#ok<AGROW>
        if isfield(e, 'stf') && ~isempty(e.stf); f_list{end+1} = lbl; end %#ok<AGROW>
    end
end
if ~isempty(s_list)
    % Routes only. NOTE this cannot tell std_only from std_only_oneTS: both
    % depress the same four routes and differ only in how many timescales each
    % carries, which used to be spelled out here as a tau_rec suffix. The
    % parameter table still names those timescales.
    R.std = strjoin(s_list, ', ');
end
if ~isempty(f_list); R.stf = strjoin(f_list, ', '); end
end

%% ------------------------------------------------------------------------
function fid = open_md(path, title, preset_name, M)
fid = fopen(path, 'w');
if fid < 0
    error('write_manuscript_tables:CannotOpen', 'Could not open %s', path);
end
fprintf(fid, '# %s\n\n', title);
fprintf(fid, '**Preset:** `%s`  ·  **Model class:** `%s`\n\n', preset_name, M.model_class);
fprintf(fid, ['> GENERATED by `write_manuscript_tables` from a model BUILT from ' ...
    'that preset under the `%s` condition. Do not edit by hand -- ' ...
    'every value here is read off the object, so it cannot disagree with what ' ...
    'the sweeps ran. Regenerate after changing the preset.\n\n'], M.condition);
fprintf(fid, '_Generated %s._\n\n', char(datetime('now')));
end

function row(fid, sym, name, val, units)
fprintf(fid, '| %s | %s | %s | %s |\n', sym, name, fmt(val), units);
end

function s = fmt(v)
if ischar(v);          s = ['`' v '`'];
elseif iscell(v);      s = ['`' strjoin(cellfun(@(x) mat2str(x, 4), v, 'UniformOutput', false), '; ') '`'];
elseif islogical(v);   s = ['`' mat2str(v) '`'];
elseif isempty(v);     s = '`[]`';
elseif isscalar(v);    s = sprintf('%g', v);
else;                  s = ['`' mat2str(v, 4) '`'];
end
end

function s = tidy(name)
% Condition display name. Uses the shared title map so the generated tables name
% a regime exactly as every figure does; the regexp path below is the fallback
% for a name the map has not been taught (an archived run may carry anything),
% and it produced things like "Sfa3 std1" for the newer regimes.
titles = srnn_condition_titles();
if titles.isKey(name)
    s = strrep(titles(name), '\tau', 'tau');   % plain text, not tex, in markdown
    return
end
s = strrep(name, '_', ' ');
s = regexprep(s, '\<sfa\>', 'SFA');
s = regexprep(s, '\<std\>', 'STD');
s = [upper(s(1)) s(2:end)];
end

function v = get_or(obj, name, default)
if isprop(obj, name); v = obj.(name); else; v = default; end
end

function d = ensure_dir(d)
if ~isfolder(d); mkdir(d); end
end
