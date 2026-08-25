function paths = write_manuscript_tables(cfg)
% WRITE_MANUSCRIPT_TABLES Generate the manuscript's equation and parameter tables.
%
%   paths = WRITE_MANUSCRIPT_TABLES()
%   paths = WRITE_MANUSCRIPT_TABLES('preset_name', p)
%
% Writes, into the Stability_Manuscript tree:
%   doc_equations_table/equation_table.md         equations + parameter table
%   doc_equations_table/adaptation_conditions.md  the four regimes, as they run
%   fig_equations/parameter_table.md              the same parameter table
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
    cfg.out_root    (1,:) char = ''
    cfg.verbose     (1,1) logical = true
end

setup_paths();
if isempty(cfg.out_root)
    % _common/ sits directly under the manuscript root.
    out_root = fileparts(fileparts(mfilename('fullpath')));
else
    out_root = cfg.out_root;
end

[preset, model_class, conditions] = srnn_param_preset(cfg.preset_name);
model = build_from_preset(cfg.preset_name, 'sfa_and_std');
M = model_facts(model, model_class, preset);

paths = {};
paths{end+1} = write_equation_table( ...
    ensure_dir(fullfile(out_root, 'doc_equations_table')), cfg.preset_name, M);
paths{end+1} = write_conditions_table( ...
    ensure_dir(fullfile(out_root, 'doc_equations_table')), cfg.preset_name, conditions, M);
paths{end+1} = write_parameter_table( ...
    ensure_dir(fullfile(out_root, 'fig_equations')), cfg.preset_name, M);

if cfg.verbose
    fprintf('[manuscript tables] preset %s (%s)\n', cfg.preset_name, model_class);
    fprintf('  %s\n', paths{:});
end
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
M.any_stf = false;
for k = 1:numel(M.routes)
    M.n_b_max = max(M.n_b_max, numel(M.routes(k).tau_rec));
    M.any_stf = M.any_stf || M.routes(k).has_stf;
end
M.has_noise = any(M.sigma_u_noise(:) > 0);
M.preset_fields = fieldnames(preset);
end

function R = route_facts(model)
% One entry per depressing/facilitating route, from the compiled synapse config.
R = struct('label', {}, 'tau_rec', {}, 'tau_rel', {}, ...
           'has_stf', {}, 'tau_dec', {}, 'tau_fac', {}, 'G', {});
sc = model.synapse_config;
if ~isstruct(sc); return; end
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
            'tau_dec', [], 'tau_fac', [], 'G', []);
        if has_std
            r.tau_rec = e.std.tau_rec(:)';
            r.tau_rel = e.std.tau_rel(:)';
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

% The dendritic equation's FORM depends on whether the preset is stochastic.
if M.has_noise
    fprintf(fid, ['$$\n\\mathrm{d}x_i = \\frac{-x_i + u_i + ' ...
        '\\sum_j w_{ij}\\, s_j\\, r_j}{\\tau_d}\\,\\mathrm{d}t ' ...
        '+ \\frac{\\sigma_u}{\\tau_d}\\,\\mathrm{d}W_i\n$$\n\n']);
    fprintf(fid, ['Additive Wiener noise enters **only** $x$, which keeps the ' ...
        'diffusion constant: Ito and Stratonovich coincide, the Milstein term ' ...
        'vanishes, the variational equation is untouched, and the noise ' ...
        'cancels in Benettin''s trajectory difference. $\\sigma_u$ is ' ...
        'input-referred, so it is comparable to the stimulus amplitude.\n\n']);
else
    fprintf(fid, ['$$\n\\frac{\\mathrm{d}x_i}{\\mathrm{d}t} = ' ...
        '\\frac{-x_i + u_i + \\sum_j w_{ij}\\, s_j\\, r_j}{\\tau_d}\n$$\n\n']);
end

% c/K, not c: the model normalises the adaptation sum by its timescale count, so
% c is the TOTAL adaptation budget and the steady state is c*r whatever K is.
fprintf(fid, ['$$\nr_i = \\phi\\!\\left( x_i - \\frac{c}{K} ' ...
    '\\sum_{k=1}^{K} a_{ik} \\right)\n$$\n\n']);
fprintf(fid, ['The rate $r_i$ is the **pre-depression** output of the ' ...
    'nonlinearity. Depression and facilitation enter as the presynaptic factor ' ...
    '$s_j$ in the recurrent sum, so both mechanisms are driven by the raw rate ' ...
    '$r_i$, not by $s_i r_i$.\n\n']);

fprintf(fid, '$$\n\\frac{\\mathrm{d}a_{ik}}{\\mathrm{d}t} = \\frac{-a_{ik} + r_i}{\\tau_{a_k}}\n$$\n\n');

if M.n_b_max > 0
    fprintf(fid, ['$$\n\\frac{\\mathrm{d}b_{im}}{\\mathrm{d}t} = ' ...
        '\\frac{1-b_{im}}{\\tau_{\\mathrm{rec},m}} - ' ...
        '\\frac{b_{im}\\, r_i}{\\tau_{\\mathrm{rel},m}}, \\qquad m = 1,\\dots,%d\n$$\n\n'], ...
        M.n_b_max);
end
if M.any_stf
    fprintf(fid, ['$$\n\\frac{\\mathrm{d}g_{i}}{\\mathrm{d}t} = ' ...
        '\\frac{1-g_i}{\\tau_{\\mathrm{dec}}} + ' ...
        '\\frac{(G - g_i)\\, r_i}{\\tau_{\\mathrm{fac}}}\n$$\n\n']);
end

% The synaptic factor, spelled for what this preset actually has.
parts = {};
if M.n_b_max == 1; parts{end+1} = 'b_i'; end
if M.n_b_max > 1;  parts{end+1} = sprintf('\\prod_{m=1}^{%d} b_{im}', M.n_b_max); end
if M.any_stf;      parts{end+1} = 'g_i'; end
if isempty(parts); parts = {'1'}; end
fprintf(fid, '$$\ns_i = %s\n$$\n\n', strjoin(parts, ' \\, '));
if M.n_b_max > 1
    % Quote the actual depth this preset reaches. The ratio tau_rec/tau_rel is
    % the ONLY combination that sets the steady state, so it is what the numbers
    % below are computed from.
    ratios = [];
    for k = 1:numel(M.routes)
        if ~isempty(M.routes(k).tau_rec)
            ratios = M.routes(k).tau_rec ./ M.routes(k).tau_rel;
            break
        end
    end
    fprintf(fid, ['The depression timescales enter as a **product**, which is ' ...
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

fprintf(fid, ['The sparse weights are drawn element-wise and scaled by the ' ...
    'synaptic gain $g$:\n\n$$\nw_{ij} = g\\, S_{ij}\\left( \\tilde{\\sigma}_{ij}\\, ' ...
    '\\xi_{ij} + \\tilde{\\mu}_{ij} \\right), \\qquad \\xi_{ij} \\sim ' ...
    '\\mathcal{N}(0,1), \\qquad S_{ij} \\sim \\mathrm{Bernoulli}(\\alpha)\n$$\n\n']);
fprintf(fid, ['with $(\\tilde{\\mu}, \\tilde{\\sigma})$ indexed by ' ...
    '**(postsynaptic, presynaptic)** cell type, and given as multiples of ' ...
    '$F = 1/\\sqrt{n\\,\\alpha(2-\\alpha)}$.\n\n']);

write_param_rows(fid, M);
end

%% ------------------------------------------------------------------------
function p = write_parameter_table(dir_out, preset_name, M)
p = fullfile(dir_out, 'parameter_table.md');
fid = open_md(p, 'SRNN model parameter table', preset_name, M);
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
write_param_rows(fid, M);
end

function write_param_rows(fid, M)
fprintf(fid, '## Parameters as run\n\n');
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
row(fid, '$c$',          'SFA coupling per type',     M.c,        '--');
for k = 1:numel(M.routes)
    r = M.routes(k);
    if ~isempty(r.tau_rec)
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
    'paired.\n\n'], numel(conditions));

if M.is_pairs
    fprintf(fid, ['| Condition | SFA timescales $\\tau_a$ (s) | Depressing routes ' ...
        '| Facilitating routes |\n']);
    fprintf(fid, '|---|---|---|---|\n');
    for k = 1:numel(conditions)
        cd_ = conditions{k};
        R = cond_routes(cd_);
        % The TIMESCALES, not a count. Conditions carry tau_a and n_a is derived
        % from it on the model, so there is no n_a field here to read -- and the
        % values are what a reader actually wants: two regimes can share a count
        % and differ in which timescales they use.
        fprintf(fid, '| %s | %s | %s | %s |\n', ...
            tidy(cd_.name), tau_a_str(cd_), R.std, R.stf);
    end
    fprintf(fid, '\n');
    fprintf(fid, ['**Why this class.** `SRNNCellTypePairs` names each synaptic ' ...
        'route individually, so depression can be put on one route and not ' ...
        'another. `SRNNModel2` cannot express that: its depression count is ' ...
        'per **presynaptic population**, so it can say "all outgoing E ' ...
        'synapses depress" but never distinguish E to E from E to I.\n\n']);
else
    fprintf(fid, '| Condition | $n_{a_E}$ | $n_{b_E}$ |\n|---|---|---|\n');
    for k = 1:numel(conditions)
        cd_ = conditions{k};
        fprintf(fid, '| %s | %d | %d |\n', tidy(cd_.name), cd_.n_a_E, cd_.n_b_E);
    end
    fprintf(fid, '\n');
end

fprintf(fid, ['**Implementation note.** When a mechanism is switched off its ' ...
    'state variables are excluded from the state vector and from the Jacobian ' ...
    'entirely, rather than being integrated at zero. This prevents spurious ' ...
    'zero eigenvalues from disabled dynamics.\n\n']);
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
        lbl = sprintf('%s to %s', pres{a}, posts{b});
        if isfield(e, 'std') && ~isempty(e.std); s_list{end+1} = lbl; end %#ok<AGROW>
        if isfield(e, 'stf') && ~isempty(e.stf); f_list{end+1} = lbl; end %#ok<AGROW>
    end
end
if ~isempty(s_list)
    % The routes ALONE cannot tell std_only from std_only_oneTS: both depress
    % the same four routes and differ only in how many timescales each carries.
    % Naming the timescales is what makes those two rows distinguishable.
    R.std = sprintf('%s (tau_rec %s)', strjoin(s_list, ', '), std_taus(sc));
end
if ~isempty(f_list); R.stf = strjoin(f_list, ', '); end
end

function s = std_taus(sc)
% The depression recovery timescales, assuming every depressing route shares
% them -- which every preset in this repo does. Falls back to naming the routes
% separately if one ever does not, rather than quietly reporting the first.
pres = fieldnames(sc);
found = {};
for a = 1:numel(pres)
    posts = fieldnames(sc.(pres{a}));
    for b = 1:numel(posts)
        e = sc.(pres{a}).(posts{b});
        if isfield(e, 'std') && ~isempty(e.std)
            found{end+1} = mat2str(e.std.tau_rec, 4); %#ok<AGROW>
        end
    end
end
u = unique(found);
if isscalar(u); s = u{1}; else; s = strjoin(u, ' / '); end
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
    'that preset under the `sfa_and_std` condition. Do not edit by hand -- ' ...
    'every value here is read off the object, so it cannot disagree with what ' ...
    'the sweeps ran. Regenerate after changing the preset.\n\n']);
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
