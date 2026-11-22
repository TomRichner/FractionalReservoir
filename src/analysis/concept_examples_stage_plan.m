function plan = concept_examples_stage_plan(cfg)
% CONCEPT_EXAMPLES_STAGE_PLAN The networks a bundle's concept_examples stage builds and runs, validated before any compute.
%
%   plan = CONCEPT_EXAMPLES_STAGE_PLAN(cfg)
%
% The pre-flight half of the opt-in stage 'concept_examples' (Plan 2,
% package B; Tom, 2026-11-21: "We need the reader to really understand the
% concepts visually with figure 1 before they read the rest of the
% paper"). Config in, plan out: nothing is computed and no file is
% touched. run_concept_examples_stage carries the plan out and
% cx_load reads the result back.
%
% WHAT IT IS FOR. Figure 1A (fig_main1_concept_networks) and the
% connectivity gallery S3 (fig_supp3_connectivity_gallery). Both draw
% networks of the paper's own sparse E/I ensemble (the bundle's preset,
% 500 neurons, indegree 100) WITHOUT ADAPTATION (the preset's
% 'no_adaptation' condition), and both draw the eigenvalues of the
% Jacobian WITH EVERY GAIN EQUAL: the piecewise nonlinearity has slope 1
% on its linear part, so at gain 1
%     J = (-I + W) / tau_d,   eigenvalues (lambda(W) - 1) / tau_d   in 1/s,
% the bulk a disk of radius R / tau_d about -1/tau_d and the imaginary axis
% the stability boundary. R is the structural bulk radius of the drawn
% matrix (ospec_weight_stats' R_struct, from the realized block variances).
%
% TWO LISTS
%   examples   the four columns of Figure 1A. Each is BUILT, its spectrum
%              and weight statistics computed, AND SIMULATED for a short
%              stretch without adaptation (the preset's noise and
%              nonlinearity), keeping a handful of E and I neurons, the
%              population mean rates and lambda_1 (Benettin).
%   gallery    the structures of S3. Built and their spectra computed;
%              NOT simulated.
% Each network of both lists also gets a representative weight image: the
% submatrix of W on n_image neurons (half E, half I, drawn from a private
% stream), so its connection probability and block statistics are those of
% the network itself (weight_image_sample).
%
% A NETWORK is a row {key, name, level_of_chaos, weight_seed}:
%   key            what the figures and the table call it
%   name           ospec_network's grammar ('reference', 'EEx1p5',
%                  'lrA2p5W1p0', 'EEx1p6-EIx1p33-IEx1p33', ...; no noise
%                  token)
%   level_of_chaos the class's gain on every weight (1: the preset's
%                  network). Recurrence strength: the bulk radius and the
%                  block-mean outliers scale with it.
%   weight_seed    [] : the stage's seed (plan.seed). [a b] : this draw
%                  (a network chosen for its draw: the local outlier).
%
% THE DEFAULT EXAMPLES (found at 500 neurons on 2026-11-21; the numbers in
% the comments are seed [7 8], no adaptation)
%   decay        reference, level 0.15   bulk inside the boundary (rightmost
%                -2.7 1/s), activity settles to a fixed point, lambda_1 -3.7
%   irregular    reference, level 1      the paper's network: the bulk
%                crosses (rightmost +39 1/s), activity irregular, lambda_1
%                +2.1 to +2.3
%   runaway      EEx1p5, level 1         one real eigenvalue at +252 1/s;
%                every E neuron at the ceiling within 0.3 s, lambda_1 -5.5
%   oscillation  lrA2p5W1p0, level 1     an added rotating rank-two
%                structure: a pair at 111 +- 48i 1/s; activity oscillates
%                at about 0.6 to 0.9 Hz, lambda_1 about 0
% Column 4 is NOT the strengthened E-I loop of the figure plan. Without
% adaptation every E-I loop tried (EEx1p6-EIx1p33-IEx1p33 at levels 0.15
% and 1, and a grid of E-to-E 1.6 to 3, I-to-E 0.7 to 1.5, E-to-I 1.33 to
% 2, I-to-I 0.5 and 1 at level 0.15) pins: E silent with I tonic, or E
% and I at the ceiling; at level 1 the loop is irregular, not periodic
% (lambda_1 +4.9). The two-population mean field of this nonlinearity
% (setpoint 0.2, so a neuron at zero input fires at 0.3) gave no limit
% cycle anywhere on a wider grid either: the E-silent state is stable
% whenever the I-to-E weight times the tonic I rate exceeds the distance
% to E's floor, and avoiding it and the saturated corner leaves no room
% for an unstable focus. The rotating low-rank addition oscillates because
% its patterns are mean-free over the neurons: no population is pushed to
% the floor or the ceiling.
%
% THE DEFAULT GALLERY (S3; all at level 1, seed [7 8] unless stated):
% balanced (reference); E-to-E up (EEx1p5); I-to-E down (EIx0p9); I-to-I
% up (IIx1p25, a real eigenvalue at +93 1/s); E-I loop
% (EEx1p6-EIx1p33-IEx1p33); local outlier (reference on weight seed
% [27 28], a real eigenvalue at +54 1/s); added real low-rank mode
% (lrA1p5W0p0); added rotating low-rank mode (lrA1p4W1p0).
% EVERY ONE IS A NETWORK OF THE outlier_structures STAGE
% (outlier_structures_stage_plan's default list, same name and weight
% seed), eight that span its three groups: the block means (three single
% routes and the E-I loop), the local outliers (one, and the reference
% network they are compared with) and the added low-rank modes (one real,
% one rotating). So what S3 shows is what S7 panel B tests
% (fig_supp7_grouped draws the gallery's networks there). Until 2026-11-22
% the third structure was 'E-to-E down' (EEx0p9, a complex pair to the left
% of the boundary at -21 +- 81i 1/s), which that stage does not run: S3
% showed a structure S7 said nothing about. A gallery given through
% cfg.concept_examples_gallery may still name any network.
%
% THE CONFIG
%   cfg.preset_name   the bundle's two-type SRNNCellTypePairs preset; it
%                     must state the condition 'no_adaptation'
%   cfg.fs, cfg.ode_solver   the bundle's (stage_fs, ospec_bundle_solver)
%   cfg.verbose
%   cfg.concept_examples_<knob>   one field per knob; an unknown one is refused
%
% THE KNOBS (cfg.concept_examples_<name>; default)
%   examples     the m x 4 cell above (m >= 1)
%   gallery      the m x 4 cell above (may be empty: {})
%   seed         [7 8]   the weight seed of rows with weight_seed []
%   T            12      seconds simulated per example, from the class's
%                        initial state
%   dt_sample    0.01    seconds between kept samples
%   n_trace      5       neurons kept per type (the first n_trace of each
%                        type: the order of the neurons is random already)
%   lya          true    lambda_1 by the class's Benettin estimate over
%                        [lya_from, T]
%   lya_from     4       seconds; the warmup before it is min(5, lya_from)
%   n_image      30      neurons of the weight image (even: half E, half I)
%   image_seed   1       the private stream the image's neurons are drawn from
%   fs           the bundle's cfg.fs, else 500 Hz
%   n_override   0       the network's size for tests (0: the preset's;
%                        else at least 20, with indegree 0.2 n and the
%                        smaller network's own weight unit, as the
%                        outlier stages do)
%
% THE PLAN
%   plan.folder     'concept_examples'
%   plan.file       'concept_examples.mat'
%   plan.examples, plan.gallery   1 x m struct: key, name, label,
%                   level_of_chaos, seed (1 x 2), target
%   plan.<knob>     every knob, resolved; plan.preset_name, plan.condition
%                   ('no_adaptation'), plan.ode_solver, plan.verbose
%
% See also: run_concept_examples_stage, cx_load, cx_network, weight_image_sample,
%           fig_main1_concept_networks, fig_supp3_connectivity_gallery,
%           ospec_network, outlier_structures_stage_plan

arguments
    cfg (1,1) struct
end

PREFIX = 'concept_examples_';
KNOBS = {'examples', 'gallery', 'seed', 'T', 'dt_sample', 'n_trace', 'lya', 'lya_from', 'n_image', ...
    'image_seed', 'fs', 'n_override'};
CONDITION = 'no_adaptation';

preset_name = field_or_empty(cfg, 'preset_name');
if isempty(preset_name)
    bad_config('The concept_examples stage runs on cfg.preset_name, which this bundle does not set.');
end
names = fieldnames(cfg);
names = names(startsWith(names, PREFIX));
unknown = setdiff(names, strcat(PREFIX, KNOBS));
if ~isempty(unknown)
    bad_config('Unknown setting(s) %s. The knobs of the concept_examples stage are %s%s.', ...
        strjoin(unknown', ', '), PREFIX, strjoin(KNOBS, [', ' PREFIX]));
end
knob = @(name, default) override(cfg, [PREFIX name], default);

%% The preset
[preset, model_class, conditions] = srnn_param_preset(preset_name);   % errors on an unknown name
if ~strcmp(model_class, 'SRNNCellTypePairs')
    bad_config('The concept_examples stage runs SRNNCellTypePairs presets; ''%s'' is written for %s.', ...
        preset_name, model_class);
end
if ~(isfield(preset, 'mu_tilde_relative') && isequal(size(preset.mu_tilde_relative), [2 2]))
    bad_config(['The network names of the concept_examples stage scale the block means of a two-type ' ...
        'network; preset ''%s'' does not state a 2 x 2 mu_tilde_relative.'], preset_name);
end
cond_names = cellfun(@(c) c.name, conditions, 'UniformOutput', false);
if ~ismember(CONDITION, cond_names)
    bad_config('The concept_examples stage runs the network without adaptation; preset ''%s'' has no condition ''%s''.', ...
        preset_name, CONDITION);
end

%% The knobs
P = struct();
P.seed       = knob('seed', [7 8]);
P.T          = knob('T', 12);
P.dt_sample  = knob('dt_sample', 0.01);
P.n_trace    = knob('n_trace', 5);
P.lya        = knob('lya', true);
P.lya_from   = knob('lya_from', 4);
P.n_image    = knob('n_image', 30);
P.image_seed = knob('image_seed', 1);
P.fs         = knob('fs', stage_fs(field_or_empty(cfg, 'fs'), 500, 'concept_examples_stage_plan'));
P.n_override = knob('n_override', 0);

if ~(isnumeric(P.seed) && isequal(size(P.seed), [1 2]) && all(P.seed == round(P.seed)) && all(P.seed >= 0))
    bad_config('cfg.%sseed must be a 1 x 2 pair of whole numbers (rng_seeds).', PREFIX);
end
for k = {'T', 'dt_sample', 'fs', 'lya_from'}
    v = P.(k{1});
    if ~(is_real_scalar(v) && v > 0)
        bad_config('cfg.%s%s must be a positive scalar.', PREFIX, k{1});
    end
end
if ~(islogical(P.lya) && isscalar(P.lya)) && ~(isnumeric(P.lya) && isscalar(P.lya) && ismember(P.lya, [0 1]))
    bad_config('cfg.%slya must be true or false.', PREFIX);
end
P.lya = logical(P.lya);
if P.lya && P.lya_from >= P.T
    bad_config('cfg.%slya_from (%g s) must be before the end of the run, T = %g s.', PREFIX, P.lya_from, P.T);
end
if ~(is_whole(P.n_trace) && P.n_trace >= 1)
    bad_config('cfg.%sn_trace must be a whole number of at least 1.', PREFIX);
end
if ~(is_whole(P.n_image) && P.n_image >= 2 && mod(P.n_image, 2) == 0)
    bad_config('cfg.%sn_image must be an even whole number of at least 2 (half E, half I).', PREFIX);
end
if ~is_whole(P.image_seed) || P.image_seed < 0
    bad_config('cfg.%simage_seed must be a whole number of at least 0.', PREFIX);
end
if ~(is_whole(P.n_override) && (P.n_override == 0 || P.n_override >= 20))
    bad_config('cfg.%sn_override must be 0 (the preset''s size) or a whole number of at least 20.', PREFIX);
end
n_net = preset.n;
if P.n_override > 0; n_net = P.n_override; end
if P.n_image > n_net || P.n_trace > n_net / 2
    bad_config('cfg.%sn_image (%d) and 2 n_trace (%d) cannot exceed the network''s %d neurons.', ...
        PREFIX, P.n_image, 2 * P.n_trace, n_net);
end
if P.dt_sample < 1 / P.fs - 1e-12 || abs(P.dt_sample * P.fs - round(P.dt_sample * P.fs)) > 1e-9
    bad_config('cfg.%sdt_sample (%g s) must be a whole number of integration steps 1/fs = %g s.', PREFIX, P.dt_sample, 1 / P.fs);
end
if P.T * P.fs < 3
    bad_config('cfg.%sT must cover at least three integration steps.', PREFIX);
end

%% The two lists
P.examples = network_list(knob('examples', default_examples()), 'examples', false, preset, P.seed, PREFIX);
P.gallery  = network_list(knob('gallery', default_gallery()), 'gallery', true, preset, P.seed, PREFIX);

verbose = 'minimal';
if isfield(cfg, 'verbose') && ~isempty(cfg.verbose); verbose = cfg.verbose; end
[ode_solver, problem] = ospec_bundle_solver(cfg, 'concept_examples_stage_plan');
if ~isempty(problem); bad_config('%s', problem); end
plan = P;
plan.ode_solver  = ode_solver;
plan.folder      = 'concept_examples';
plan.file        = 'concept_examples.mat';
plan.preset_name = preset_name;
plan.condition   = CONDITION;
plan.verbose     = verbose;
end

%% ------------------------------------------------------------------------
function spec = default_examples()
spec = { ...
    'decay',       'reference',  0.15, []; ...
    'irregular',   'reference',  1,    []; ...
    'runaway',     'EEx1p5',     1,    []; ...
    'oscillation', 'lrA2p5W1p0', 1,    []};
end

function spec = default_gallery()
spec = { ...
    'balanced',      'reference',              1, []; ...
    'EE_up',         'EEx1p5',                 1, []; ...
    'IE_down',       'EIx0p9',                 1, []; ...
    'II_up',         'IIx1p25',                1, []; ...
    'EI_loop',       'EEx1p6-EIx1p33-IEx1p33', 1, []; ...
    'local_outlier', 'reference',              1, [27 28]; ...
    'lowrank_real',  'lrA1p5W0p0',             1, []; ...
    'lowrank_rot',   'lrA1p4W1p0',             1, []};
end

function L = network_list(spec, which, may_be_empty, preset, seed, PREFIX)
L = struct('key', {}, 'name', {}, 'label', {}, 'level_of_chaos', {}, 'seed', {}, 'target', {});
if may_be_empty && (isempty(spec) || (iscell(spec) && size(spec, 1) == 0))
    return
end
if ~(iscell(spec) && ismatrix(spec) && size(spec, 2) == 4 && size(spec, 1) >= 1)
    bad_config('cfg.%s%s must be an m x 4 cell, one row {key, name, level_of_chaos, weight_seed} per network.', PREFIX, which);
end
for k = 1:size(spec, 1)
    [key, name, loc, ws] = spec{k, :};
    if ~((ischar(key) || (isstring(key) && isscalar(key))) && strlength(key) > 0 && isvarname(char(key)))
        bad_config('cfg.%s%s, row %d: the key must be a valid MATLAB name.', PREFIX, which, k);
    end
    if ~((ischar(name) || (isstring(name) && isscalar(name))) && strlength(name) > 0)
        bad_config('cfg.%s%s, row %d: the name must be text.', PREFIX, which, k);
    end
    name = char(name);
    if any(~cellfun(@isempty, regexp(strsplit(name, '-'), '^noise\d+$', 'once')))
        bad_config('cfg.%s%s, row %d (%s): no noise token in a name.', PREFIX, which, k, name);
    end
    try
        net = ospec_network(name, preset.mu_tilde_relative);
    catch err
        bad_config('cfg.%s%s, row %d: %s', PREFIX, which, k, err.message);
    end
    if ~(is_real_scalar(loc) && loc > 0)
        bad_config('cfg.%s%s, row %d (%s): level_of_chaos must be a positive scalar.', PREFIX, which, k, name);
    end
    if isempty(ws)
        ws = seed;
    elseif ~(isnumeric(ws) && isequal(size(ws), [1 2]) && all(ws == round(ws)) && all(ws >= 0))
        bad_config('cfg.%s%s, row %d (%s): the weight seed is [] or a 1 x 2 pair of whole numbers.', PREFIX, which, k, name);
    end
    L(end + 1) = struct('key', char(key), 'name', name, 'label', net.label, 'level_of_chaos', loc, ...
        'seed', double(ws), 'target', net.target); %#ok<AGROW>
end
if numel(unique({L.key})) ~= numel(L)
    bad_config('cfg.%s%s names a key twice.', PREFIX, which);
end
end

function v = override(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
if isfield(s, name) && iscell(s.(name)) && isempty(s.(name)); v = s.(name); end   % an empty gallery is a choice
end

function tf = is_real_scalar(v)
tf = isscalar(v) && isnumeric(v) && isreal(v) && isfinite(v);
end

function tf = is_whole(v)
tf = is_real_scalar(v) && v == round(v);
end

function v = field_or_empty(s, name)
v = '';
if isfield(s, name); v = s.(name); end
end

function bad_config(varargin)
error('concept_examples_stage_plan:BadConfig', varargin{:});
end
