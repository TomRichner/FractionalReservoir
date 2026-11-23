function F = outlier_fig_data(cfg, caller, kind, want)
% OUTLIER_FIG_DATA What the four outlier figures draw, read from ONE run's outlier_spectra stage.
%
%   F = outlier_fig_data(cfg, 'fig_outlier_plane', 'plane')
%   F = outlier_fig_data(cfg, 'fig_outlier_dose', 'dose', 'networks', 'dose')
%
% THE ONLY PLACE the outlier figures (fig_outlier_plane, _dose,
% _slow_variables, _dwell; FD-045) read the stage's files (FD-040; the
% contract is the header of ospec_load_run). A figure never touches a run's
% struct R: this function loads each run it needs, keeps the few numbers
% the figure of that `kind` draws or tabulates (the local functions
% reduce_*, which name every field of R that is used), and lets R go. So a
% renamed field of the contract is a change to this file alone.
%
% It reads one run directory, <cfg.run_dir>/outlier_spectra/ospec_index.mat,
% or the index file given as cfg.data_file, and errors when neither is
% there (resolve_data_file: no fallback to another location).
%
% INPUT
%   cfg      the figure's cfg: run_dir, data_file, conditions, seed
%   caller   the figure's name, for error identifiers
%   kind     'plane' | 'dose' | 'slow' | 'dwell'
%   networks (name-value) cellstr of network names as the index has them,
%            in the order to draw; or 'dose': the reference and every
%            E-to-E level (EEx<a>p<b>) the index holds, by level
%   x_thr    (name-value; 'slow') the dendritic fraction that makes an
%            eigenvalue of the full Jacobian a network mode (0.1)
%
% OUTPUT F
%   file, folder      the index file and the stage's folder
%   settings          the plan the stage ran
%   networks          1 x n_net cellstr; net_labels the same as axis text
%   cond_names        1 x n_cond cellstr, in the order drawn (select_conditions)
%   colors, titles    per condition (manuscript_style, sentence case)
%   seeds             n_seed x 2, in the order of the index (the stage lists
%                     the seeds outermost, so row 1 is the first seed run)
%   i_seed            the row of `seeds` that is DRAWN where one seed is
%                     drawn: cfg.seed ([] = row 1; a first seed number or
%                     the 1 x 2 pair)
%   D                 n_net x n_cond x n_seed cell: the reduced run, or []
%                     where the index has no finished run
%   window_txt        the measured window(s) for a caption, from the runs
%                     ('20 s', or '200 s (reference), 20 s (EEx1p5)': the
%                     networks need not share one window, FD-057)
%   n_failed          runs of the selection the index marks as not finished
%   st                manuscript_style()
%
% Errors (identifier <caller>:...): NoNetwork (a named network is not in
% the index), NoRuns (nothing finished in the selection), NoSeed (cfg.seed
% is not a seed of the selection), NeedsDense / NeedsFull (kind 'slow' on
% runs without the dense full-Jacobian spectra; a figure may ask).
%
% See also: ospec_read_index, ospec_load_run, ospec_measures,
%           ospec_excursions, fig_outlier_plane, fig_outlier_dose,
%           fig_outlier_slow_variables, fig_outlier_dwell

arguments
    cfg (1,1) struct
    caller (1,:) char
    kind (1,:) char {mustBeMember(kind, {'plane', 'dose', 'slow', 'dwell'})}
    want.networks = {'reference', 'EEx1p5'}
    want.x_thr (1,1) double = 0.1
end

F = struct();
F.file = resolve_data_file(cfg.data_file, cfg.run_dir, ...
    {fullfile(cfg.run_dir, 'outlier_spectra')}, 'ospec_index.mat', ...
    'Run the opt-in stage outlier_spectra (run_outlier_spectra_stage) first');
F.folder = fileparts(F.file);
[runs, F.settings] = ospec_read_index(F.folder);

%% Networks
have = unique(runs.network, 'stable');
if ischar(want.networks) && strcmp(want.networks, 'dose')
    level = cellfun(@ee_level, have);
    keep = ~isnan(level);
    [~, order] = sort(level(keep));
    nets = have(keep); nets = nets(order);
else
    nets = cellstr(want.networks);
    missing = setdiff(nets, have);
    if ~isempty(missing)
        error([caller ':NoNetwork'], 'No runs of network ''%s'' in %s; it holds: %s.', ...
            missing{1}, F.file, strjoin(have, ', '));
    end
end
nets = reshape(nets, 1, []);
if isempty(nets)
    error([caller ':NoNetwork'], 'No reference or E-to-E network in %s; it holds: %s.', F.file, strjoin(have, ', '));
end
sel = ismember(runs.network, nets);
F.n_failed = nnz(sel & ~runs.ok);
runs = runs(sel & runs.ok, :);
if isempty(runs)
    error([caller ':NoRuns'], 'No finished run of %s in %s.', strjoin(nets, ', '), F.file);
end
F.networks = nets;
F.net_labels = cellfun(@net_label, nets, 'UniformOutput', false);

%% Conditions, seeds
F.cond_names = select_conditions(unique(runs.condition, 'stable'), cfg.conditions);
F.st = manuscript_style();
F.colors = cellfun(@(c) F.st.condition_color(c), F.cond_names, 'UniformOutput', false);
F.titles = grouped_condition_label(F.cond_names, 'sentence');
runs = runs(ismember(runs.condition, F.cond_names), :);
F.seeds = unique([runs.seed1, runs.seed2], 'rows', 'stable');
F.i_seed = 1;
if ~isempty(cfg.seed)
    s = cfg.seed(:)';
    if isscalar(s)
        k = find(F.seeds(:, 1) == s, 1);
    else
        k = find(F.seeds(:, 1) == s(1) & F.seeds(:, 2) == s(2), 1);
    end
    if isempty(k)
        error([caller ':NoSeed'], 'Seed %s is not among the seeds of these runs: %s.', mat2str(s), mat2str(F.seeds));
    end
    F.i_seed = k;
end

%% The runs, reduced
n_net = numel(nets); n_cond = numel(F.cond_names); n_seed = size(F.seeds, 1);
F.D = cell(n_net, n_cond, n_seed);
for k = 1:height(runs)
    w = strcmp(nets, runs.network{k});
    c = strcmp(F.cond_names, runs.condition{k});
    s = F.seeds(:, 1) == runs.seed1(k) & F.seeds(:, 2) == runs.seed2(k);
    R = ospec_load_run(F.folder, runs(k, :));
    switch kind
        case 'plane'; d = reduce_plane(R, caller);
        case 'dose';  d = reduce_dose(R);
        case 'slow';  d = reduce_slow(R, want.x_thr, caller);
        case 'dwell'; d = reduce_dwell(R);
    end
    d.network = R.network; d.condition = R.condition; d.seed = R.seed; d.n = R.n;
    d.window_s = R.window_s; d.n_cal = numel(R.cal.t);
    F.D{w, c, s} = d;
end
F.window_txt = window_text(F.D, nets);
end

function s = window_text(D, nets)
% The measured windows of the loaded runs, for a caption: '20 s' when every
% run has one length, else '200 s (reference), 20 s (EEx1p1, ...)'. The
% networks need not share one window (FD-057: only long_networks get the
% long one), so a caption must not take the first run's for all.
win = nan(1, numel(nets));
for w = 1:numel(nets)
    r = D(w, :, :); r = r(~cellfun(@isempty, r));
    if ~isempty(r); win(w) = max(cellfun(@(d) d.window_s, r)); end
end
u = unique(win(~isnan(win)), 'stable');
if isscalar(u)
    s = sprintf('%g s', u);
    return
end
parts = arrayfun(@(v) sprintf('%g s (%s)', v, strjoin(nets(win == v), ', ')), sort(u, 'descend'), 'UniformOutput', false);
s = strjoin(parts, ', ');
end

%% ------------------------------------------------------------------------
function d = reduce_plane(R, caller)
% The frozen-gain spectra at the calibration states as RATES (1/s), each
% state's own theoretical radius, the null spectra, the uniform-gain
% rightmost eigenvalue, and the neurons' operating points.
tau = R.tau_d;
if ~isequal(R.cell_type_names(:)', {'E', 'I'})
    error([caller ':Types'], 'The figure colors E and I; the run has types %s.', strjoin(R.cell_type_names, ', '));
end
d.tau_d = tau;
d.center = -1 / tau;
d.z = (R.cal.ev_M - 1) ./ tau;                          % n x n_cal
d.rad = R.cal.R_eff ./ tau;                             % 1 x n_cal: each state's own theoretical radius
d.kappa = R.kappa;                                      % the null-calibrated reach, as a factor on the radius
d.z_uniform = ([R.null.uniform_re] - 1) ./ tau;         % 1 x n_cal: rightmost real part with uniform gains
nul = cat(3, R.null.ev_null);                           % n x n_draws x n_cal
d.z_null = (nul - 1) ./ tau;
d.ratio_theory_med = median(R.ratio_theory, 'omitnan'); % over the whole window
d.t = R.cal.t;
d.u = R.cal.u; d.rate = R.cal.rate; d.slope = R.cal.slope;
d.type = R.cal.cell_type(:);
d.S_a = R.S_a;
end

function d = reduce_dose(R)
% One row of measures (ospec_measures), the ones the dose figure plots or
% tabulates.
m = ospec_measures(R);
r = m.row;
keep = {'structural_ratio', 'ratio_theory_med', 'ratio_theory_p05', 'ratio_theory_p95', 'rho_uniform_theory', ...
    'reined_theory', 'frac_on_E', 'frac_on_I', 'gain_ratio_EI', 'kappa', 'lambda_1', 'rate_top_med', 'frac_unstable'};
for k = 1:numel(keep)
    d.(keep{k}) = r.(keep{k});
end
end

function d = reduce_slow(R, x_thr, caller)
% The dense full-Jacobian spectra at the calibration states (1/s), the
% rightmost frozen-gain eigenvalue as a rate, and its continuation.
if ~R.has_dense || isempty(R.dense)
    error([caller ':NeedsDense'], ['The run %s, %s, seed [%d %d] has no dense full-Jacobian spectra: the stage ' ...
        'computes them for its full networks only (cfg.outlier_spectra_full_networks).'], ...
        R.network, R.condition, R.seed(1), R.seed(2));
end
D = R.dense;
lam = double([D.lam]); xf = double([D.xfrac]); ov = double([D.ov]);   % N x n_cal
d.tau_d = R.tau_d;
d.x_thr = x_thr;
d.net = lam; d.net(xf < x_thr) = NaN;                   % the network modes; NaN the single-neuron slow ones
d.n_net = sum(xf >= x_thr, 1);
d.frozen = ([D.lam1_M] - 1) ./ R.tau_d;                 % 1 x n_cal complex
d.cont = [D.cont]; d.cont_ov = [D.cont_ov]; d.cont_xfrac = [D.cont_xfrac];
d.n_pos = [D.n_pos]; d.abscissa = [D.abscissa];
% The continuation under another dendritic fraction: is it the same eigenvalue?
thr_alt = [0.02 0.3];
d.thr_alt = thr_alt;
d.same_alt = false(numel(thr_alt), size(lam, 2));
for q = 1:size(lam, 2)
    for j = 1:numel(thr_alt)
        o = ov(:, q); o(xf(:, q) < thr_alt(j)) = -Inf;
        [~, i] = max(o);
        z = complex(real(lam(i, q)), abs(imag(lam(i, q))));
        d.same_alt(j, q) = abs(z - d.cont(q)) <= 1e-4 * max(1, abs(d.cont(q)));
    end
end
end

function d = reduce_dwell(R)
% The excursions of the rightmost eigenvalue beyond zero: frozen gains
% (the series at R.dt) and the full Jacobian's rightmost network mode (at
% R.full.dt; an eigenvalue beyond zero is a network mode, ospec_load_run).
d.frozen = ospec_excursions(R.t, R.rate_top);
d.has_full = R.has_full_series && ~isempty(R.full);
d.full = [];
if d.has_full
    d.full = ospec_excursions(R.full.t, real(R.full.net));
end
end

function v = ee_level(name)
% The E-to-E multiple of a network name: 1 for 'reference', a.b for
% 'EEx<a>p<b>', NaN for anything else (another structure, a variant).
if strcmp(name, 'reference')
    v = 1;
    return
end
tok = regexp(name, '^EEx(\d+)p(\d+)$', 'tokens', 'once');
if isempty(tok); v = NaN; else; v = str2double([tok{1} '.' tok{2}]); end
end

function s = net_label(name)
tok = regexp(name, '^EEx(\d+)p(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    s = sprintf('E-to-E %s %s.%s', char(215), tok{1}, tok{2});   % the sign itself, as Figure 5 and S3 write it (TeX's times set it with wide spaces)
elseif strcmp(name, 'reference')
    s = 'Reference';
else
    s = strrep(name, '_', '\_');
end
end
