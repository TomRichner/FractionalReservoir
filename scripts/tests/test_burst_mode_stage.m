%% test_burst_mode_stage.m
% The opt-in stage burst_mode and the figures drawn from it (Plan 2,
% package E): Figure 8 A to C (fig_burst_mode), supplementary S10
% (fig_supp10_grouped) and, where the clinical panels are on this machine,
% Figure 8 in the new order (ied_stim_fig_combined with burst_mode_dir).
%
% Checks:
%   1. THE PLAN: every run mode is planned (3 networks at essential, 9 from
%      fast); the first seed is [42 42]; the bundle's fs and integrator are
%      taken; refused (BadConfig): an unknown knob, an unknown run mode,
%      dt_full not a multiple of dt_sample, transient_s >= hold_dur.
%   2. THE STAGE at short size (two networks, 8 s per level, the bursting
%      network's own 50 neurons): index and files written; the frozen arm
%      restarts on the run (mismatch 0) and holds every adaptive variable
%      (drift 0); M = I + tau_d J_xx exactly and = W diag(d) to 1e-12; the
%      full-Jacobian mode is real-or-upper-half (Im >= 0) and the dense grid
%      only in the first network; bursts before the step in the first
%      network, each with a coinciding full-Jacobian crossing; the frozen
%      arm saturates (mean E rate during the step above 0.5) and does not
%      burst; the loader refuses a run directory without the stage.
%   3. THE FIGURES draw from that output: PNG and _table.md written; S10
%      with a synthetic ied_stim summary (two networks); the new Figure 8
%      when the clinical folder of ied_stim_config exists (skipped with a
%      note otherwise), from the 2026-11-17 model run of the main checkout
%      or, without it, skipped.
%
% Asserts; prints one PASS line. About a minute and a half without a pool.
% Assumes setup_paths has run and a pool allocation is set (pool_set(0)
% will do).
%
% See also: burst_mode_stage_plan, run_burst_mode_stage, fig_burst_mode,
%           fig_supp10_grouped, ied_stim_fig_combined

t_all = tic;

%% 1. The plan
base = struct('run_mode', 'essential', 'fs', 500, 'ode_solver', 'sra1', 'verbose', 'near-none');
pl = burst_mode_stage_plan(base);
assert(pl.n_seeds == 3 && isequal(pl.seeds(1, :), [42 42]) && strcmp(pl.condition, 'sfa3_std1'));
assert(pl.fs == 500 && strcmp(pl.ode_solver, 'sra1') && strcmp(pl.folder, 'burst_mode'));
assert(pl.runs(1).keep_ev && ~pl.runs(2).keep_ev);
for m = run_mode_names()
    c = base; c.run_mode = m{1};
    p = burst_mode_stage_plan(c);
    assert(p.n_seeds == 3 + 6 * ~strcmp(m{1}, 'essential'), 'mode %s', m{1});
end
refused = {setfield(base, 'burst_mode_foo', 1), setfield(base, 'run_mode', 'huge'), ...
    setfield(base, 'burst_mode_dt_full', 0.015), setfield(base, 'burst_mode_transient_s', 60)}; %#ok<SFLD>
for k = 1:numel(refused)
    assert(throws(@() burst_mode_stage_plan(refused{k}), 'burst_mode_stage_plan:BadConfig'), 'refusal %d', k);
end

%% 2. The stage at short size
cfg = base;
cfg.burst_mode_n_seeds = 2;
cfg.burst_mode_hold_dur = 8;
cfg.burst_mode_transient_s = 2;
cfg.burst_mode_frozen_pre_s = 1;
pl = burst_mode_stage_plan(cfg);
run_dir = tempname; mkdir(run_dir);
cleanup = onCleanup(@() rmdir(run_dir, 's'));
run_burst_mode_stage(cfg, run_dir, pl);
B = bm_load_runs(run_dir, '', 'all');
assert(height(B.runs) == 2 && all(B.runs.ok));
for k = 1:2
    R = B.R{k};
    assert(R.check.restart_mismatch == 0 && R.check.held_drift == 0, 'frozen arm, network %d', k);
    assert(R.check.M_vs_I_plus_tau_Jxx == 0 && R.check.M_vs_W_diag_d < 1e-12, 'M check, network %d', k);
    assert(all(imag(R.full.net) >= 0));
    assert(xor(k == 1, all(R.full.is_coarse)), 'dense full-Jacobian grid in the first network only');
    assert(isempty(R.intact.ev) == (k > 1));
    % The no-stimulation continuation (Figure 8 G): the intact run before the step.
    assert(R.check.nostim_pre_step_mismatch == 0 && R.nostim.t(1) < R.t_step && all(R.nostim.u < R.dc_level / 2), ...
        'no-stimulation continuation, network %d', k);
end
R = B.R{1};
assert(B.runs.n_bursts_before(1) >= 1, 'the first network bursts before the step');
assert(R.meas.full.before.n_bursts_with_cross == R.meas.full.before.n_bursts, 'every burst preceded by a crossing');
assert(B.runs.rate_E_during_frozen_arm(1) > 0.5 && B.runs.n_bursts_during_frozen_arm(1) == 0, 'frozen arm saturates');
assert(throws(@() bm_load_runs(tempname), ''), 'the loader refuses a run without the stage');

%% 3. The figures
out_dir = fullfile(run_dir, 'figs');
o = fig_burst_mode('run_dir', run_dir, 'out_dir', fullfile(out_dir, 'a'), 'visible', false);
% Polish of 2026-11-22: E and I in the manuscript's colors; the two modes
% told apart by color and width. Plan 3 (Figure 8 G): one key, above the
% time course, with the no-stimulation continuation dashed; the
% adaptation-frozen arm is not drawn.
st = manuscript_style();
assert(isequal([st.E_color; st.I_color], SRNNCellTypePairs.type_colors(2)), ...
    'the manuscript''s E and I colors are the model''s two anchors');
lg = findall(o.figs, 'Type', 'legend');
assert(isscalar(lg), 'one key on G');
names = lg.String; h = lg.PlotChildren; ax = lg.Axes;
assert(isequal(names, {'E', 'I', 'Gains frozen', 'Full Jacobian', 'Without stimulation'}), 'the key of G');
assert(isequal(h(1).Color, st.E_color) && isequal(h(2).Color, st.I_color), 'E and I traces in the manuscript_style colors');
assert(~isequal(h(3).Color, h(4).Color) && h(3).LineWidth >= 2 * h(4).LineWidth, ...
    'gains frozen and full Jacobian differ in color and in width');
assert(strcmp(h(5).LineStyle, '--'), 'the continuation without stimulation is dashed');
set(lg, 'Units', 'inches');
top = max(arrayfun(@(a) a.Position(2) + a.Position(4), findall(o.figs, 'Type', 'axes')'));
assert(lg.Position(2) >= top - 1e-6, 'the key sits above the time course');
% Plan 3 review (the lead's decision 16): after the step only the E rate and
% the full Jacobian's mode continue dashed, in one quiet gray.
dashed = findall(ax, 'Type', 'line', 'LineStyle', '--');
assert(numel(dashed) == 2, 'dashed in the rate axes: E without stimulation, and the key''s entry');
assert(isequal(dashed(1).Color, dashed(2).Color) && isequal(h(5).Color, dashed(1).Color), ...
    'the continuation and its key entry in one gray');
mode_ax = findall(o.figs, 'Type', 'axes');
mode_ax = mode_ax(arrayfun(@(a) startsWith(strjoin(cellstr(a.YLabel.String), ' '), 'Re'), mode_ax));
n_mode_dash = sum(arrayfun(@(a) numel(findall(a, 'Type', 'line', 'LineStyle', '--')), mode_ax));
assert(n_mode_dash == 1, 'one dashed mode (the full Jacobian''s) without stimulation, found %d', n_mode_dash);
close(o.figs);
assert(isfile(fullfile(out_dir, 'a', 'Fig_Burst_Mode.png')) && isfile(fullfile(out_dir, 'a', 'Fig_Burst_Mode_table.md')));
ied = fullfile(run_dir, 'ied'); mkdir(ied);
rows = table([42; 42; 73; 73], [1; 2; 1; 2], [false; true; false; true], [0.10; 0.15; 0.12; 0.17], ...
    [0.05; 0.08; 0.06; 0.09], [0.30; 0.37; 0.31; 0.40], ...
    'VariableNames', {'seed', 'block', 'stim', 'mean_rate', 'sfa_feedback', 'depression'});
save(fullfile(ied, 'summary.mat'), 'rows');
o = fig_supp10_grouped('run_dir', run_dir, 'ied_run_dir', ied, 'out_dir', fullfile(out_dir, 'b'), 'visible', false);
% Panel B: no trace leaves its axes (the depression trace used to dip under the floor).
axs = findall(o.figs, 'Type', 'axes');
n_dep = 0;
for k = 1:numel(axs)
    ln = findall(axs(k), 'Type', 'line');
    y = arrayfun(@(l) [min(l.YData), max(l.YData)], ln, 'UniformOutput', false);
    y = vertcat(y{:});
    if isempty(y) || ~startsWith(strjoin(cellstr(axs(k).YLabel.String), ' '), 'Depr'); continue; end
    n_dep = n_dep + 1;
    assert(axs(k).YLim(1) < min(y(:, 1)) && axs(k).YLim(2) >= max(y(:, 2)), 'S10 (B): the depression traces are inside their axes');
end
assert(n_dep == 1, 'S10 (B): one depression panel checked, found %d', n_dep);
% The common canvas (Plan 2): canvas_w wide, type by role.
st = manuscript_style();
o.figs.Units = 'inches';
assert(abs(o.figs.Position(3) - st.canvas_w) < 1e-6, 'S10 is canvas_w wide');
A = canvas_audit(o.figs);
assert(all(A.ok), 'S10: type by role; off: %s', strjoin(A.name(~A.ok), '; '));
close(o.figs);
assert(isfile(fullfile(out_dir, 'b', 'Fig_Supp10_Grouped.png')) && isfile(fullfile(out_dir, 'b', 'Fig_Supp10_Grouped_table.md')));
assert(throws(@() fig_supp10_grouped('run_dir', run_dir, 'ied_run_dir', tempname, 'save', false, 'visible', false), ''), ...
    'S10 errors without an ied_stim run');

c8 = ied_stim_config('stim_blocks_20261117');
legacy = fullfile(fileparts(fileparts(which('setup_paths'))), 'FractionalReservoir', 'data', 'IED_model', 'stim_blocks_20261117');
note = '';
if isfolder(fullfile(c8.clinical_dir, 'panels')) && isfile(fullfile(legacy, 'summary.mat'))
    mdir = fullfile(out_dir, 'model');
    fig_ied_stim_blocks('data_dir', legacy, 'out_dir', mdir, 'which', 'main');
    o = ied_stim_fig_combined('data_dir', legacy, 'out_dir', fullfile(out_dir, 'c'), ...
        'model_fig', fullfile(mdir, 'fig_ied_stim_blocks.fig'), 'clinical_dir', c8.clinical_dir, ...
        'burst_mode_dir', run_dir);
    info = imfinfo(o.files{1});
    st = manuscript_style();
    dpi = info.XResolution * 0.0254;                 % the PNG's resolution is in pixels per meter
    assert(abs(info.Width / dpi - st.canvas_w) <= st.canvas_tol && info.Height / dpi <= st.canvas_max_h + 0.01, ...
        'Figure 8 is on the common canvas: %.2f x %.2f in', info.Width / dpi, info.Height / dpi);
    assert(all(cellfun(@isfile, o.files)));
    % Titles of the clinical and modular panels: sentence case, no counts.
    f8 = openfig(o.files{2}, 'invisible');
    titles = arrayfun(@(a) string(strjoin(cellstr(a.Title.String), ' ')), findall(f8, 'Type', 'axes'));
    close(f8);
    assert(all(ismember(["No stimulation", "Stimulation", "Model, no stimulation", "Model, stimulation"], titles)) && ...
        ~any(contains(titles, ["(", "Median", "No Stimulation"])), 'Figure 8: sentence-case titles without counts');
    assert(contains(fileread(o.files{4}), 'Titles restyled on the composite'), 'the provenance lists the original titles');
else
    note = ' (new Figure 8 skipped: no clinical panels or no 2026-11-17 model run on this machine)';
end
fprintf('PASS test_burst_mode_stage%s (%.0f s)\n', note, toc(t_all));

function tf = throws(f, id)
tf = false;
try
    f();
catch err
    tf = isempty(id) || strcmp(err.identifier, id);
end
end
