function out = fig_stim_engages_adaptation(cfg)
% FIG_STIM_ENGAGES_ADAPTATION Bursting: a DC step engages SFA and STD.
%
%   out = FIG_STIM_ENGAGES_ADAPTATION()
%   out = FIG_STIM_ENGAGES_ADAPTATION('rng_seeds', [42 42])
%
% A small, sparse network driven by a DC staircase. Below the step the network
% is quiet; above it, adaptation and depression are engaged and the population
% BURSTS. Two figures:
%
%   1. the model's own time-series summary (stimulus, x, r, synaptic output,
%      and the adaptation/depression states)
%   2. the power spectrum of the mean dendritic potential, one trace per DC
%      level, showing how the drive reshapes the burst spectrum
%
% and, optionally, a third: population synchrony (chi^2) per DC level.
%
% PORTED to SRNNCellTypePairs behind the 'bursting_pairs' preset. Two changes
% beyond the class swap:
%
%   * THE NETWORK IS NOW A PRESET. n = 50, indegree = 10, f = [0.7 0.3],
%     piecewise S_a = 0.9 / S_c = 0.5 and the STD routes were ~120 lines of
%     loose assignments at the top of a script. They are physics, so they belong
%     in srnn_param_preset, where the generated parameter tables can see them.
%   * TRACES COME FROM plot_data, NOT FROM S_out. The original indexed the raw
%     state vector by hand --
%         x_cols = (nE*n_a_E + nI*n_a_I + nE*n_b_E + nI*n_b_I) + (1:n)
%     -- which is SRNNModel2's packing and simply does not hold for
%     SRNNCellTypePairs (whose b-states are per ROUTE, not per population). It
%     also silently breaks whenever a condition changes an adaptation count.
%     plot_data is keyed by cell type and route, so the same code works for
%     either class and cannot go out of step with the packing. plot_deci = 1
%     keeps it at full sample rate, which the PSD needs.
%
% See also: srnn_param_preset, dc_staircase_stimulus, build_from_preset

arguments
    cfg.preset_name  (1,:) char    = 'bursting_pairs'
    cfg.out_dir      (1,:) char    = ''
    cfg.rng_seeds    (1,2) double  = [42 42]
    cfg.dc_levels    (1,:) double  = [0.0 0.4]
    cfg.hold_dur     (1,1) double  = 50
    cfg.ramp_dur     (1,1) double  = 0
    cfg.noise_intensity (1,1) double = 0.001
    cfg.psd_settle   (1,1) double  = 15
    cfg.psd_win_len_s (1,1) double = 10
    cfg.fs           (1,1) double  = 400
    cfg.include_chi2 (1,1) logical = true
    cfg.save         (1,1) logical = true
    cfg.visible      (1,1) logical = true
    cfg.run_dir      (1,:) char    = ''   % unused; accepted for a uniform call
end

setup_paths();
out_dir      = default_out_dir(cfg.out_dir, mfilename('fullpath'));
project_root = fileparts(which('setup_paths'));
st           = manuscript_style();

dc_levels = cfg.dc_levels;
nL        = numel(dc_levels);
T_range   = [0, nL * cfg.hold_dur];

%% Stimulus: a DC staircase plus a white-noise probe
% The staircase sets the operating point; the noise probes how the network
% filters its input, which the per-level PSD reads out. This is stimulus
% PROTOCOL, not network physics, which is why it is here and not in the preset
% (and why the preset carries no input_config: it would need a function handle).
input_config = struct();
input_config.dc_levels       = dc_levels;
input_config.hold_dur        = cfg.hold_dur;
input_config.ramp_dur        = cfg.ramp_dur;
input_config.noise_intensity = cfg.noise_intensity;
input_config.intrinsic_drive = [];                 % required by the class
input_config.generator       = @dc_staircase_stimulus;

%% Run
model = build_from_preset(cfg.preset_name, 'sfa_and_std', ...
    'input_config',     input_config, ...
    'u_ex_scale',       1.0, ...
    'fs',               cfg.fs, ...
    'T_range',          T_range, ...
    'rng_seeds',        cfg.rng_seeds, ...
    'lya_method',       'none', ...
    'store_full_state', true, ...
    'plot_deci',        1);        % full rate: the PSD needs every sample
model.run();

[fig1, ~] = model.plot();
set(fig1, 'Name', 'bursting time series');

p = model.plot_data;
t_full = p.t(:);
x_all  = concat_types(p.x);       % n x nt, every neuron
r_all  = concat_types(p.r);
type1  = model.cell_type_names{1};
x_E    = p.x.(type1);             % the adapting population, for chi^2
r_E    = p.r.(type1);
x_mean = mean(x_all, 1)';         % nt x 1, mean dendritic potential

%% Per-level PSD of the mean dendritic potential
psd_f = logspace(log10(0.3), log10(100), 100);
cmap1 = parula(6);
cmap  = [cmap1(1,:); cmap1(5,:)];
show_levels = unique([1, nL], 'stable');
show_labels = arrayfun(@(k) level_label(k, nL), show_levels, 'UniformOutput', false);

fig2 = figure('Color', 'w', 'Name', 'bursting PSD by DC level', ...
    'Position', [200 200 620 420]);
ax2 = axes(fig2); hold(ax2, 'on');
for ii = 1:numel(show_levels)
    k   = show_levels(ii);
    sel = level_window(t_full, k, cfg.hold_dur, cfg.psd_settle);
    seg = x_mean(sel) - mean(x_mean(sel));
    win_len  = min(round(cfg.psd_win_len_s * cfg.fs), numel(seg));
    noverlap = round(0.5 * win_len);
    [pxx, fpx] = pwelch(seg, hamming(win_len), noverlap, psd_f, cfg.fs);
    loglog(ax2, fpx, pxx, 'LineWidth', st.line_lw, 'Color', cmap(ii, :), ...
        'DisplayName', show_labels{ii});
end
hold(ax2, 'off'); box(ax2, 'off');
set(ax2, 'XScale', 'log', 'YScale', 'log', 'FontSize', st.tick_fs);
xlabel(ax2, 'frequency (Hz)', 'FontSize', st.label_fs);
ylabel(ax2, 'PSD of mean x', 'FontSize', st.label_fs);
legend(ax2, 'Location', 'southwest', 'Box', 'off', 'Interpreter', 'none');

figs = [fig1, fig2];
fig_tags = {'bursting_timeseries', 'bursting_psd'};

%% Population synchrony, per DC level
chi2_report = '';
if cfg.include_chi2
    [fig3, chi2_report] = chi2_panel(t_full, x_E, r_E, dc_levels, ...
        cfg.hold_dur, cfg.psd_settle, st);
    figs(end+1) = fig3;
    fig_tags{end+1} = 'bursting_synchrony';
end

if ~cfg.visible; set(figs, 'Visible', 'off'); end

%% Save
out = struct('figs', figs, 'files', {{}}, 'source', ['preset: ' cfg.preset_name]);
if cfg.save
    for k = 1:numel(figs)
        save_figure_stable(out_dir, fig_tags{k}, figs(k));
        out.files = [out.files, existing_outputs(out_dir, fig_tags{k})];
    end
    capture_git_provenance(out_dir, project_root);

    write_figure_readme(out_dir, struct( ...
        'tag',    'fig_stim_engages_adaptation', ...
        'title',  'Stability_Manuscript figure: stimulus engages adaptation (bursting)', ...
        'script', 'fig_stim_engages_adaptation.m', ...
        'what',   ['A small, sparse network driven by a DC staircase. Below the ' ...
                   'step it is quiet; above it, adaptation and depression are ' ...
                   'engaged and the population BURSTS. Figure 1 is the model''s ' ...
                   'time-series summary; figure 2 the power spectrum of the ' ...
                   'mean dendritic potential, one trace per DC level; figure 3 ' ...
                   '(optional) population synchrony per level.'], ...
        'how',    ['Built from the bursting_pairs preset with a DC-staircase ' ...
                   'stimulus plus a white-noise probe. The staircase sets the ' ...
                   'operating point; the noise probes how the network filters ' ...
                   'its input, which the per-level PSD reads out. Each level ' ...
                   'skips its first settling seconds before the analysis window.'], ...
        'source', struct('preset', cfg.preset_name, 'rng_seeds', cfg.rng_seeds, ...
                         'dc_levels', dc_levels, 'hold_dur', cfg.hold_dur), ...
        'settings', figure_settings(model), ...
        'figures', {out.files}, ...
        'sections', struct( ...
            'heading', {'traces come from plot_data', 'synchrony'}, ...
            'body',    {readme_plotdata(), chi2_report})));
end
end

%% ------------------------------------------------------------------------
function M = concat_types(s)
% Concatenate a per-cell-type state struct into one n x nt matrix.
f = fieldnames(s);
M = s.(f{1});
for i = 2:numel(f)
    M = [M; s.(f{i})]; %#ok<AGROW>
end
end

function sel = level_window(t, k, hold_dur, settle)
% The settled part of DC level k.
lo = (k-1)*hold_dur + settle;
hi = k*hold_dur;
sel = t > lo & t <= hi;
end

function s = level_label(k, nL)
if k == 1
    s = 'no-stim';
elseif k == nL
    s = 'stim';
else
    s = sprintf('level %d', k);
end
end

function [fig3, report] = chi2_panel(t, x_E, r_E, dc_levels, hold_dur, settle, st)
% Population synchrony per DC level, on the adapting population.
%
% For a population matrix M (n x nt) with population mean R(t):
%     chi2 = var_t(R) / mean_i( var_t(M_i) )
% chi2 = 1 is perfect synchrony; chi2 = 1/n is perfect asynchrony. It is the
% mean pairwise correlation plus a 1/n floor,
%     chi2 = 1/n + (1 - 1/n)*rho_bar,
% so rho_bar is the floor-corrected, n-independent form. Both are reported.
%
% Computed on x rather than r: in the burst regime r clips at the
% nonlinearity's ceiling, and saturation distorts variance. x is unbounded and
% shows the burst cleanly. r is reported alongside for comparison.
nL = numel(dc_levels);
chi2_x = nan(1, nL); rho_x = nan(1, nL);
chi2_r = nan(1, nL); rho_r = nan(1, nL);
n = size(x_E, 1);
for k = 1:nL
    sel = level_window(t, k, hold_dur, settle);
    [chi2_x(k), rho_x(k)] = pop_chi2(x_E(:, sel));
    [chi2_r(k), rho_r(k)] = pop_chi2(r_E(:, sel));
end

fig3 = figure('Color', 'w', 'Name', 'bursting synchrony', ...
    'Position', [240 240 620 420]);
ax = axes(fig3); hold(ax, 'on');
plot(ax, dc_levels, rho_x, '-o', 'LineWidth', st.line_lw, 'DisplayName', 'x');
plot(ax, dc_levels, rho_r, '--s', 'LineWidth', st.line_lw, 'DisplayName', 'r');
hold(ax, 'off'); box(ax, 'off');
set(ax, 'FontSize', st.tick_fs);
xlabel(ax, 'DC level (input units)', 'FontSize', st.label_fs);
ylabel(ax, 'mean pairwise correlation', 'FontSize', st.label_fs);
legend(ax, 'Location', 'northwest', 'Box', 'off');

report = sprintf(['Population synchrony on the adapting population (n = %d), ' ...
    'measured over the settled part of each hold. chi2 = var_t(mean_i M) / ' ...
    'mean_i var_t(M_i): 1 is perfect synchrony, 1/n perfect asynchrony. The ' ...
    'reported rho_bar = (chi2 - 1/n)/(1 - 1/n) removes that floor. Computed on ' ...
    'x rather than r because r clips at the nonlinearity ceiling in the burst ' ...
    'regime and saturation distorts variance. Values: DC = %s, rho_bar(x) = %s, ' ...
    'rho_bar(r) = %s.'], n, mat2str(dc_levels, 3), mat2str(rho_x, 3), mat2str(rho_r, 3));
end

function [chi2, rho] = pop_chi2(M)
var_floor = 1e-24;
R  = mean(M, 1);
num = var(R);
den = mean(var(M, 0, 2));
if den < var_floor
    chi2 = NaN; rho = NaN; return
end
chi2 = num / den;
n = size(M, 1);
rho = (chi2 - 1/n) / (1 - 1/n);
end

function s = readme_plotdata()
s = ['TRACES COME FROM model.plot_data, not from the raw state vector. The ' ...
     'original indexed S_out by hand -- x_cols = (nE*n_a_E + nI*n_a_I + ' ...
     'nE*n_b_E + nI*n_b_I) + (1:n) -- which encodes SRNNModel2''s state packing ' ...
     'and does not hold for SRNNCellTypePairs, whose b-states are per ROUTE ' ...
     'rather than per population. It also broke silently whenever a condition ' ...
     'changed an adaptation count. plot_data is keyed by cell type and route, ' ...
     'so the same code serves either class and cannot go out of step with the ' ...
     'packing. plot_deci = 1 keeps it at the full sample rate, which the PSD ' ...
     'requires.'];
end
