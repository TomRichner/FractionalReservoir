function md_path = write_run_parameters_md(run_dir)
% WRITE_RUN_PARAMETERS_MD Human-readable parameters.md for a run_all directory.
%
% Usage:
%   md_path = write_run_parameters_md(fullfile(project_root, 'data', ...
%                 'param_space', 'run_all_aug_14_26_17_25'));
%
% Writes <run_dir>/parameters.md describing what the run actually used: the
% preset and run mode, every parameter the preset set, every parameter in
% effect (preset, run mode, class default alike) as run in the SFA+STD
% condition, the four adaptation conditions, and the per-analysis sweep axes
% and timings.
%
% READS THE RUN, NOT THE SOURCE. Every value comes from the run's own saved
% artifacts -- run_manifest.mat and each sub-analysis's psa_object.mat. That
% matters because srnn_param_preset is edited over time, so a directory's named
% preset stops matching today's source; the current preset is consulted ONLY to
% diff against what was recorded, and any disagreement is flagged as drift.
%
% The preset as it ran is recoverable because every sub-script sets
%   psa.model_defaults = merge_struct(preset_defaults, cfg.model)
% and saveobj persists model_defaults. Subtracting the fields of cfg.model --
% reconstructed by re-calling analysis_run_config with the recorded run_mode --
% leaves the preset exactly as it was at run time.
%
% NEVER THROWS on odd or missing data: a missing manifest, an unreadable psa, a
% legacy run with an empty resolved_defaults, an unknown preset name or run mode
% all degrade to a line in the Caveats section, because this is called at the
% end of run_all_analyses where a failure must not obscure finished compute.
%
% See also: run_all_analyses, srnn_param_preset, analysis_run_config,
%           srnn_adaptation_conditions, ParamSpaceAnalysis2

arguments
    run_dir (1,:) char
end

if ~isfolder(run_dir)
    error('write_run_parameters_md:NoSuchDir', ...
        'Run directory does not exist: %s', run_dir);
end

caveats = {};

%% Manifest and provenance
[man, c] = read_manifest(run_dir);   caveats = [caveats, c];
prov     = read_provenance_txt(run_dir);

run_mode    = getfield_or(man, 'run_mode', '');
preset_name = getfield_or(man, 'preset_name', '');
model_class = getfield_or(man, 'model_class', '');

%% Per-analysis collection
[entries, c] = collect_analyses(run_dir, run_mode);   caveats = [caveats, c];

if isempty(entries)
    caveats{end+1} = ['No sub-analysis directory contained a `psa_object.mat`, ' ...
        'so no parameters could be recovered.'];
end

% The manifest is written before any analysis runs, so prefer it; fall back to
% the analyses themselves when it is missing.
if isempty(model_class) && ~isempty(entries)
    model_class = entries(1).psa.model_class;
    caveats{end+1} = sprintf(['Model class was taken from the saved analyses ' ...
        '(`%s`) because the manifest did not record one.'], model_class);
end

%% Preset as it ran, and drift against the current source
[rec_preset, c] = merged_recorded_preset(entries);   caveats = [caveats, c];
[cur_preset, cur_ok, c] = current_preset(preset_name);   caveats = [caveats, c];

%% Write
md_path = fullfile(run_dir, 'parameters.md');
fid = fopen(md_path, 'w', 'n', 'UTF-8');
if fid < 0
    error('write_run_parameters_md:CannotWrite', ...
        'Could not open for writing: %s', md_path);
end
closer = onCleanup(@() fclose(fid));

[~, run_name] = fileparts(strip_trailing_sep(run_dir));

write_header(fid, run_name, preset_name, run_mode, model_class);
write_provenance(fid, man, prov);
write_preset_section(fid, rec_preset, cur_preset, cur_ok, preset_name, entries);
write_analyses_section(fid, entries);
write_conditions_section(fid, entries, model_class);
caveats = [caveats, write_parameters_section(fid, entries, rec_preset, model_class)];
write_caveats(fid, caveats);

fprintf('Wrote %s\n', md_path);
end

%% ------------------------------------------------------------------------
%% Reading the run
%% ------------------------------------------------------------------------

function [man, caveats] = read_manifest(run_dir)
% READ_MANIFEST run_manifest.mat -> struct, or an empty struct plus a caveat.
caveats = {};
man = struct();
f = fullfile(run_dir, 'run_manifest.mat');
if ~isfile(f)
    caveats{end+1} = ['No `run_manifest.mat` in this directory, so the preset ' ...
        'name, run mode and git provenance are unknown.'];
    return
end
try
    s = load(f);
    if isfield(s, 'run_manifest')
        man = s.run_manifest;
    else
        caveats{end+1} = '`run_manifest.mat` contained no `run_manifest` variable.';
    end
catch ME
    caveats{end+1} = sprintf('Could not read `run_manifest.mat`: %s', ME.message);
end
end

function prov = read_provenance_txt(run_dir)
% READ_PROVENANCE_TXT The `key: value` lines of git_provenance.txt.
%
% Only the fields the manifest does NOT carry are of interest here (hostname,
% platform, matlab), but parsing the lot is simpler than picking.
prov = struct();
f = fullfile(run_dir, 'git_provenance.txt');
if ~isfile(f), return; end
try
    txt = fileread(f);
catch
    return
end
lines = regexp(txt, '\r?\n', 'split');
for i = 1:numel(lines)
    % "hostname: R5611351 | platform: PCWIN64 | matlab: 2026a" -- one line can
    % carry several pairs, so split on '|' before splitting on ':'.
    parts = strsplit(lines{i}, '|');
    for j = 1:numel(parts)
        % Keys may contain spaces ("Captured at:"), so normalise to a valid
        % fieldname rather than skipping the line.
        tok = regexp(strtrim(parts{j}), '^([A-Za-z][A-Za-z_ ]*):\s*(.*)$', 'tokens', 'once');
        if isempty(tok), continue; end
        key = strrep(strtrim(tok{1}), ' ', '_');
        if ~isfield(prov, key)
            prov.(key) = strtrim(tok{2});
        end
    end
end
end

function [entries, caveats] = collect_analyses(run_dir, run_mode)
% COLLECT_ANALYSES One struct per sub-analysis directory that holds a psa.
%
% A directory without psa_object.mat is not an analysis. That is what drops
% replot_*, which holds only figures regenerated from sweeps already listed
% here: replot_sensitivity repoints a loaded psa's output_dir and saves figures,
% never calling save_object, so the absence is expected rather than a fault. The
% prefix is recognised only to say so -- the psa_object.mat test still does the
% dropping, so a replot folder that somehow gained one would be reported.
caveats = {};
entries = struct([]);

d = dir(run_dir);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.', '..'}));
[~, order] = sort({d.name});
d = d(order);

for i = 1:numel(d)
    mat = fullfile(run_dir, d(i).name, 'psa_object.mat');
    if ~isfile(mat)
        if startsWith(d(i).name, 'replot_')
            caveats{end+1} = sprintf(['Skipped `%s/` -- a replot of analyses ' ...
                'already listed above, holding only regenerated figures. Replot ' ...
                'folders carry no `psa_object.mat` by design, so its parameters ' ...
                'are those of the sweeps it was plotted from.'], d(i).name); %#ok<AGROW>
        else
            caveats{end+1} = sprintf(['Skipped `%s/` -- no `psa_object.mat`, so ' ...
                'no parameters could be recovered from it (incomplete or ' ...
                'interrupted run?).'], d(i).name); %#ok<AGROW>
        end
        continue
    end

    [psa, err] = load_psa(mat);
    if isempty(psa)
        caveats{end+1} = sprintf('Skipped `%s/` -- could not load its psa: %s', ...
            d(i).name, err); %#ok<AGROW>
        continue
    end

    e = struct();
    e.dir_name = d(i).name;
    e.psa = psa;
    e.analysis = analysis_from_prefix(psa.folder_prefix);
    if isempty(e.analysis)
        caveats{end+1} = sprintf(['`%s/` has an unrecognised folder_prefix ' ...
            '(`%s`), so its run-mode timings could not be attributed.'], ...
            d(i).name, psa.folder_prefix); %#ok<AGROW>
    end

    % cfg is what run_mode contributed. Reconstructible from saved data alone:
    % analysis_run_config depends on the preset only through sigma_u_noise,
    % which lives in model_defaults.
    e.cfg = [];
    if ~isempty(e.analysis) && ~isempty(run_mode)
        try
            e.cfg = analysis_run_config(e.analysis, run_mode, psa.model_defaults);
        catch ME
            caveats{end+1} = sprintf(['Could not reconstruct the run-mode ' ...
                'settings for `%s/`: %s'], d(i).name, ME.message); %#ok<AGROW>
        end
    end

    e.preset = recorded_preset(psa, e.cfg);
    e.short  = short_label(e);
    e.sfa_and_std = pick_condition(psa, 'sfa_and_std');
    e.cond_fields = condition_fields(psa);

    if isempty(fieldnames(psa.resolved_defaults))
        caveats{end+1} = sprintf(['`%s/` predates `resolved_defaults`, so its ' ...
            'parameters are known only from `model_defaults`.'], d(i).name); %#ok<AGROW>
    end

    if isempty(entries), entries = e; else, entries(end+1) = e; end %#ok<AGROW>
end
end

function [psa, err] = load_psa(mat_path)
% LOAD_PSA Plain load + pick by class.
%
% Deliberately NOT ParamSpaceAnalysis2.from_dir, which pulls in every 4 MB
% per-condition result file when the saved object carries no results.
psa = [];
err = '';
try
    s = load(mat_path);
    f = fieldnames(s);
    for i = 1:numel(f)
        if isa(s.(f{i}), 'ParamSpaceAnalysis2')
            psa = s.(f{i});
            return
        end
    end
    err = 'no ParamSpaceAnalysis2 object inside';
catch ME
    err = ME.message;
end
end

function a = analysis_from_prefix(prefix)
% ANALYSIS_FROM_PREFIX folder_prefix -> the analysis_run_config key.
switch prefix
    case '1D_sensitivity',  a = 'sensitivity';
    case 'tau_sensitivity', a = 'tau_sensitivity';
    case 'param_space',     a = 'param_space';
    otherwise,              a = '';
end
end

function p = recorded_preset(psa, cfg)
% RECORDED_PRESET The preset as it actually ran: model_defaults minus cfg.model.
%
% Whatever remains that the current preset does not carry is a sub-script
% override (e.g. run_sensitivity_analysis adds tau_b_E_rec on SRNNModel2), which
% is worth surfacing rather than hiding.
p = psa.model_defaults;
if isempty(cfg) || ~isfield(cfg, 'model'), return; end
f = fieldnames(cfg.model);
for i = 1:numel(f)
    if isfield(p, f{i})
        p = rmfield(p, f{i});
    end
end
end

function [p, caveats] = merged_recorded_preset(entries)
% MERGED_RECORDED_PRESET Union of the per-analysis recorded presets.
%
% They should agree; where they do not it is a sub-script override applying to
% some analyses only, which becomes a caveat rather than a silent first-wins.
caveats = {};
p = struct();
if isempty(entries), return; end

p = entries(1).preset;
differing = {};
for k = 2:numel(entries)
    q = entries(k).preset;
    names = union(fieldnames(p), fieldnames(q));
    for i = 1:numel(names)
        n = names{i};
        if ~isfield(q, n)
            differing{end+1} = n; %#ok<AGROW>
        elseif ~isfield(p, n)
            p.(n) = q.(n);
            differing{end+1} = n; %#ok<AGROW>
        elseif ~same_value(p.(n), q.(n))
            differing{end+1} = n; %#ok<AGROW>
        end
    end
end
differing = unique(differing);
if ~isempty(differing)
    caveats{end+1} = sprintf(['The analyses did not all receive the same preset ' ...
        'fields -- these differ between them and are usually sub-script ' ...
        'overrides: %s. The per-analysis value is shown in the parameter list.'], ...
        strjoin(cellfun(@(s) ['`' s '`'], differing, 'UniformOutput', false), ', '));
end
end

function [p, ok, caveats] = current_preset(preset_name)
% CURRENT_PRESET What srnn_param_preset returns for this name TODAY.
caveats = {};
p = struct();
ok = false;
if isempty(preset_name)
    % Runs predating the preset mechanism have a manifest but no preset name.
    caveats{end+1} = ['This run recorded no preset name, so the *What the preset ' ...
        'set* table shows the overrides that reached the sweep without saying ' ...
        'which preset they came from.'];
    return
end
try
    p = srnn_param_preset(preset_name);
    ok = true;
catch ME
    caveats{end+1} = sprintf(['Preset `%s` could not be looked up in the ' ...
        'current source (%s), so no drift check was possible.'], ...
        preset_name, ME.message);
end
end

function cond = pick_condition(psa, name)
% PICK_CONDITION The named condition struct, or [] when the run omitted it.
cond = [];
for i = 1:numel(psa.conditions)
    c = psa.conditions{i};
    if isfield(c, 'name') && strcmp(c.name, name)
        cond = c;
        return
    end
end
end

function names = condition_fields(psa)
% CONDITION_FIELDS Every field any condition actually sets, minus 'name'.
names = {};
for i = 1:numel(psa.conditions)
    names = union(names, setdiff(fieldnames(psa.conditions{i}), {'name'}));
end
end

function s = short_label(e)
% SHORT_LABEL Compact name for an analysis, used in per-analysis sub-lists.
switch e.analysis
    case 'sensitivity'
        swept = setdiff(e.psa.grid_params, {'reps'});
        if isempty(swept)
            s = 'sensitivity';
        else
            s = sprintf('sensitivity(%s)', strjoin(swept, ', '));
        end
    case 'tau_sensitivity', s = 'tau_sensitivity';
    case 'param_space',     s = 'param_space';
    otherwise,              s = e.dir_name;
end
end

%% ------------------------------------------------------------------------
%% Markdown sections
%% ------------------------------------------------------------------------

function write_header(fid, run_name, preset_name, run_mode, model_class)
fprintf(fid, '%s\n', ['# Run parameters — `' run_name '`']);
fprintf(fid, '\n');
fprintf(fid, '%s\n', sprintf('**Preset:** `%s`  ·  **Run mode:** `%s`  ·  **Model class:** `%s`', ...
    or_unknown(preset_name), or_unknown(run_mode), or_unknown(model_class)));
fprintf(fid, '\n');
fprintf(fid, '%s\n', ['Generated by `write_run_parameters_md` from this run''s own saved artifacts ' ...
    '(`run_manifest.mat` and each sub-analysis''s `psa_object.mat`). Every value below is the ' ...
    'one **actually used**, not what the current source would produce today.']);
fprintf(fid, '\n');
fprintf(fid, '%s\n', ['Parameters are listed **as run in the `sfa_and_std` condition** (SFA and ' ...
    'short-term depression both on). The other three adaptation regimes differ from it only in ' ...
    'the condition-owned fields listed under *Adaptation conditions*.']);
fprintf(fid, '\n');
end

function write_provenance(fid, man, prov)
fprintf(fid, '%s\n', '## Run provenance');
fprintf(fid, '\n');
fprintf(fid, '%s\n', '| | |');
fprintf(fid, '%s\n', '|---|---|');
row = @(k, v) fprintf(fid, '%s\n', ['| ' k ' | ' escape_cell(v) ' |']);

row('Started',     getfield_or(prov, 'Captured_at', getfield_or(man, 'timestamp', '(unknown)')));
row('Folder stamp', getfield_or(man, 'timestamp', '(unknown)'));
commit = getfield_or(man, 'git_commit_short', getfield_or(man, 'git_commit', '(unknown)'));
dirty  = getfield_or(man, 'git_dirty', []);
if islogical(dirty) || isnumeric(dirty)
    if ~isempty(dirty) && dirty
        commit = [commit ' (working tree **dirty**)'];
    elseif ~isempty(dirty)
        commit = [commit ' (clean)'];
    end
end
row('Git commit',  commit);
row('Git branch',  getfield_or(man, 'git_branch', '(unknown)'));
row('Hostname',    getfield_or(prov, 'hostname', '(unknown)'));
row('Platform',    getfield_or(prov, 'platform', '(unknown)'));
row('MATLAB',      getfield_or(man, 'matlab_version', getfield_or(prov, 'matlab', '(unknown)')));
row('Figure policy', getfield_or(man, 'master_save_figs', '(unknown)'));
fprintf(fid, '\n');
end

function write_preset_section(fid, rec_preset, cur_preset, cur_ok, preset_name, entries)
fprintf(fid, '%s\n', '## What the preset set');
fprintf(fid, '\n');

names = sort_names(fieldnames(rec_preset));
if isempty(names)
    fprintf(fid, '%s\n', 'The preset set nothing (no `model_defaults` were recorded).');
    fprintf(fid, '\n');
    return
end

fprintf(fid, '%s\n', ['These are the model overrides that reached the sweep, recovered as ' ...
    '`model_defaults` minus the run-mode settings. Anything a preset does not carry — ' ...
    'timings, sweep sizes, adaptation counts — appears in the sections below instead.']);
fprintf(fid, '\n');

% Drift is checked against today's srnn_param_preset. A preset edited after the
% run will disagree, and the recorded column is the authoritative one.
drift = {};
if cur_ok
    for i = 1:numel(names)
        n = names{i};
        if ~isfield(cur_preset, n) || ~same_value(rec_preset.(n), cur_preset.(n))
            drift{end+1} = n; %#ok<AGROW>
        end
    end
end

if ~isempty(drift)
    fprintf(fid, '%s\n', sprintf(['> ⚠ **Preset drift.** %d field(s) differ from ' ...
        '`srnn_param_preset(''%s'')` as it stands today. The *as run* column is ' ...
        'authoritative; the preset definition has changed since this run.'], ...
        numel(drift), preset_name));
    fprintf(fid, '\n');
end

fprintf(fid, '%s\n', '| Parameter | Value as run | Current `srnn_param_preset` |');
fprintf(fid, '%s\n', '|---|---|---|');
for i = 1:numel(names)
    n = names{i};
    if ~cur_ok
        cur_txt = '*(not checked)*';
    elseif ~isfield(cur_preset, n)
        cur_txt = sprintf('*(absent — %s)*', override_origin(n, entries));
    elseif same_value(rec_preset.(n), cur_preset.(n))
        cur_txt = '✓ same';
    else
        % Two structs with the same fields render identically one-line, so a
        % flagged struct has to say which leaves actually moved.
        cur_txt = ['⚠ ' fmt_diff(rec_preset.(n), cur_preset.(n))];
    end
    fprintf(fid, '%s\n', ['| `' n '` | `' escape_cell(fmt_value(rec_preset.(n))) '` | ' cur_txt ' |']);
end
fprintf(fid, '\n');
end

function s = override_origin(name, entries)
% OVERRIDE_ORIGIN Whether a field the current preset lacks came from every
% analysis or only some -- "only some" means a sub-script added it.
have = arrayfun(@(e) isfield(e.preset, name), entries);
if all(have)
    s = 'not in the preset today';
else
    s = sprintf('script override, %d of %d analyses', sum(have), numel(have));
end
end

function write_analyses_section(fid, entries)
fprintf(fid, '%s\n', '## Analyses in this run');
fprintf(fid, '\n');
if isempty(entries)
    fprintf(fid, '%s\n', 'None recovered.');
    fprintf(fid, '\n');
    return
end

fprintf(fid, '%s\n', '| Directory | Type | Swept | Levels | Reps | Solver | fs | `T_range` | LLE window | Conds | Grid pts |');
fprintf(fid, '%s\n', '|---|---|---|---|---|---|---|---|---|---|---|');
for k = 1:numel(entries)
    e = entries(k);
    psa = e.psa;

    swept = setdiff(psa.grid_params, {'reps'}, 'stable');
    swept_txt = strjoin(cellfun(@(s) ['`' s '`'], swept, 'UniformOutput', false), ', ');
    if isempty(swept_txt), swept_txt = '—'; end

    reps_txt = '—';
    if isfield(psa.explicit_vectors, 'reps')
        reps_txt = sprintf('%d', numel(psa.explicit_vectors.reps));
    end

    [solver, fs_txt, tr_txt, lya_txt] = timing_cells(e);

    npts = '—';
    if ~isempty(psa.num_combinations)
        npts = sprintf('%d', psa.num_combinations);
    end

    fprintf(fid, '%s\n', ['| `' e.dir_name '` | ' or_unknown(e.analysis) ' | ' swept_txt ...
        ' | ' sprintf('%d', psa.n_levels) ' | ' reps_txt ' | ' solver ' | ' fs_txt ...
        ' | ' tr_txt ' | ' lya_txt ' | ' sprintf('%d', numel(psa.conditions)) ...
        ' | ' npts ' |']);
end
fprintf(fid, '\n');

% The axis descriptors do not fit the table, so they follow it.
fprintf(fid, '%s\n', 'Swept axes in detail:');
fprintf(fid, '\n');
for k = 1:numel(entries)
    e = entries(k);
    fprintf(fid, '%s\n', ['- `' e.dir_name '`']);
    if isempty(e.psa.grid_params)
        fprintf(fid, '%s\n', '  - *(no grid axes)*');
    end
    for i = 1:numel(e.psa.grid_params)
        n = e.psa.grid_params{i};
        fprintf(fid, '%s\n', ['  - `' n '` — ' fmt_axis(e, n)]);
    end
end
fprintf(fid, '\n');
end

function [solver, fs_txt, tr_txt, lya_txt] = timing_cells(e)
% TIMING_CELLS The run-mode-owned settings, read back from what the model got.
%
% Read from resolved_defaults rather than from cfg, so the table reports what
% the model was actually constructed with; cfg only says what run_mode asked for.
rd = e.psa.resolved_defaults;
solver = value_or(rd, 'ode_solver', e.psa.model_defaults, '(unknown)');
fs_txt = value_or(rd, 'fs',         e.psa.model_defaults, '(unknown)');
tr_txt = value_or(rd, 'T_range',    e.psa.model_defaults, '(unknown)');
% An absent lya_T_interval in cfg.model is analysis_run_config saying "leave it
% to the class default", which is worth naming rather than showing a bare value.
lya_txt = value_or(rd, 'lya_T_interval', e.psa.model_defaults, '');
if ~isempty(e.cfg) && isfield(e.cfg, 'model') && ~isfield(e.cfg.model, 'lya_T_interval')
    if isempty(lya_txt)
        lya_txt = '*(class default)*';
    else
        lya_txt = [lya_txt ' *(class default)*'];
    end
elseif isempty(lya_txt)
    lya_txt = '(unknown)';
end
end

function s = value_or(rd, name, fallback, default_txt)
if isfield(rd, name)
    s = ['`' escape_cell(fmt_value(rd.(name))) '`'];
elseif isfield(fallback, name)
    s = ['`' escape_cell(fmt_value(fallback.(name))) '`'];
else
    s = default_txt;
end
end

function write_conditions_section(fid, entries, model_class)
fprintf(fid, '%s\n', '## Adaptation conditions');
fprintf(fid, '\n');
if isempty(entries)
    fprintf(fid, '%s\n', 'None recovered.');
    fprintf(fid, '\n');
    return
end

% The conditions are the same object across analyses except where one runs a
% subset (tau_sensitivity runs sfa_and_std alone), so show the richest set.
[~, widest] = max(arrayfun(@(e) numel(e.psa.conditions), entries));
conds = entries(widest).psa.conditions;
fields = setdiff(condition_fields(entries(widest).psa), {'name'});
fields = sort_names(fields);

fprintf(fid, '%s\n', ['A condition owns whichever model fields it sets, so those fields are ' ...
    'excluded from the frozen defaults and appear here instead. Each is applied on top of ' ...
    'everything in the parameter list below.']);
fprintf(fid, '\n');

hdr = ['| Condition | ' strjoin(cellfun(@(f) ['`' f '`'], fields, 'UniformOutput', false), ' | ') ' |'];
sep = ['|---|' repmat('---|', 1, numel(fields))];
fprintf(fid, '%s\n', hdr);
fprintf(fid, '%s\n', sep);
for i = 1:numel(conds)
    c = conds{i};
    cells = cell(1, numel(fields));
    for j = 1:numel(fields)
        if isfield(c, fields{j})
            cells{j} = ['`' escape_cell(fmt_value(c.(fields{j}))) '`'];
        else
            cells{j} = '*(not set)*';
        end
    end
    fprintf(fid, '%s\n', ['| `' c.name '` | ' strjoin(cells, ' | ') ' |']);
end
fprintf(fid, '\n');

if strcmp(model_class, 'SRNNModel2')
    fprintf(fid, '%s\n', ['`n_a_E` is the number of SFA timescales on the excitatory population; ' ...
        '`n_b_E` = 1 puts short-term depression on **every** outgoing excitatory synapse.']);
    fprintf(fid, '\n');
end

% synapse_config nests too deeply for a table cell, so expand it once.
for i = 1:numel(conds)
    c = conds{i};
    if isfield(c, 'synapse_config') && ~isempty(fieldnames(c.synapse_config))
        fprintf(fid, '%s\n', ['`synapse_config` for `' c.name '` (identical for every ' ...
            'condition that has depressing routes):']);
        fprintf(fid, '\n');
        fmt_value_block(fid, 'synapse_config', c.synapse_config, 0, '');
        fprintf(fid, '\n');
        break
    end
end

% A subset run is worth naming: tau_sensitivity deliberately runs one condition.
for k = 1:numel(entries)
    if numel(entries(k).psa.conditions) < numel(conds)
        names = cellfun(@(c) ['`' c.name '`'], entries(k).psa.conditions, 'UniformOutput', false);
        fprintf(fid, '%s\n', sprintf('`%s` ran only %s.', entries(k).dir_name, strjoin(names, ', ')));
        fprintf(fid, '\n');
    end
end
end

function caveats = write_parameters_section(fid, entries, rec_preset, model_class)
% WRITE_PARAMETERS_SECTION The unified list: every parameter, with its source.
caveats = {};
fprintf(fid, '%s\n', '## All parameters used');
fprintf(fid, '\n');
if isempty(entries)
    fprintf(fid, '%s\n', 'None recovered.');
    fprintf(fid, '\n');
    return
end

[ref, borrowed] = default_reference_model(model_class, rec_preset);
if isempty(ref)
    caveats{end+1} = sprintf(['`%s` could not be constructed with no arguments, ' ...
        'so *class default* and *derived at construction* could not be told apart; ' ...
        'both are reported as *not set explicitly*.'], model_class);
elseif ~isempty(borrowed)
    caveats{end+1} = sprintf(['`%s` has required constructor arguments, so the ' ...
        'reference used to recognise class defaults had to borrow %s from the ' ...
        'preset. Scalar **aliases** of those (`f_E` for `f`, `mu_EE_relative` and ' ...
        '`sigma_EE_relative` for the tilde blocks, and so on) therefore read as ' ...
        '*class default* below when the preset is what actually set them — see the ' ...
        'preset table above for the values it carried.'], model_class, ...
        strjoin(cellfun(@(s) ['`' s '`'], borrowed, 'UniformOutput', false), ', '));
end

fprintf(fid, '%s\n', ['Everything the model ran with, whatever set it. The source tag is the ' ...
    '**last** thing to write the value, in the order the sweep applies them: grid axis > ' ...
    'condition > run mode > preset > class default. Where the analyses disagree — the swept ' ...
    'axes and the timings — the per-analysis values are broken out.']);
fprintf(fid, '\n');

names = sort_names(universe(entries));
for i = 1:numel(names)
    write_parameter(fid, names{i}, entries, ref);
end
fprintf(fid, '\n');
end

function names = universe(entries)
% UNIVERSE Every parameter name any analysis knows about.
names = {};
for k = 1:numel(entries)
    psa = entries(k).psa;
    names = union(names, fieldnames(psa.resolved_defaults));
    names = union(names, psa.grid_params);
    names = union(names, entries(k).cond_fields);
    names = union(names, fieldnames(psa.model_defaults));
end
end

function write_parameter(fid, name, entries, ref)
% WRITE_PARAMETER One list item, collapsed across analyses where they agree.
recs = arrayfun(@(e) param_record(e, name, ref), entries);
recs = recs(~strcmp({recs.kind}, 'absent'));
if isempty(recs)
    return
end

% Group by what would be printed: same rendering and same source => one line.
keys = arrayfun(@(r) [r.kind '|' r.text], recs, 'UniformOutput', false);
[uk, ia] = unique(keys, 'stable');

if isscalar(uk)
    r = recs(ia(1));
    if strcmp(r.kind, 'grid')
        fprintf(fid, '%s\n', ['- **' name '** — ' r.text ' — *' r.source '*']);
    elseif isstruct(r.value) && isscalar(r.value) && ~isempty(fieldnames(r.value))
        fmt_value_block(fid, ['**' name '**'], r.value, 0, r.source);
    else
        fprintf(fid, '%s\n', ['- **' name '** = `' r.text '` — *' r.source '*']);
    end
    return
end

fprintf(fid, '%s\n', ['- **' name '** — differs by analysis:']);
for g = 1:numel(uk)
    sel = strcmp(keys, uk{g});
    r = recs(ia(g));
    who = strjoin(unique({recs(sel).where}, 'stable'), ', ');
    if strcmp(r.kind, 'grid')
        fprintf(fid, '%s\n', ['  - ' r.text ' — in ' who ' — *' r.source '*']);
    else
        fprintf(fid, '%s\n', ['  - `' r.text '` — in ' who ' — *' r.source '*']);
    end
end
end

function r = param_record(e, name, ref)
% PARAM_RECORD What NAME was, and what set it, for one analysis.
%
% Precedence mirrors run_single_job's last-write-wins constructor call: grid
% beats condition beats run mode beats preset. run mode beats preset because the
% sub-scripts merge as merge_struct(preset_defaults, cfg.model).
psa = e.psa;
r = struct('kind', 'absent', 'source', '', 'text', '', 'value', [], 'where', e.short);

if ismember(name, psa.grid_params) && strcmp(name, 'reps')
    % reps is added as a grid axis rather than as a property of the sweep, so it
    % is a repetition COUNT here -- and it collides with the model's own scalar
    % `reps` property, which is what the analyses with no reps axis report.
    r.kind = 'grid';
    r.source = 'sweep axis (repetitions per grid point)';
    r.text = fmt_axis(e, name);
    return
end

if ismember(name, psa.grid_params)
    r.kind = 'grid';
    r.source = 'grid axis';
    r.text = fmt_axis(e, name);
    return
end

if ~isempty(e.sfa_and_std) && isfield(e.sfa_and_std, name)
    r.kind = 'condition';
    r.source = 'condition `sfa_and_std`';
    r.value = e.sfa_and_std.(name);
    r.text = fmt_value(r.value);
    return
end

% Value always comes from resolved_defaults (what the model was built with);
% model_defaults is the fallback for runs that predate the freeze.
if isfield(psa.resolved_defaults, name)
    r.value = psa.resolved_defaults.(name);
elseif isfield(psa.model_defaults, name)
    r.value = psa.model_defaults.(name);
else
    return  % absent
end
r.text = fmt_value(r.value);

if ~isempty(e.cfg) && isfield(e.cfg, 'model') && isfield(e.cfg.model, name)
    r.kind = 'run_mode';
    r.source = sprintf('run mode `%s`', e.cfg.run_mode);
elseif isfield(e.preset, name)
    r.kind = 'preset';
    r.source = 'preset';
elseif isempty(ref)
    r.kind = 'default';
    r.source = 'not set explicitly';
elseif isprop(ref, name) && same_value(ref.(name), r.value)
    r.kind = 'default';
    r.source = 'class default';
else
    r.kind = 'derived';
    r.source = 'derived at construction';
end
end

function [m, borrowed] = default_reference_model(model_class, rec_preset)
% DEFAULT_REFERENCE_MODEL A model carrying nothing but its class defaults.
%
% SRNNCellTypePairs has required constructor arguments (it is general over cell
% types, so there is no sensible default for them), so a bare feval fails. Retry
% with only those arguments, taken from the preset: everything else then still
% reads as the class default it is. BORROWED names the arguments that had to
% come from the preset -- their scalar aliases (f_E from f, mu_EE_relative from
% mu_tilde_relative) then compare equal to the reference and would otherwise be
% reported as class defaults when the preset is what really set them.
m = [];
borrowed = {};
try
    m = feval(model_class);
    return
catch
end
required = {'n_cellTypes', 'cell_type_names', 'f', 'mu_tilde_relative', 'sigma_tilde_relative'};
args = {};
for i = 1:numel(required)
    if isfield(rec_preset, required{i})
        args = [args, {required{i}, rec_preset.(required{i})}]; %#ok<AGROW>
        borrowed{end+1} = required{i}; %#ok<AGROW>
    end
end
if isempty(args)
    borrowed = {};
    return
end
try
    m = feval(model_class, args{:});
catch
    m = [];
    borrowed = {};
end
end

function write_caveats(fid, caveats)
if isempty(caveats), return; end
fprintf(fid, '%s\n', '## Caveats');
fprintf(fid, '\n');
for i = 1:numel(caveats)
    fprintf(fid, '%s\n', ['- ' caveats{i}]);
end
fprintf(fid, '\n');
end

%% ------------------------------------------------------------------------
%% Formatting
%% ------------------------------------------------------------------------

function s = fmt_axis(e, name)
% FMT_AXIS How one grid axis was swept, in words.
psa = e.psa;

if isfield(psa.vector_param_config, name)
    v = psa.vector_param_config.(name);
    s = sprintf(['swept %s element %s → %s over %d levels (%s level spacing); ' ...
        '%d elements, fixed end %s, %s spacing'], ...
        v.vary_element, num2str(v.vary_range(1), '%.6g'), num2str(v.vary_range(2), '%.6g'), ...
        psa.n_levels, v.level_spacing, v.n_elements, ...
        num2str(v.fixed_value, '%.6g'), v.spacing);
    return
end

if isfield(psa.explicit_vectors, name)
    v = psa.explicit_vectors.(name);
    s = sprintf('swept over %d explicit values %s → %s', numel(v), ...
        num2str(min(v), '%.6g'), num2str(max(v), '%.6g'));
    return
end

if isfield(psa.param_ranges, name)
    rg = psa.param_ranges.(name);
    s = sprintf('swept %s → %s, %d levels', num2str(rg(1), '%.6g'), ...
        num2str(rg(2), '%.6g'), psa.n_levels);
    if ismember(name, psa.integer_params)
        s = [s ' (rounded to integers)'];
    end
    return
end

s = sprintf('swept over %d levels', psa.n_levels);
end

function s = fmt_value(v)
% FMT_VALUE One-line rendering of any parameter value.
%
% Anything that cannot be said on one line is summarised here; the caller
% switches to fmt_value_block when the nesting is worth showing.
if ischar(v)
    s = ['''' v ''''];
elseif isstring(v)
    if isscalar(v), s = ['"' char(v) '"']; else, s = sprintf('<%s string>', dims(v)); end
elseif islogical(v)
    if isscalar(v)
        if v, s = 'true'; else, s = 'false'; end
    else
        s = mat2str(v);
    end
elseif isnumeric(v)
    if isempty(v)
        s = '[]';
    elseif isscalar(v)
        s = num2str(v, '%.6g');
    elseif numel(v) <= 25
        s = mat2str(v, 6);
    else
        s = sprintf('<%s %s, %s to %s>', dims(v), class(v), ...
            num2str(min(v(:)), '%.6g'), num2str(max(v(:)), '%.6g'));
    end
elseif iscell(v)
    if isempty(v)
        s = '{}';
    elseif iscellstr(v)
        s = ['{' strjoin(cellfun(@(c) ['''' c ''''], v(:)', 'UniformOutput', false), ', ') '}'];
    elseif numel(v) <= 6
        s = ['{' strjoin(cellfun(@fmt_value, v(:)', 'UniformOutput', false), ', ') '}'];
    else
        s = sprintf('<%s cell>', dims(v));
    end
elseif isstruct(v)
    if isscalar(v)
        f = fieldnames(v);
        if isempty(f)
            s = 'struct with no fields';
        else
            s = ['struct: ' strjoin(f(:)', ', ')];
        end
    else
        s = sprintf('<%s struct>', dims(v));
    end
elseif isa(v, 'function_handle')
    s = func2str(v);
else
    s = sprintf('<%s %s>', dims(v), class(v));
end
end

function s = fmt_diff(a, b)
% FMT_DIFF How B differs from A, for the drift column.
%
% For a scalar struct the one-line rendering is just its field list, so two
% structs that differ only in a leaf would read as identical while flagged as
% drifting. Walk them and name the leaves instead.
if isstruct(a) && isstruct(b) && isscalar(a) && isscalar(b)
    parts = struct_diff_leaves(a, b, '');
    if ~isempty(parts)
        s = strjoin(parts, '; ');
        return
    end
end
s = ['`' fmt_value(b) '`'];
end

function parts = struct_diff_leaves(a, b, prefix)
% STRUCT_DIFF_LEAVES `path: was → now` for every leaf that moved.
parts = {};
names = union(fieldnames(a), fieldnames(b));
for i = 1:numel(names)
    n = names{i};
    path = [prefix n];
    if ~isfield(b, n)
        parts{end+1} = ['`' path '` removed']; %#ok<AGROW>
    elseif ~isfield(a, n)
        parts{end+1} = ['`' path '` added = `' fmt_value(b.(n)) '`']; %#ok<AGROW>
    elseif same_value(a.(n), b.(n))
        continue
    elseif isstruct(a.(n)) && isstruct(b.(n)) && isscalar(a.(n)) && isscalar(b.(n))
        parts = [parts, struct_diff_leaves(a.(n), b.(n), [path '.'])]; %#ok<AGROW>
    else
        parts{end+1} = ['`' path '`: `' fmt_value(a.(n)) '` → `' fmt_value(b.(n)) '`']; %#ok<AGROW>
    end
end
end

function fmt_value_block(fid, label, v, indent, source)
% FMT_VALUE_BLOCK A nested struct as an indented markdown sub-list.
%
% Tables cannot hold synapse_config's route nesting and a disp() dump is
% unreadable, so nested structs get real list nesting instead.
pad = repmat(' ', 1, indent);
if isempty(source)
    fprintf(fid, '%s\n', [pad '- ' label]);
else
    fprintf(fid, '%s\n', [pad '- ' label ' — *' source '*']);
end
if indent >= 8
    fprintf(fid, '%s\n', [pad '  - …']);
    return
end
f = fieldnames(v);
for i = 1:numel(f)
    x = v.(f{i});
    if isstruct(x) && isscalar(x) && ~isempty(fieldnames(x))
        fmt_value_block(fid, ['`' f{i} '`'], x, indent + 2, '');
    else
        fprintf(fid, '%s\n', [pad '  - `' f{i} '` = `' fmt_value(x) '`']);
    end
end
end

function tf = same_value(a, b)
% SAME_VALUE Equality that copes with function handles nested in structs.
%
% isequaln compares two distinct-but-identical anonymous handles as unequal, so
% handles are compared by their text and structs are walked field by field.
if isa(a, 'function_handle') || isa(b, 'function_handle')
    tf = isa(a, 'function_handle') && isa(b, 'function_handle') && ...
        strcmp(func2str(a), func2str(b));
    return
end
if isstruct(a) && isstruct(b) && isscalar(a) && isscalar(b)
    fa = sort(fieldnames(a));
    fb = sort(fieldnames(b));
    tf = isequal(fa, fb);
    if ~tf, return; end
    for i = 1:numel(fa)
        if ~same_value(a.(fa{i}), b.(fa{i}))
            tf = false;
            return
        end
    end
    return
end
if iscell(a) && iscell(b)
    tf = isequal(size(a), size(b));
    if ~tf, return; end
    for i = 1:numel(a)
        if ~same_value(a{i}, b{i})
            tf = false;
            return
        end
    end
    return
end
tf = isequaln(a, b);
end

%% ------------------------------------------------------------------------
%% Small utilities
%% ------------------------------------------------------------------------

function v = getfield_or(s, name, default)
if isstruct(s) && isscalar(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
    if isnumeric(v) || islogical(v)
        return
    end
    if isstring(v), v = char(v); end
else
    v = default;
end
end

function s = or_unknown(s)
if isempty(s), s = '(unknown)'; end
end

function s = dims(v)
s = strjoin(arrayfun(@(k) sprintf('%d', k), size(v), 'UniformOutput', false), 'x');
end

function names = sort_names(names)
[~, order] = sort(lower(names));
names = names(order);
names = names(:)';
end

function s = escape_cell(s)
% ESCAPE_CELL A markdown table cell cannot contain a bare pipe.
s = strrep(s, '|', '\|');
end

function p = strip_trailing_sep(p)
while ~isempty(p) && (p(end) == filesep || p(end) == '/')
    p(end) = [];
end
end
