function combined_dir = combine_runs(run_dirs, out_dir)
% COMBINE_RUNS Pool multiple run_all_<dt> runs of the SAME config into single plots.
%
%   combined_dir = COMBINE_RUNS(run_dirs)
%   combined_dir = COMBINE_RUNS(run_dirs, out_dir)
%
% Inputs:
%   run_dirs : cell array of paths to run_all_<dt> folders to combine. Each must
%              have been produced with a DISTINCT network seed offset (see below)
%              so the pooled data are independent draws, not duplicates.
%   out_dir  : optional output folder for the combined figures
%              (default: <parent-of-run_dirs{1}>/combined_runs_<dt>).
%
% Output:
%   combined_dir : the folder the combined figures were written to.
%
% What it pools (all by concatenating metric samples):
%   - 1D sensitivity : reps per swept level          (plot_sensitivity + pool_with)
%   - tau sensitivity: reps per swept level          (plot_sensitivity + pool_with)
%   - param space    : the whole-grid LLE distribution(plot + pool_with)
%   - DC LLE         : per-seed level means           (stack + recompute mean/std)
%
% Output layout under combined_dir:
%   combined_manifest.{mat,txt}  which runs were pooled + their settings
%   figures/                     tau, param-space, DC figures
%   sensitivity/figures/         per-param 1D sensitivity figures PLUS the stacked
%                                sensitivity_LLE_combined.{fig,png} (all 1D
%                                sweeps in one figure, like replot_sensitivity).
%
% Safety: runs are pooled only if their configs match (ParamSpaceAnalysis2.
% same_config / the DC config struct + dc_levels) AND they used DISTINCT seed
% offsets. Two runs with the same offset would be identical data; that is an
% error. NOTE: run_all_analyses auto-increments the offset only WITHIN a MATLAB
% session (run_index 0,1,2,...). A fresh session restarts at offset 0, so to
% pool runs from different sessions set `run_index = <N>` in the base workspace
% before running run_all_analyses so it uses offset N*1e6.
%
% No simulations are re-run.

    if nargin < 1 || isempty(run_dirs)
        error('combine_runs:NoRuns', 'Provide a cell array of run_all_<dt> folders.');
    end
    if ischar(run_dirs) || isstring(run_dirs)
        run_dirs = cellstr(run_dirs);
    end
    run_dirs = run_dirs(:)';
    for i = 1:numel(run_dirs)
        if ~isfolder(run_dirs{i})
            error('combine_runs:BadDir', 'Not a folder: %s', run_dirs{i});
        end
    end
    if numel(run_dirs) < 2
        warning('combine_runs:OneRun', ...
            'Only %d run given; nothing to pool (will still plot it).', numel(run_dirs));
    end

    setup_paths();

    if nargin < 2 || isempty(out_dir)
        dt_str = lower(datestr(now, 'mmm_dd_yy_HH_MM')); %#ok<TNOW1,DATST>
        out_dir = fullfile(fileparts(run_dirs{1}), sprintf('combined_runs_%s', dt_str));
    end
    if ~isfolder(out_dir); mkdir(out_dir); end
    combined_dir = out_dir;
    fig_dir = fullfile(out_dir, 'figures');
    fprintf('Combining %d runs -> %s\n\n', numel(run_dirs), out_dir);

    % Load run manifests and verify compatibility. The manifest fingerprints the
    % nonlinearity (activation + S_a/S_c), which same_config cannot see, so this
    % is what guards against pooling e.g. a piecewise run with a logistic one.
    run_info = check_manifests(run_dirs);

    % Record which runs were combined (provenance) into the output dir.
    write_combined_manifest(out_dir, run_dirs, run_info);

    % --- 1D sensitivity: pool per swept param, then stack into one figure (like
    % replot_sensitivity + assemble_sensitivity_figure). Uses its OWN subfolder
    % so the assembly doesn't pick up the tau figures, which share the same
    % "LLE Sensitivity - ..." figure name. ---
    sens_dir = fullfile(out_dir, 'sensitivity');
    combine_psa_sweep(run_dirs, '1D_sensitivity_*', fullfile(sens_dir, 'figures'), 'sensitivity', [-2, 2]);
    try
        assemble_sensitivity_figure(sens_dir, 'LLE');
        assemble_sensitivity_figure(sens_dir, 'mean_rate');
    catch ME
        warning('combine_runs:AssembleFailed', ...
            'Sensitivity assembly skipped: %s', ME.message);
    end

    % --- tau sensitivity ---
    combine_psa_sweep(run_dirs, 'tau_sensitivity_*', fig_dir, 'tau', [-1.5, 1.5]);

    % --- Param space: whole-grid distribution (plot) ---
    combine_param_space(run_dirs, fig_dir);

    % --- DC LLE: stack per-seed level means ---
    combine_dc(run_dirs, fig_dir);

    fprintf('Done. Combined figures in:\n  %s\n', out_dir);
    fprintf('  (stacked 1D sensitivity: %s)\n', fullfile(sens_dir, 'figures'));
end

% ======================================================================
%  LOCAL FUNCTIONS
% ======================================================================
function info = check_manifests(run_dirs)
    % Load each run_manifest.mat, print a summary, require the present manifests
    % to agree on the nonlinearity fingerprint, and RETURN a per-run provenance
    % struct array (for the combined manifest). Missing manifests (e.g. runs
    % predating this feature) are warned about, not fatal.
    manifests = {}; present = [];  missing = [];
    info = struct('dir', {}, 'name', {}, 'has_manifest', {}, 'run_mode', {}, ...
        'seed_offset', {}, 'activation', {}, 'git_commit_short', {});
    fprintf('Run manifests:\n');
    for i = 1:numel(run_dirs)
        [~, nm] = fileparts(run_dirs{i});
        entry = struct('dir', run_dirs{i}, 'name', nm, 'has_manifest', false, ...
            'run_mode', '', 'seed_offset', NaN, 'activation', '', 'git_commit_short', '');
        mf = fullfile(run_dirs{i}, 'run_manifest.mat');
        if exist(mf, 'file')
            S = load(mf); m = S.run_manifest;
            manifests{end+1} = m; present(end+1) = i; %#ok<AGROW>
            entry.has_manifest = true;
            entry.run_mode = getf(m, 'run_mode', '?');
            entry.seed_offset = getf(m, 'master_seed_offset', NaN);
            entry.activation = getf(m, 'activation', '');
            entry.git_commit_short = getf(m, 'git_commit_short', '');
            fprintf('  run %d: run_mode=%s seed_offset=%g activation=%s\n', ...
                i, entry.run_mode, entry.seed_offset, getf(m, 'activation', '<none>'));
        else
            missing(end+1) = i; %#ok<AGROW>
            fprintf('  run %d: (no run_manifest.mat)\n', i);
        end
        info(i) = entry; %#ok<AGROW>
    end
    fprintf('\n');

    if ~isempty(missing)
        warning('combine_runs:NoManifest', ...
            ['Run(s) %s have no run_manifest.mat; their nonlinearity cannot be ', ...
             'verified. Proceeding -- ensure they used the same activation.'], ...
            mat2str(missing));
    end

    flds = {'activation', 'activation_derivative', 'S_a', 'S_c'};
    for j = 2:numel(manifests)
        for f = 1:numel(flds)
            k = flds{f};
            if ~isequaln(getf(manifests{1}, k, []), getf(manifests{j}, k, []))
                error('combine_runs:NonlinearityMismatch', ...
                    ['Runs %d and %d disagree on ''%s'' -- refusing to pool ', ...
                     'different nonlinearities.\n  run %d: %s\n  run %d: %s'], ...
                    present(1), present(j), k, ...
                    present(1), mfval(getf(manifests{1}, k, [])), ...
                    present(j), mfval(getf(manifests{j}, k, [])));
            end
        end
    end
end

function combine_psa_sweep(run_dirs, pattern, fig_dir, tag, lle_hist_range)
    % Group matching PSA subfolders across runs by swept parameter, validate,
    % and pool via plot_sensitivity('pool_with', ...).
    groups = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:numel(run_dirs)
        listing = dir(fullfile(run_dirs{i}, pattern));
        listing = listing([listing.isdir]);
        for k = 1:numel(listing)
            src = fullfile(listing(k).folder, listing(k).name);
            pf = fullfile(src, 'psa_object.mat');
            if ~exist(pf, 'file'); continue; end
            S = load(pf); fns = fieldnames(S); psa = S.(fns{1});
            swept = setdiff(psa.grid_params, {'reps'}, 'stable');
            if isempty(swept); continue; end
            key = strjoin(swept, '+');
            if ~isKey(groups, key); groups(key) = {}; end
            g = groups(key); g{end+1} = psa; groups(key) = g; %#ok<AGROW>
        end
    end

    ks = keys(groups);
    for ki = 1:numel(ks)
        psas = groups(ks{ki});
        [ok, reason] = validate_psa_group(psas);
        if ~ok
            warning('combine_runs:incompatible', ...
                'Skipping %s "%s": %s', tag, ks{ki}, reason);
            continue;
        end
        base = psas{1}; others = psas(2:end);
        base.output_dir = fig_dir;   % not used for saving here, but harmless
        base.plot_sensitivity('metric', 'LLE', 'hist_range', lle_hist_range, 'pool_with', others);
        base.plot_sensitivity('metric', 'mean_rate', 'pool_with', others);
        name = sprintf('combined_%s_%s', tag, matlab.lang.makeValidName(ks{ki}));
        save_some_figs_to_folder_2(fig_dir, name, [], {'fig', 'png'});
        close all;
        fprintf('  [%s] pooled %d runs for "%s"\n', tag, numel(psas), ks{ki});
    end
end

function combine_param_space(run_dirs, fig_dir)
    psas = {};
    for i = 1:numel(run_dirs)
        listing = dir(fullfile(run_dirs{i}, 'param_space_*'));
        listing = listing([listing.isdir]);
        for k = 1:numel(listing)
            pf = fullfile(listing(k).folder, listing(k).name, 'psa_object.mat');
            if ~exist(pf, 'file'); continue; end
            S = load(pf); fns = fieldnames(S); psas{end+1} = S.(fns{1}); %#ok<AGROW>
        end
    end
    if isempty(psas); return; end
    [ok, reason] = validate_psa_group(psas);
    if ~ok
        warning('combine_runs:incompatible', 'Skipping param_space: %s', reason);
        return;
    end
    base = psas{1}; others = psas(2:end);
    base.plot('metric', 'LLE', 'pool_with', others);
    base.plot('metric', 'mean_rate', 'pool_with', others);
    save_some_figs_to_folder_2(fig_dir, 'combined_param_space', [], {'fig', 'png'});
    close all;
    fprintf('  [param_space] pooled %d runs\n', numel(psas));
end

function combine_dc(run_dirs, fig_dir)
    R = {};
    for i = 1:numel(run_dirs)
        listing = dir(fullfile(run_dirs{i}, 'dc_lle_nSeeds_*'));
        listing = listing([listing.isdir]);
        for k = 1:numel(listing)
            rf = fullfile(listing(k).folder, listing(k).name, 'dc_lle_results.mat');
            if ~exist(rf, 'file'); continue; end
            S = load(rf); R{end+1} = S.dc_lle_results; %#ok<AGROW>
        end
    end
    if isempty(R); return; end

    base = R{1};
    offsets = dc_offset(base);
    pooled = base.per_seed_level_mean;
    for j = 2:numel(R)
        if ~isequaln(R{j}.dc_levels, base.dc_levels)
            warning('combine_runs:DCLevels', 'Skipping DC pooling: dc_levels differ.'); return;
        end
        if ~isequaln(R{j}.config, base.config)
            warning('combine_runs:DCConfig', 'Skipping DC pooling: config structs differ.'); return;
        end
        oj = dc_offset(R{j});
        if ismember(oj, offsets)
            warning('combine_runs:DCOffset', ...
                'Skipping DC pooling: two runs share seed_offset=%g (identical seeds).', oj);
            return;
        end
        offsets(end+1) = oj; %#ok<AGROW>
        pooled = [pooled; R{j}.per_seed_level_mean]; %#ok<AGROW>
    end

    level_mean = mean(pooled, 1)';
    level_std  = std(pooled, 0, 1)';
    n_total = size(pooled, 1);

    figure('Name', 'DC LLE (combined runs)');
    confplot(base.dc_levels, level_mean, level_std, level_std, [0 0 0.8; 0.7 0.8 1.0]);
    yline(0, '--k', 'edge of chaos', 'HandleVisibility', 'off');
    xlabel('DC level (input units)');
    ylabel('mean local \lambda_1  (\pm std across seeds)');
    title(sprintf('DC LLE, combined  (n_{seeds}=%d over %d runs)', n_total, numel(R)));
    grid off;
    save_some_figs_to_folder_2(fig_dir, 'combined_dc_lle', [], {'fig', 'png'});
    close all;
    fprintf('  [dc_lle] pooled %d runs (%d seeds total)\n', numel(R), n_total);
end

function [ok, reason] = validate_psa_group(psas)
    % Same config across all + distinct network_seed_offset.
    ok = true; reason = '';
    if numel(psas) < 2; return; end   % single run: nothing to check
    base = psas{1};
    offsets = base.network_seed_offset;
    for j = 2:numel(psas)
        [tf, r] = base.same_config(psas{j});
        if ~tf
            ok = false; reason = sprintf('config mismatch (run %d): %s', j, r); return;
        end
        offsets(end+1) = psas{j}.network_seed_offset; %#ok<AGROW>
    end
    if numel(unique(offsets)) < numel(offsets)
        ok = false;
        reason = ['runs share a network_seed_offset (identical seeds would ', ...
            'double-count). Give each run a distinct offset (see combine_runs help).'];
    end
end

function off = dc_offset(r)
    if isfield(r, 'seed_offset'); off = r.seed_offset; else; off = 0; end
end

function write_combined_manifest(out_dir, run_dirs, run_info)
    % Save which runs were pooled + their key settings, so a combined_runs
    % folder is self-documenting. Writes both a .mat (machine-readable) and a
    % .txt (human-readable).
    combined_manifest = struct();
    combined_manifest.created = char(datetime('now'));
    combined_manifest.n_runs = numel(run_dirs);
    combined_manifest.source_runs = run_dirs(:);   % cell of paths
    combined_manifest.run_info = run_info;          % per-run provenance struct array
    save(fullfile(out_dir, 'combined_manifest.mat'), 'combined_manifest');

    fid = fopen(fullfile(out_dir, 'combined_manifest.txt'), 'w');
    if fid < 0; return; end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'Combined runs created: %s\n', combined_manifest.created);
    fprintf(fid, 'Number of source runs: %d\n\n', numel(run_dirs));
    for i = 1:numel(run_info)
        fprintf(fid, 'run %d: %s\n', i, run_info(i).name);
        fprintf(fid, '  path: %s\n', run_info(i).dir);
        if run_info(i).has_manifest
            fprintf(fid, '  run_mode=%s  seed_offset=%g  git=%s\n', ...
                run_info(i).run_mode, run_info(i).seed_offset, run_info(i).git_commit_short);
            fprintf(fid, '  activation=%s\n', run_info(i).activation);
        else
            fprintf(fid, '  (no run_manifest.mat -- settings/nonlinearity unverified)\n');
        end
    end
end

function v = getf(s, f, d)
    if isfield(s, f); v = s.(f); else; v = d; end
end

function str = mfval(v)
    if ischar(v); str = v; else; str = mat2str(v); end
end
