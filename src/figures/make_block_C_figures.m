function out = make_block_C_figures(data, out_dir, cfg)
% MAKE_BLOCK_C_FIGURES Every block C figure and the statistics tables for one result, in one call.
%
%   out = MAKE_BLOCK_C_FIGURES(data, out_dir)
%   out = MAKE_BLOCK_C_FIGURES(data, out_dir, 'variants', {f_3to2, f_2to3})
%   out = MAKE_BLOCK_C_FIGURES(data, out_dir, 'sets', {'three', {'no_adaptation', 'sfa3_stdc2'}})
%
% For looking at one result outside a bundle, a run still in progress
% included, with the report's color scales (square-root memory capacity,
% symmetric-log lambda_1). A bundle draws its own selection, per grid of
% its run directory and with linear color scales, through
% make_all_paper_figures and its registry (mu7chain_gainbias_config).
%
% data      the result of a gain sweep (against the bias, or whatever second
%           axis the sweep ran; the figures and tables name it), any of: its mc_gain_sweep_results.mat;
%           its mc_gain_sweep folder or the run directory above it (the
%           consolidated file when it exists, else the job files, read with
%           read_mc_gain_sweep_jobs); a jobs folder; or anything
%           read_mc_gain_sweep takes.
% out_dir   where to write; one subfolder per condition set.
%
% Options:
%   variants    results of the gain-only runs on variant networks, drawn as
%               further rows of fig_mc_gain_curves; each entry takes any
%               form data takes. Default none.
%   row_labels  that figure's row labels, the main result's first (default:
%               fig_mc_gain_curves' own, from each result's preset).
%   sets        name / condition-list pairs, one subfolder each. Default:
%               the paper's three conditions, then the four with product MTS.
%   reference   the paired-comparison reference of the statistics
%               (default sfa3_stdc2).
%   visible     figures on screen (default false).
%   verbose     'minimal' (default) | 'verbose' | 'near-none'.
%
% A CONDITION A RESULT DOES NOT HOLD IS LEFT OUT of that set, with a note
% (out.notes, and the console), not an error, so a result with fewer
% conditions than the sets name is still drawn. A set left with no
% condition is skipped. For the gain curves the conditions are those held
% by EVERY result drawn (the main and the variants).
%
% Per set, in <out_dir>/<set>/: Fig_MC_gain_bias_maps, Fig_MC_gain_curves,
% Fig_MC_robustness, Fig_MC_vs_lambda1 (save_figure_stable's formats) and
% MC_gain_bias_statistics.{csv,md} with _by_bias.csv (named by the second
% axis: _by_self_feedback.csv for a self-feedback sweep) and _by_gain.csv.
%
% out: source, notes (cellstr), and one field per set that was drawn (dir,
% conditions, files, the statistics table S).
%
% See also: fig_mc_gain_bias_maps, fig_mc_gain_curves, fig_mc_robustness,
%           fig_mc_vs_lambda1, mc_gain_bias_statistics, read_mc_gain_sweep_jobs

arguments
    data
    out_dir (1,:) char
    cfg.variants   (1,:) cell    = {}
    cfg.row_labels (1,:) cell    = {}
    cfg.sets       (1,:) cell    = {'three_conditions', {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}, ...
                                    'four_conditions',  {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'}}
    cfg.reference  (1,:) char    = 'sfa3_stdc2'
    cfg.visible    (1,1) logical = false
    cfg.verbose    (1,:) char    = 'minimal'
end

setup_paths();
if mod(numel(cfg.sets), 2) ~= 0
    error('make_block_C_figures:BadSets', 'sets is name / condition-list pairs.');
end
D = load_result(data);
D_variants = cellfun(@load_result, cfg.variants, 'UniformOutput', false);
if ~isfolder(out_dir); mkdir(out_dir); end

out = struct('source', D.source, 'notes', {{}});
for s = 1:2:numel(cfg.sets)
    name = cfg.sets{s};
    wanted = cfg.sets{s + 1};
    conds = wanted(ismember(wanted, D.conditions));
    if numel(conds) < numel(wanted)
        out.notes{end+1} = sprintf('%s: the result holds no %s; left out.', name, ...
            strjoin(setdiff(wanted, conds, 'stable'), ', '));
        vprintf(cfg.verbose, 'minimal', '[make_block_C_figures] %s\n', out.notes{end});
    end
    if isempty(conds)
        out.notes{end+1} = sprintf('%s: no condition left; set skipped.', name);
        continue
    end
    set_dir = fullfile(out_dir, name);
    if ~isfolder(set_dir); mkdir(set_dir); end
    common = {'save', true, 'visible', cfg.visible, 'out_dir', set_dir, 'verbose', cfg.verbose};

    files = drawn(fig_mc_gain_bias_maps('data', D, 'conditions', conds, common{:}), set_dir);

    conds_curves = conds;
    for k = 1:numel(D_variants)
        conds_curves = conds_curves(ismember(conds_curves, D_variants{k}.conditions));
    end
    if numel(conds_curves) < numel(conds)
        out.notes{end+1} = sprintf('%s: gain curves without %s (not in every variant''s result).', ...
            name, strjoin(setdiff(conds, conds_curves, 'stable'), ', '));
    end
    if ~isempty(conds_curves)
        files = [files, drawn(fig_mc_gain_curves('data', [{D}, D_variants], ...
            'row_labels', cfg.row_labels, 'conditions', conds_curves, common{:}), set_dir)]; %#ok<AGROW>
    end

    files = [files, drawn(fig_mc_robustness('data', D, 'conditions', conds, common{:}), set_dir)]; %#ok<AGROW>
    files = [files, drawn(fig_mc_vs_lambda1('data', D, 'conditions', conds, common{:}), set_dir)]; %#ok<AGROW>

    [S, stat_files] = mc_gain_bias_statistics('data', D, 'conditions', conds, ...
        'reference', cfg.reference, 'out_dir', set_dir, 'verbose', cfg.verbose);
    files = [files, stat_files]; %#ok<AGROW>
    out.(name) = struct('dir', set_dir, 'conditions', {conds}, 'files', {files}, 'S', S);
    vprintf(cfg.verbose, 'minimal', '[make_block_C_figures] %s: %d conditions, %d files -> %s\n', ...
        name, numel(conds), numel(files), set_dir);
end
end

%% ------------------------------------------------------------------------
function files = drawn(fig_out, set_dir)
% A figure's saved files as full paths (existing_outputs gives bare names),
% and its window closed.
files = fig_out.files;
bare = ~contains(files, filesep);
files(bare) = fullfile(set_dir, files(bare));
close(fig_out.figs);
end

function D = load_result(x)
% A result in memory, a results .mat, or a folder: the consolidated file
% when there is one, else the job files.
if ~(ischar(x) || isstring(x))
    D = read_mc_gain_sweep(x);
    return
end
x = char(x);
if isfile(x)
    D = read_mc_gain_sweep(x);
    return
end
if ~isfolder(x)
    error('make_block_C_figures:NoData', 'Neither a file nor a folder: %s', x);
end
candidates = {fullfile(x, 'mc_gain_sweep', 'mc_gain_sweep_results.mat'), ...
              fullfile(x, 'mc_gain_sweep_results.mat')};
found = candidates(cellfun(@isfile, candidates));
if ~isempty(found)
    D = read_mc_gain_sweep(found{1});
else
    D = read_mc_gain_sweep_jobs(x);
end
end
