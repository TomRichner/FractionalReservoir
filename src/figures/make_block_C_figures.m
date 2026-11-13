function out = make_block_C_figures(data, out_dir, cfg)
% MAKE_BLOCK_C_FIGURES Every block C figure and the statistics tables for one run, in one call.
%
%   out = MAKE_BLOCK_C_FIGURES(data, out_dir)
%   out = MAKE_BLOCK_C_FIGURES(data, out_dir, 'ei_files', {f_3to2, f_2to3})
%   out = MAKE_BLOCK_C_FIGURES(data, out_dir, 'sets', {'three', {'no_adaptation', 'sfa3_stdc2'}})
%
% data      the gain-by-bias result, any of: its mc_gain_sweep_results.mat;
%           its mc_gain_sweep\ folder or the run directory above it (the
%           consolidated file when it exists, else the job files, read with
%           read_mc_gain_sweep_jobs, so a run in progress can be drawn); a
%           jobs\ folder; or a struct from read_mc_gain_sweep.
% out_dir   where to write; one subfolder per condition set.
%
% Options:
%   ei_files    results of the gain-only runs on the E:I variants (C3b,
%               C3c), drawn as further rows of fig_mc_gain_curves; each
%               entry takes any form data takes. Default none.
%   ei_labels   their row labels (default from the preset names).
%   sets        name / condition-list pairs, one subfolder each. Default
%               {'three_conditions', {no_adaptation, sfa1_std1, sfa3_stdc2},
%                'four_conditions',  {... , sfa3_std2}} (the paper's three,
%               then the four with product MTS; Plan C section 5).
%   reference   the paired-comparison reference of the statistics
%               (default sfa3_stdc2).
%   visible     figures on screen (default false).
%   verbose     'minimal' (default) | 'verbose' | 'near-none'.
%
% A CONDITION A RESULT DOES NOT HOLD IS LEFT OUT of that set, with a note
% (out.notes, and the console), not an error: the pilot has two of the four
% conditions, and the figures and tables are still wanted for it. A set
% left with no condition is skipped. For the gain curves the conditions
% are those held by EVERY result drawn (the main and the E:I ones).
%
% Per set, in <out_dir>/<set>/: Fig_MC_gain_bias_maps, Fig_MC_gain_curves,
% Fig_MC_robustness (save_figure_stable's formats) and
% MC_gain_bias_statistics.{csv,md} with _by_bias.csv and _by_gain.csv.
%
% out: struct with per-set fields (files, conditions, statistics table S),
% notes (cellstr) and source.
%
% See also: fig_mc_gain_bias_maps, fig_mc_gain_curves, fig_mc_robustness,
%           mc_gain_bias_statistics, read_mc_gain_sweep_jobs

arguments
    data
    out_dir (1,:) char
    cfg.ei_files                  = {}
    cfg.ei_labels                 = {}
    cfg.sets      (1,:) cell      = {'three_conditions', {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'}, ...
                                     'four_conditions',  {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'}}
    cfg.reference (1,:) char      = 'sfa3_stdc2'
    cfg.visible   (1,1) logical   = false
    cfg.verbose   (1,:) char      = 'minimal'
end

setup_paths();
D = load_any(data);
ei = cfg.ei_files; if ~iscell(ei); ei = {ei}; end
D_ei = cellfun(@load_any, ei, 'UniformOutput', false);
if ~isfolder(out_dir); mkdir(out_dir); end

out = struct('source', D.source, 'notes', {{}});
assert(mod(numel(cfg.sets), 2) == 0, 'make_block_C_figures:BadSets', 'sets is name / list pairs.');
for s = 1:2:numel(cfg.sets)
    name = cfg.sets{s}; want = cfg.sets{s + 1};
    have = want(ismember(want, D.conditions));
    if numel(have) < numel(want)
        out.notes{end+1} = sprintf('%s: the result holds no %s; left out.', name, ...
            strjoin(want(~ismember(want, D.conditions)), ', '));
        vprintf(cfg.verbose, 'minimal', '[make_block_C_figures] %s\n', out.notes{end});
    end
    if isempty(have)
        out.notes{end+1} = sprintf('%s: no condition left; set skipped.', name);
        continue
    end
    d = fullfile(out_dir, name);
    if ~isfolder(d); mkdir(d); end
    common = {'save', true, 'visible', cfg.visible, 'out_dir', d, 'verbose', cfg.verbose};
    o = fig_mc_gain_bias_maps('data', D, 'conditions', have, common{:});
    files = o.files; close(o.figs);

    curves_data = [{D}, D_ei];
    have_c = have;
    for k = 1:numel(D_ei); have_c = have_c(ismember(have_c, D_ei{k}.conditions)); end
    if numel(have_c) < numel(have)
        out.notes{end+1} = sprintf('%s: gain curves without %s (not in every E:I result).', name, ...
            strjoin(setdiff(have, have_c, 'stable'), ', '));
    end
    if ~isempty(have_c)
        labels = {};
        if ~isempty(cfg.ei_labels); labels = [{'E : I = 1 : 1'}, cfg.ei_labels]; end
        o = fig_mc_gain_curves('data', curves_data, 'row_labels', labels, 'conditions', have_c, common{:});
        files = [files, o.files]; close(o.figs); %#ok<AGROW>
    end

    o = fig_mc_robustness('data', D, 'conditions', have, common{:});
    files = [files, o.files]; close(o.figs); %#ok<AGROW>

    [S, f_stat] = mc_gain_bias_statistics('data', D, 'conditions', have, 'reference', cfg.reference, ...
        'out_dir', d, 'verbose', cfg.verbose);
    % existing_outputs gives the figures' bare names; make every entry a path
    bare = ~contains(files, filesep);
    files(bare) = fullfile(d, files(bare));
    files = [files, f_stat]; %#ok<AGROW>
    out.(name) = struct('dir', d, 'conditions', {have}, 'files', {files}, 'S', S);
    vprintf(cfg.verbose, 'minimal', '[make_block_C_figures] %s: %d conditions, %d files -> %s\n', ...
        name, numel(have), numel(files), d);
end
end

%% ------------------------------------------------------------------------
function D = load_any(x)
% A struct from read_mc_gain_sweep, a results .mat, or a folder: the
% consolidated file when there is one, else the job files.
if isstruct(x) && isfield(x, 'gains'); D = x; return; end
if istable(x) || (isstruct(x) && isfield(x, 'T')); D = read_mc_gain_sweep(x); return; end
x = char(x);
if isfile(x); D = read_mc_gain_sweep(x); return; end
if ~isfolder(x)
    error('make_block_C_figures:NoData', 'Neither a file nor a folder: %s', x);
end
for c = {fullfile(x, 'mc_gain_sweep', 'mc_gain_sweep_results.mat'), fullfile(x, 'mc_gain_sweep_results.mat')}
    if isfile(c{1}); D = read_mc_gain_sweep(c{1}); return; end
end
D = read_mc_gain_sweep_jobs(x);
end
