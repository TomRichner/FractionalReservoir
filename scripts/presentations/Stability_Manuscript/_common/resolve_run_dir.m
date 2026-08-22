function run_dir = resolve_run_dir(opts)
% RESOLVE_RUN_DIR Locate the run_all_<dt> directory a figure should plot.
%
%   run_dir = RESOLVE_RUN_DIR()                         % newest run, any preset
%   run_dir = RESOLVE_RUN_DIR('preset_name', p)         % newest run of preset p
%   run_dir = RESOLVE_RUN_DIR('run_dir', d)             % validate an explicit one
%
% THE PROBLEM THIS SOLVES. Five figure scripts hardcoded a data_root, and they
% DID NOT AGREE: fig_EI_param_space pointed at run_all_jul_06_26_22_00 (an old
% SRNNModel2 run) while the four allStd/medians scripts pointed at
% run_all_aug_14_26_17_25 (a SRNNCellTypePairs run). "Regenerate the figures"
% therefore built the manuscript out of two different runs on two different
% model classes, silently. One resolution, passed down from the master script,
% makes that impossible.
%
% ERRORS LOUDLY when nothing matches, rather than falling back to "newest
% anything" -- plotting last month's preset without saying so is exactly the
% failure this replaces.
%
% Selection is by run_manifest.mat, not by folder name: the name carries only a
% timestamp, while the manifest records the preset and the model class that
% actually ran.
%
% See also: srnn_param_preset, preset_default_values, write_run_parameters_md

arguments
    opts.run_dir     (1,:) char = ''
    opts.preset_name (1,:) char = ''
    opts.search_root (1,:) char = ''
    opts.require_manifest (1,1) logical = true
end

% --- An explicit directory: validate and return it -------------------------
if ~isempty(opts.run_dir)
    if ~isfolder(opts.run_dir)
        error('resolve_run_dir:NoSuchDir', ...
            'run_dir does not exist:\n  %s', opts.run_dir);
    end
    run_dir = opts.run_dir;
    if ~isempty(opts.preset_name)
        actual = manifest_preset(run_dir);
        if ~isempty(actual) && ~strcmp(actual, opts.preset_name)
            warning('resolve_run_dir:PresetMismatch', ...
                ['run_dir was built from preset ''%s'', but ''%s'' was requested.\n' ...
                 '  %s'], actual, opts.preset_name, run_dir);
        end
    end
    return
end

% --- Otherwise search ------------------------------------------------------
search_root = opts.search_root;
if isempty(search_root)
    search_root = fullfile(fileparts(which('setup_paths')), 'data', 'param_space');
end
if ~isfolder(search_root)
    error('resolve_run_dir:NoSearchRoot', ...
        'No run directory search root at:\n  %s', search_root);
end

listing = dir(fullfile(search_root, 'run_all_*'));
listing = listing([listing.isdir]);
if isempty(listing)
    error('resolve_run_dir:NoRuns', ...
        'No run_all_* directories under:\n  %s', search_root);
end

% Keep only runs that carry a manifest naming the requested preset.
keep    = false(1, numel(listing));
presets = cell(1, numel(listing));
for k = 1:numel(listing)
    d = fullfile(listing(k).folder, listing(k).name);
    presets{k} = manifest_preset(d);
    if isempty(presets{k})
        keep(k) = ~opts.require_manifest && isempty(opts.preset_name);
    elseif isempty(opts.preset_name)
        keep(k) = true;
    else
        keep(k) = strcmp(presets{k}, opts.preset_name);
    end
end

if ~any(keep)
    available = unique(presets(~cellfun(@isempty, presets)));
    if isempty(available)
        detail = '    (no run carries a run_manifest.mat)';
    else
        detail = sprintf('    %s\n', available{:});
    end
    error('resolve_run_dir:NoMatchingRun', ...
        ['No run under\n  %s\nwas built from preset ''%s''.\n' ...
         '  Presets available there:\n%s\n' ...
         '  Run the analyses for that preset first, or pass an explicit run_dir.'], ...
        search_root, opts.preset_name, detail);
end

listing = listing(keep);
[~, newest] = max([listing.datenum]);
run_dir = fullfile(listing(newest).folder, listing(newest).name);
end

%% ------------------------------------------------------------------------
function name = manifest_preset(run_dir)
% The preset a run recorded, or '' if it has no readable manifest.
name = '';
f = fullfile(run_dir, 'run_manifest.mat');
if ~isfile(f); return; end
try
    M = load(f, 'run_manifest');
    if isfield(M, 'run_manifest') && isfield(M.run_manifest, 'preset_name')
        name = M.run_manifest.preset_name;
    end
catch
    % Unreadable manifest is treated as absent; the caller reports it.
end
end
