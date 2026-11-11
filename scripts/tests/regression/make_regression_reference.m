function out_file = make_regression_reference(tier, opts)
% MAKE_REGRESSION_REFERENCE Save the reference digest for one tier of the regression.
%
%   out_file = make_regression_reference('rapid')
%   out_file = make_regression_reference('full')                    % from the run tagged 'reference'
%   out_file = make_regression_reference('full', 'tag', 'reference', 'ignore', patterns, 'notes', text)
%
% TIERS
%   'rapid'  runs regression_rapid_run now (about 15 s), digests it, and
%            writes scripts/tests/fixtures/regression_rapid_digest.mat.
%   'full'   digests an EXISTING run of the full bundle,
%            data/regression_fast/<tag> (and figs/regression_fast/<tag> when
%            it exists), and writes .../fixtures/regression_fast_digest.mat.
%            It does not run the bundle; regression_fast_run does, and takes
%            about 50 minutes.
%
% REGENERATING IS A DELIBERATE ACT, as for the preset fixture. There is no
% overwrite: to replace a reference, delete the .mat and run this again, and
% say in the commit message WHY the expected values moved. A reference that
% is regenerated to make a check pass is worse than no reference.
%
% MAKE IT FROM A CLEAN, COMMITTED TREE. The commit hash is stored with the
% digest; if tracked files have uncommitted changes this refuses, because the
% reference would then describe code that no commit holds.
%
% IGNORE is a cell of regular expressions for keys known to differ between
% two runs of identical code. Establish it by running the tier twice on the
% same code and comparing; it is stored in the reference so every later
% comparison applies it. NOTES is free text stored alongside (how
% reproducible the run was found to be, which commit a 'full' run was
% produced at).
%
% The .mat needs `git add -f` (*.mat is gitignored repo-wide).
%
% See also: regression_rapid_run, regression_fast_run, regression_digest,
%           regression_rapid_check, regression_check
arguments
    tier (1,:) char {mustBeMember(tier, {'rapid', 'full'})}
    opts.tag (1,:) char = 'reference'
    opts.ignore cell = {}
    opts.notes (1,:) char = ''
    opts.allow_dirty (1,1) logical = false
end
project_root = fileparts(which('setup_paths'));
names = struct('rapid', 'regression_rapid_digest.mat', 'full', 'regression_fast_digest.mat');
out_file = fullfile(project_root, 'scripts', 'tests', 'fixtures', names.(tier));
if isfile(out_file)
    error('make_regression_reference:Exists', ...
        ['%s already exists.\nReplacing the reference is deliberate: delete the file, ' ...
         'run this again, and record in the commit message why the expected values moved.'], out_file);
end

[status, commit] = system(sprintf('git -C "%s" rev-parse HEAD', project_root));
if status ~= 0; error('make_regression_reference:NoGit', 'git rev-parse failed: %s', commit); end
commit = strtrim(commit);
% Uncommitted changes to tracked files; untracked files (the run itself, this
% fixture) do not count.
[~, dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', project_root));
dirty = strtrim(dirty);
if ~isempty(dirty) && ~opts.allow_dirty
    error('make_regression_reference:DirtyTree', ...
        'Tracked files have uncommitted changes; commit them first.\n%s', dirty);
end

switch tier
    case 'rapid'
        run_dir = regression_rapid_run();
        digest = regression_digest(run_dir);
        rmdir(run_dir, 's');
        preset_name = 'regression_pairs_noStim_noise0p025_dualStd_3cond';
        run_mode = 'none (direct calls)';
        tag = '';
    case 'full'
        tag = opts.tag;
        cfg = regression_fast_config(tag);
        run_dir = fullfile(project_root, cfg.run_dir);
        fig_root = fullfile(project_root, cfg.fig_root);
        if isfolder(fig_root)
            digest = regression_digest(run_dir, 'fig_root', fig_root);
        else
            digest = regression_digest(run_dir);
        end
        preset_name = cfg.preset_name;
        run_mode = cfg.run_mode;
end
digest.ignore = opts.ignore(:)';

made_on = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
notes = opts.notes;
matlab_release = version;
machine = computer;
save(out_file, 'digest', 'commit', 'made_on', 'notes', 'matlab_release', 'machine', ...
    'preset_name', 'run_mode', 'tag', 'tier', '-v7');
d = dir(out_file);
fprintf('Wrote %s\n  %d leaves, %d files, %d figure files, %d ignore pattern(s), %.0f kB\n', ...
    out_file, numel(digest.keys), numel(digest.files), numel(digest.fig_files), ...
    numel(digest.ignore), d.bytes / 1024);
fprintf('  commit %s, %s\n', commit, made_on);
fprintf('  *.mat is gitignored -- commit with:  git add -f %s\n', strrep(out_file, '\', '/'));
end
