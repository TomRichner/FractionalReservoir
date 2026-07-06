function info = capture_git_provenance(output_dir, repo_root)
%CAPTURE_GIT_PROVENANCE Save git state (commit, branch, dirty diff) to output_dir.
%
% Writes git_provenance.txt with full commit hash, branch, remote, dirty
% file list, and (if the working tree is dirty) working_changes.patch
% containing the unstaged+staged diff against HEAD plus a list of any
% untracked files. Together these let a future run check out the exact
% commit and re-apply the working changes to reproduce results.
%
% Usage:
%   capture_git_provenance(master_output_dir, project_root)
%   info = capture_git_provenance(master_output_dir, project_root)
%
% Optional output `info` is a struct with fields commit, commit_short, branch,
% is_dirty (empty strings / false when not a git working tree) so callers can
% record the git state machine-readably without re-shelling out.
%
% If git is not available or repo_root is not a git working tree, a
% warning is issued and a stub file is written; the caller is not
% interrupted.

    if nargin < 2 || isempty(repo_root)
        repo_root = pwd;
    end

    info = struct('commit', '', 'commit_short', '', 'branch', '', 'is_dirty', false);

    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    % --no-pager is essential: when called via system(), `git log`/`git diff`
    % would otherwise spawn a pager (less) that blocks waiting for input,
    % hanging the run. rev-parse/config/status never paginate, which is why
    % only the log call hung.
    git = sprintf('git -C "%s" --no-pager', repo_root);
    prov_path = fullfile(output_dir, 'git_provenance.txt');

    [s_check, ~] = system([git ' rev-parse --is-inside-work-tree']);
    if s_check ~= 0
        warning('capture_git_provenance:NotARepo', ...
            'Not a git repo at %s; writing stub provenance file.', repo_root);
        fid = fopen(prov_path, 'w');
        fprintf(fid, 'Not a git working tree: %s\nTimestamp: %s\n', ...
            repo_root, char(datetime('now')));
        fclose(fid);
        return;
    end

    commit       = run_git(git, 'rev-parse HEAD');
    commit_short = run_git(git, 'rev-parse --short HEAD');
    branch       = run_git(git, 'rev-parse --abbrev-ref HEAD');
    remote       = run_git(git, 'config --get remote.origin.url');
    last_msg     = run_git(git, 'log -1 --pretty=%s');
    last_date    = run_git(git, 'log -1 --pretty=%cI');
    porcelain    = run_git(git, 'status --porcelain');
    untracked    = run_git(git, 'ls-files --others --exclude-standard');

    is_dirty = ~isempty(strtrim(porcelain));

    info.commit = commit;
    info.commit_short = commit_short;
    info.branch = branch;
    info.is_dirty = is_dirty;

    fid = fopen(prov_path, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'Git provenance for run output: %s\n', output_dir);
    fprintf(fid, 'Captured at: %s\n\n', char(datetime('now')));
    fprintf(fid, 'commit:      %s\n', commit);
    fprintf(fid, 'commit_short:%s\n', commit_short);
    fprintf(fid, 'branch:      %s\n', branch);
    fprintf(fid, 'remote:      %s\n', remote);
    fprintf(fid, 'last_msg:    %s\n', last_msg);
    fprintf(fid, 'last_date:   %s\n', last_date);
    fprintf(fid, 'dirty:       %s\n\n', mat2str(is_dirty));
    if is_dirty
        fprintf(fid, '--- git status --porcelain ---\n%s\n', porcelain);
        if ~isempty(strtrim(untracked))
            fprintf(fid, '\n--- untracked files (not in patch) ---\n%s\n', untracked);
        end
    end
    clear cleanup;  % flush

    if is_dirty
        patch_path = fullfile(output_dir, 'working_changes.patch');
        cmd = sprintf('%s diff HEAD > "%s"', git, patch_path);
        [s_diff, ~] = system(cmd);
        if s_diff ~= 0
            warning('capture_git_provenance:DiffFailed', ...
                'git diff HEAD failed (status %d); patch not saved.', s_diff);
        end
    end

    fprintf('Git provenance saved: %s @ %s%s\n', branch, commit_short, ...
        ternary(is_dirty, ' (DIRTY — see working_changes.patch)', ''));
end

function out = run_git(git, args)
    [status, raw] = system([git ' ' args]);
    if status ~= 0
        out = '';
    else
        out = strtrim(raw);
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
