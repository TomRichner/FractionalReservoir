function fs = stage_fs(fs_opt, default_fs, who)
% STAGE_FS The integration rate a stage runs at: the caller's, or the stage's own.
%
%   fs = STAGE_FS(fs_opt, default_fs, who)
%
% fs_opt is the 'fs' option a stage was given. 0 or empty means "the stage's
% own default" (default_fs, from its run-mode table or constant), so a stage
% called without the option behaves exactly as before the option existed.
% Anything else must be a finite positive real scalar and is returned as is.
% who names the caller for the error identifier (<who>:BadFs).
%
% Why (FD-031 item 14, Tom 2026-11-17): the bundle config states the run's
% integration rate once (cfg.fs, 500 Hz in mu7chain_config), and
% run_all_paper_analyses hands it to every stage that takes it. Before that
% each stage hard-coded its own rate and cfg.fs was read by nothing.
%
% See also: run_all_paper_analyses, resolve_run_context

if isempty(fs_opt) || (isnumeric(fs_opt) && isscalar(fs_opt) && fs_opt == 0)
    fs = default_fs;
    return
end
if ~(isnumeric(fs_opt) && isscalar(fs_opt) && isreal(fs_opt) && isfinite(fs_opt) && fs_opt > 0)
    error([who ':BadFs'], 'fs must be 0 or empty (the stage''s own rate) or a finite positive scalar.');
end
fs = double(fs_opt);
end
