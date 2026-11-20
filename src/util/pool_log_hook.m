function pool_log_hook(str)
%POOL_LOG_HOOK Append a line the shared pool module printed to this repo's transcript log.
%
%   POOL_LOG_HOOK(str)
%
% The pool module (pool_set, pool_ensure, pool_restart, pool_release,
% pool_status) lives in the shared MatlabAgenticOrg repository since
% 2026-11-19 and prints through its own helper, not through vprintf. It
% calls a function of this name, when one is on the path, with every line
% it prints. Passing the line to vlog keeps the pool's lines in
% <run_dir>/command_window.log, where they were while the module was part
% of this repository.
%
% See also: vlog, vprintf, use_matlab_agentic_org
vlog('append', str);
end
