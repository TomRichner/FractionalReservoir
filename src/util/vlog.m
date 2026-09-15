function out = vlog(action, arg)
%VLOG The transcript log every vprintf / vfail line is appended to, when one is open.
%
%   guard = VLOG('open', log_file)   % start appending; returns an onCleanup
%   VLOG('append', str)              % append str verbatim (no-op when closed)
%   p     = VLOG('path')             % '' when no log is open
%   VLOG('close')                    % stop (the guard does this)
%
% WHY THIS EXISTS (TR, 2026-09-15). The two entry points used to save their
% transcript with MATLAB's diary, which only sees output that reaches the
% command window: a run launched through the MATLAB MCP server has its
% output captured before diary does, and both command_window.log files came
% out EMPTY for every MCP-launched run (the mu5 run among them). Every stage
% and figure already prints through vprintf, and failures through vfail, so
% the log sink lives there instead: whatever launched MATLAB, the lines land
% in <run_dir>/command_window.log.
%
% What it does not capture: warnings and errors MATLAB raises on its own
% (run_all_paper_analyses appends lastwarn after each stage and the caught
% error on failure), and anything a parfor worker prints (per-job lines,
% 'verbose' only; the batch progress lines come from the client and are
% logged). Bare fprintf calls in library code bypass it; use vprintf.
%
% One log at a time, process-global like diary; opening a second replaces
% the first. The file is opened for append on every write and closed again,
% so a crash never leaves a half-written buffer.
%
% See also: vprintf, vfail, run_all_paper_analyses, make_all_paper_figures

persistent log_path
if isempty(log_path); log_path = ''; end

switch action
    case 'open'
        log_path = char(arg);
        d = fileparts(log_path);
        if ~isempty(d) && ~isfolder(d); mkdir(d); end
        fid = fopen(log_path, 'a');   % touch it so an empty log still means "opened, nothing printed"
        if fid < 0
            error('vlog:CannotOpen', 'Cannot open log file %s for appending.', log_path);
        end
        fclose(fid);
        out = onCleanup(@() vlog('close'));
    case 'append'
        if ~isempty(log_path)
            fid = fopen(log_path, 'a');
            if fid >= 0
                fwrite(fid, char(arg), 'char');
                fclose(fid);
            end
        end
        out = [];
    case 'path'
        out = log_path;
    case 'close'
        log_path = '';
        out = [];
    otherwise
        error('vlog:BadAction', 'Unknown action ''%s'' (open | append | path | close).', action);
end
end
