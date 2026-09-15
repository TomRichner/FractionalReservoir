function vfail(fmt, varargin)
%VFAIL An ungated failure line: stderr AND the transcript log (vlog), when one is open.
%
%   VFAIL(fmt, ...)
%
% The entry points used fprintf(2, ...) for failure lines so they were never
% hidden by the verbose setting; those lines were also the ones a transcript
% most needs, and a log written through diary lost them under the MCP server
% (see vlog). vfail keeps the stderr write and appends the same text to the
% log. Never gated.
%
% See also: vprintf, vlog

s = sprintf(fmt, varargin{:});
fprintf(2, '%s', s);
vlog('append', s);
end
