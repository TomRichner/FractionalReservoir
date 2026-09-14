function rank = verbose_level(v)
%VERBOSE_LEVEL Rank of a verbosity setting: 'near-none' 0, 'minimal' 1, 'verbose' 2.
%
%   rank = VERBOSE_LEVEL(v)
%
% The code base has ONE verbosity setting with three levels, carried as a char
% from paper_config through ctx / cfg into the model classes and the parfor
% workers (see CLAUDE.md, "Running things"):
%
%   'verbose'   -- everything: per model, per job, per seed. For a human at the
%                  prompt debugging one run.
%   'minimal'   -- THE DEFAULT. One line per stage or batch, one line per
%                  failure, the final summary. What an agent supervising a run
%                  over the MATLAB MCP wants: the transcript of a full pipeline
%                  run fits on a screen per stage.
%   'near-none' -- errors and the final one-line outcome of each entry point.
%
% Logical input is accepted for the properties and arguments that used to be
% logical (ParamSpaceAnalysis2.verbose, run_all_analyses 'verbose', old
% psa_object.mat files): true -> 'verbose', false -> 'minimal'. Anything else
% errors verbose_level:BadLevel.
%
% Warnings and errors are never gated by this setting: quiet means fewer
% lines, never hidden failures. Use vprintf(v, needed, fmt, ...) to print
% only at or above a level.
%
% See also: vprintf, paper_config, resolve_run_context

if islogical(v) || (isnumeric(v) && isscalar(v))
    if v
        rank = 2;
    else
        rank = 1;
    end
    return;
end
if isstring(v), v = char(v); end
if ~ischar(v)
    error('verbose_level:BadLevel', ...
        'verbose must be ''verbose'', ''minimal'' or ''near-none'' (got a %s).', class(v));
end
switch lower(strtrim(v))
    case 'verbose',   rank = 2;
    case 'minimal',   rank = 1;
    case 'near-none', rank = 0;
    otherwise
        error('verbose_level:BadLevel', ...
            'verbose must be ''verbose'', ''minimal'' or ''near-none'' (got ''%s'').', v);
end
end
