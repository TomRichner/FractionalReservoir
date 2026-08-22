function cleanup = with_graphics_defaults(varargin)
% WITH_GRAPHICS_DEFAULTS Set graphics-root defaults, restoring them on scope exit.
%
%   cleanup = WITH_GRAPHICS_DEFAULTS('DefaultAxesFontSize', 14, ...)
%   cleanup = WITH_GRAPHICS_DEFAULTS(struct('DefaultTextInterpreter','tex'))
%
% Returns an onCleanup object. Assign it to a variable that STAYS IN SCOPE for
% as long as the defaults are wanted; when it is cleared or the function
% returns, the previous values are put back -- including on the error path.
%
% WHY THIS IS NEEDED. `set(groot, 'DefaultX', v)` is process-global and
% permanent. A plotting function that sets defaults and does not restore them
% changes every figure drawn afterwards in that MATLAB session, including
% figures belonging to entirely different scripts.
%
% That is not hypothetical here. scripts/memory_capacity/plot_memory_capacity.m
% set DefaultTextInterpreter = 'none' and never restored it, so after any
% memory-capacity plot the session rendered '\lambda_1' as literal backslash
% text -- which is exactly what the sensitivity and param-space sheets use for
% their axis labels. Verified live: immediately after a memory-capacity run,
% get(groot,'DefaultTextInterpreter') returned 'none' against a factory 'tex'.
% Because it depends on the ORDER figures are drawn in, it survives any amount
% of testing one figure at a time and only appears once they are batched.
%
% Prefer passing style explicitly to each object (see manuscript_style). Use
% this only for plotters that read the root defaults instead of taking
% arguments.
%
% See also: manuscript_style, with_manuscript_defaults, onCleanup

if numel(varargin) == 1 && isstruct(varargin{1})
    s = varargin{1};
    names = fieldnames(s)';
    vals  = struct2cell(s)';
else
    if mod(numel(varargin), 2) ~= 0
        error('with_graphics_defaults:BadArgs', ...
            'Expected name/value pairs or a single struct.');
    end
    names = varargin(1:2:end);
    vals  = varargin(2:2:end);
end

prev = cell(size(names));
applied = false(size(names));
for k = 1:numel(names)
    try
        prev{k} = get(groot, names{k});
        set(groot, names{k}, vals{k});
        applied(k) = true;
    catch ME
        warning('with_graphics_defaults:CannotSet', ...
            'Skipping ''%s'': %s', names{k}, ME.message);
    end
end

cleanup = onCleanup(@() restore_defaults(names(applied), prev(applied)));
end

function restore_defaults(names, prev)
for k = 1:numel(names)
    try
        set(groot, names{k}, prev{k});
    catch
        % Restoration runs on the error path too; an exception here would mask
        % whatever actually went wrong, so it is deliberately swallowed.
    end
end
end
