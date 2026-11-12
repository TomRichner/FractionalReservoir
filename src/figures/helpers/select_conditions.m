function [names, idx] = select_conditions(available, requested)
% SELECT_CONDITIONS The conditions a figure draws, in the order it draws them.
%
%   [names, idx] = SELECT_CONDITIONS(available, requested)
%
%   available  cell of condition names the run (or preset) holds, in its own
%              order. A char or string names one condition.
%   requested  cell (or string array) of condition names to draw, or empty
%              ({}, '', "", []) meaning "all of them". A char or string names
%              one condition.
%
%   names      1 x N cellstr, the names to draw, left to right.
%   idx        1 x N double, their positions in `available`:
%              names == available(idx) always holds.
%
% THE RULES (fixed for every figure; a figure never re-implements them):
%
%   1. EMPTY REQUEST = EVERYTHING, UNCHANGED. An empty `requested` returns
%      `available` in its own order and idx = 1:numel(available), so a bundle
%      that sets no cfg.conditions draws exactly what it drew before the option
%      existed.
%   2. REQUESTED ORDER WINS. Otherwise the result follows the order of
%      `requested`, not of `available`: that order is the left-to-right order
%      of a figure's columns (and the order of overlaid curves and legend
%      entries). It is what lets the chained MTS take the product MTS's slot.
%   3. A NAME THE RUN DOES NOT HOLD IS AN ERROR,
%      'select_conditions:NoSuchCondition', whose message lists what is
%      available. A figure must never silently draw fewer, or other,
%      conditions than it was asked for (the same reason resolve_data_file
%      has no fallback).
%   4. A DUPLICATE IN `requested` IS AN ERROR,
%      'select_conditions:DuplicateCondition'. Drawing a condition twice is
%      never what was meant.
%   5. Names are matched exactly (case-sensitive), as everywhere else that
%      keys on condition names (srnn_condition_titles, manuscript_style).
%
% Malformed input (a non-text entry, a duplicate in `available`) raises
% 'select_conditions:BadInput'.
%
% HOW A FIGURE USES IT (the pattern; see fig_lle_vs_rate for a worked one):
%
%   arguments
%       ...
%       cfg.conditions = {}   % names to draw, in order; empty = all the run has
%   end
%   ...
%   available  = <the condition names the data holds, in its order>;
%   [cond_names, idx] = select_conditions(available, cfg.conditions);
%   % index every per-condition array by idx, or key by cond_names; take
%   % titles and colours from the NAME (manuscript_style / srnn_condition_titles),
%   % never from a saved positional title list.
%
% make_all_paper_figures passes cfg.conditions as the name-value pair
% 'conditions' to every registry entry that declares it.
%
% See also: make_all_paper_figures, fig_lle_vs_rate, srnn_condition_titles,
%           manuscript_style, resolve_data_file

available = as_names(available, 'available');
if numel(unique(available)) ~= numel(available)
    error('select_conditions:BadInput', ...
        'available holds a condition name twice: {%s}.', strjoin(available, ', '));
end

if is_empty_request(requested)
    names = available;
    idx   = 1:numel(available);
    return
end

requested = as_names(requested, 'requested');

[u, ~, j] = unique(requested, 'stable');
if numel(u) ~= numel(requested)
    counts = accumarray(j(:), 1);
    error('select_conditions:DuplicateCondition', ...
        'A condition is requested more than once (%s); requested: {%s}.', ...
        strjoin(u(counts > 1), ', '), strjoin(requested, ', '));
end

[found, idx] = ismember(requested, available);
if ~all(found)
    error('select_conditions:NoSuchCondition', ...
        ['Requested condition(s) not held by this run: {%s}.\n' ...
         'Available, in the run''s order: {%s}.'], ...
        strjoin(requested(~found), ', '), strjoin(available, ', '));
end

names = requested;
idx   = reshape(idx, 1, []);
end

%% ------------------------------------------------------------------------
function tf = is_empty_request(r)
% {}, '', "", [] and a 0-element string array all mean "everything".
tf = isempty(r) || (isstring(r) && all(strlength(r) == 0) && isscalar(r));
end

function c = as_names(x, what)
% Normalise to a 1 x N cellstr of char rows.
if ischar(x)
    c = {x};
elseif isstring(x)
    c = cellstr(x(:)');
elseif iscell(x)
    c = x(:)';
    for k = 1:numel(c)
        if isstring(c{k}) && isscalar(c{k}); c{k} = char(c{k}); end
        if ~(ischar(c{k}) && (isrow(c{k}) || isempty(c{k})))
            error('select_conditions:BadInput', ...
                '%s{%d} is not a condition name (class %s).', what, k, class(c{k}));
        end
    end
elseif isempty(x)
    c = {};
else
    error('select_conditions:BadInput', ...
        '%s must be a cell of condition names, a char or a string (got %s).', what, class(x));
end
c = reshape(c, 1, []);
end
