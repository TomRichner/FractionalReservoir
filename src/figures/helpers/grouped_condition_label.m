function out = grouped_condition_label(x, kind)
% GROUPED_CONDITION_LABEL Condition labels for the grouped main figures, by NAME.
%
%   labels = GROUPED_CONDITION_LABEL(names, kind)
%   names  = GROUPED_CONDITION_LABEL(titles, 'name_of_title')
%
%   names   cell of condition names (snake_case keys), or one char name.
%   kind    which wording the grouped figures paint:
%     'title'     srnn_condition_titles, as every other figure titles a
%                 condition: 'Multiple-Timescale Adaptation'.
%     'sentence'  the same words in sentence case (the column headings of
%                 groups 2 and 3): 'Multiple-timescale adaptation'. Only a
%                 capital followed by lower-case letters, after a space or a
%                 hyphen, is lowered, so acronyms ('SFA 3\tau') survive.
%     'tick'      sentence case without the trailing word 'adaptation'
%                 except in 'No adaptation' (group 5's margin ticks):
%                 'Single-timescale', 'Multiple-timescale (product)'.
%     'two_line'  the title split before ' Adaptation' (group 5's gain
%                 titles): {'Multiple-Timescale', 'Adaptation (product)'};
%                 'No Adaptation' and titles without the word stay one line.
%                 Returns a cell whose entries are char or 1x2 cellstr.
%     'abbrev'    the memory-horizon ticks of group 6: 'No Adaptation', '1TS',
%                 'MTS' (the chain), 'prodMTS' (the product MTS); any other
%                 name falls back to manuscript_style's condition_short, then
%                 to the name itself.
%   labels  1 x N cell, in the order of names.
%
%   'name_of_title' inverts srnn_condition_titles: given display titles (as a
%   figure painted them; a cell entry or a multi-line title is joined with
%   spaces) it returns the condition names. A title that no name carries, or
%   that more than one name carries, is an error
%   ('grouped_condition_label:UnknownTitle' / ':AmbiguousTitle'): the grouped
%   styles use it to learn which condition an inner figure drew in a panel,
%   and a guess there would be a mislabel.
%
% WHY ONE HELPER. Until FD-006 the grouped figures and their styles carried
% literal three-name lists and literal three-label lists, positionally
% aligned with the columns ({'No adaptation','Single-timescale',
% 'Multiple-timescale'}). With a selection the columns are whatever was
% selected, in its order, so every label is looked up from the condition
% NAME drawn there. The wording variants above are the ones the manuscript's
% grouped figures already use; they are derived from srnn_condition_titles
% so that a retitle there (as sfa3_std2 -> '(product)' on 2026-11-12) moves
% the grouped figures too. 'abbrev' is the one short map kept here, because
% '1TS' / 'MTS' appear nowhere else.
%
% See also: srnn_condition_titles, titles_by_name, manuscript_style,
%           select_conditions, fig_grouped_main

if ischar(x) || isstring(x); x = cellstr(x); end
x = reshape(x, 1, []);

if strcmp(kind, 'name_of_title')
    out = cellfun(@name_of_title, x, 'UniformOutput', false);
    return
end

titles = titles_by_name(x);
switch kind
    case 'title'
        out = titles;
    case 'sentence'
        out = cellfun(@sentence_case, titles, 'UniformOutput', false);
    case 'tick'
        out = cellfun(@tick_of, titles, 'UniformOutput', false);
    case 'two_line'
        out = cellfun(@two_line, titles, 'UniformOutput', false);
    case 'abbrev'
        out = cellfun(@abbrev, x, 'UniformOutput', false);
    otherwise
        error('grouped_condition_label:BadKind', 'Unknown label kind ''%s''.', kind);
end
end

%% ------------------------------------------------------------------------
function s = sentence_case(t)
s = regexprep(t, '(?<=\S[\s-])([A-Z])(?=[a-z])', '${lower($1)}');
end

function s = tick_of(t)
s = sentence_case(t);
if ~startsWith(s, 'No ')
    s = regexprep(s, ' adaptation', '', 'once');
end
end

function s = two_line(t)
k = strfind(t, ' Adaptation');
if isempty(k) || strcmp(t, 'No Adaptation')
    s = t;
else
    s = {t(1:k(1)-1), t(k(1)+1:end)};
end
end

function s = abbrev(name)
known = {'no_adaptation', 'No Adaptation'; ...
         'sfa1_std1',     '1TS'; ...
         'sfa3_stdc2',    'MTS'; ...
         'sfa3_std2',     'prodMTS'};
j = find(strcmp(known(:, 1), name), 1);
if ~isempty(j)
    s = known{j, 2};
    return
end
st = manuscript_style();
if isKey(st.condition_short, name)
    s = st.condition_short(name);
else
    s = name;
end
end

function name = name_of_title(t)
if iscell(t) || isstring(t); t = strjoin(cellstr(t), ' '); end
t = regexprep(char(t), '\\color\[rgb\]\{[^}]*\}', '');   % a label this file's callers coloured
t = strtrim(regexprep(t, '\s+', ' '));
map = srnn_condition_titles();
keys_ = keys(map); vals = values(map);
hit = keys_(strcmp(vals, t));
if isempty(hit)
    error('grouped_condition_label:UnknownTitle', ...
        'No condition is titled ''%s'' in srnn_condition_titles.', t);
elseif numel(hit) > 1
    error('grouped_condition_label:AmbiguousTitle', ...
        'The title ''%s'' belongs to several conditions: %s.', t, strjoin(hit, ', '));
end
name = hit{1};
end
