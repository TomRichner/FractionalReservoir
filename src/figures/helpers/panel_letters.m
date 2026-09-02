function labels = panel_letters(n)
% PANEL_LETTERS Parenthesised panel labels: (a), (b), ... (z), (aa), (ab), ...
%
%   labels = PANEL_LETTERS(28)   % {'(a)', ..., '(z)', '(aa)', '(ab)'}
%
% Every figure that letters its panels used
%
%   arrayfun(@(ch) sprintf('(%c)', ch), char('a' + (0:n-1)), ...)
%
% which walks straight off the end of the alphabet: 'a' + 26 is '{', so panel 27
% is labelled "({)" and panel 28 "(|)". That was invisible while sheets had at
% most 26 panels, and appeared the moment the sweeps went from 4 adaptation
% regimes to 7 -- the mu sheet is 4 parameters x 7 conditions = 28 panels.
%
% It did not error. It wrote "({)" into the figure, and MATLAB only complained
% later, when savefig tried to serialise a text object whose string is not valid
% tex ("Error in state of SceneNode"), and again when openfig read it back. The
% figure still saved; the label was simply wrong.
%
% See also: AddLetters2Plots

arguments
    n (1,1) double {mustBeInteger, mustBeNonnegative}
end

labels = cell(1, n);
for k = 1:n
    labels{k} = sprintf('(%s)', base26(k));
end
end

function s = base26(k)
% 1 -> a, 26 -> z, 27 -> aa, 28 -> ab, ... (spreadsheet-column style).
s = '';
while k > 0
    r = mod(k - 1, 26);
    s = [char('a' + r), s]; %#ok<AGROW>
    k = floor((k - 1) / 26);
end
end
