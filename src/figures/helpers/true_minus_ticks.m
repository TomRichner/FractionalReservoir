function info = true_minus_ticks(fig)
% TRUE_MINUS_TICKS Write every numeric tick label of a finished figure with the true minus sign.
%
%   info = true_minus_ticks(fig)
%
% MATLAB writes its own tick labels with the hyphen-minus ('-', U+002D),
% while labels a figure sets by hand were written with the true minus
% (U+2212) in some places, so one figure could show both glyphs (Tom,
% 2026-11-23: "probably should be minus"). This rewrites the numeric tick
% labels of every axes ruler and color bar of FIG with U+2212 in place of
% the hyphen, the minus of an exponent included ('10^{-2}').
%
% It is the last step before a figure is saved (composite_save,
% fig_grouped_main and the composites that call save_figure_stable
% themselves): labels it sets are fixed, so a layout step after it could
% leave them out of step with their ticks. What it does not change:
%   - which ticks are shown: where it writes the labels of a ruler it also
%     fixes that ruler's tick values at the values it has now, so a later
%     resize cannot move them;
%   - a label that is not a number, or a number in a power of ten, such as
%     a condition name, "E-to-E x 1.5" or a date: only labels made of an
%     optional sign, digits, a decimal point, an exponent (1e-3, 10^{-3},
%     2\times10^{-3}) and an optional percent sign are touched;
%   - a ruler with the 'latex' interpreter (LaTeX sets its own minus);
%   - a ruler whose labels MATLAB writes itself with a common exponent
%     (x 10^{-3} beside the axis): fixed labels would lose that exponent,
%     so such a ruler is left as it is and listed in info.skipped.
% A ruler with no hyphen in any numeric label is not touched at all, so a
% figure without negative ticks is exported exactly as before.
%
% info: struct with fields changed (number of rulers rewritten) and
% skipped (cell of descriptions of rulers left because of an exponent).

minus = char(8722);
info = struct('changed', 0, 'skipped', {{}});
if isempty(fig) || ~isgraphics(fig); return; end

axs = findall(fig, 'Type', 'axes');
for k = 1:numel(axs)
    ax = axs(k);
    for p = {'XAxis', 'YAxis', 'ZAxis'}
        rulers = ax.(p{1});
        for r = 1:numel(rulers)                    % yyaxis: two y rulers
            info = fix_ruler(rulers(r), 'TickValues', info, minus, sprintf('%s of axes %d', p{1}, k));
        end
    end
end
pax = findall(fig, 'Type', 'polaraxes');
for k = 1:numel(pax)
    for p = {'ThetaAxis', 'RAxis'}
        info = fix_ruler(pax(k).(p{1}), 'TickValues', info, minus, sprintf('%s of polar axes %d', p{1}, k));
    end
end
cbs = findall(fig, 'Type', 'colorbar');
for k = 1:numel(cbs)
    info = fix_ruler(cbs(k), 'Ticks', info, minus, sprintf('color bar %d', k));
end
end

function info = fix_ruler(r, tick_prop, info, minus, what)
if strcmpi(r.TickLabelInterpreter, 'latex'); return; end
labels = r.TickLabels;
if isempty(labels); return; end
was_char = ischar(labels);
labels = cellstr(labels);
new = labels;
for i = 1:numel(labels)
    if contains(labels{i}, '-') && is_numeric_label(labels{i})
        new{i} = strrep(labels{i}, '-', minus);
    end
end
if isequal(new, labels); return; end
if strcmp(r.TickLabelsMode, 'auto') && has_exponent(r)
    info.skipped{end + 1} = what;
    return
end
if strcmp(r.([tick_prop 'Mode']), 'auto')
    r.(tick_prop) = r.(tick_prop);                 % fix the ticks where they are
end
if was_char; new = char(new); end
r.TickLabels = new;
info.changed = info.changed + 1;
end

function tf = is_numeric_label(s)
s = strtrim(s);
num = '(\d+\.?\d*|\.\d+)';
plain = ['^[-+]?' num '([eE][-+]?\d+)?\s*%?$'];
power = ['^[-+]?(' num '\s*(\\times|\\cdot|' char(215) '|' char(183) ')\s*)?10\^(\{\s*[-+]?\d+\s*\}|[-+]?\d)$'];
tf = ~isempty(regexp(s, plain, 'once')) || ~isempty(regexp(s, power, 'once'));
end

function tf = has_exponent(r)
tf = false;
if isprop(r, 'Exponent')
    tf = r.Exponent ~= 0;
elseif isprop(r, 'Ruler') && isprop(r.Ruler, 'Exponent')   % color bar
    tf = r.Ruler.Exponent ~= 0;
end
end
