function apply_percent_axis(ax, d, base_label, label_fs)
% APPLY_PERCENT_AXIS Relabel an x-axis as percent departure from a default.
%
%   APPLY_PERCENT_AXIS(ax, d, base_label, label_fs)
%
% Ticks are placed at data positions d*(1 + p/100) and labelled p%, so the
% underlying plotted data is untouched -- only the ruler changes. 0% is always
% included, since "the preset's own network" is the reference the rest is read
% against.
%
% For a NEGATIVE default (the inhibitory mu blocks) increasing percent means a
% more negative value, i.e. leftward in data coordinates. XDir is reversed there
% so percent still ascends left-to-right and all the mu panels read alike.
%
% Extracted verbatim from the local subfunction of
% Fig_sensitivity_analysis_allStd.m so fig_sensitivity_medians can share it.
%
% See also: preset_default_values, mark_default_value
    xl = xlim(ax);
    p  = sort(([xl(1), xl(2)] / d - 1) * 100);

    % Tick step: the finest that keeps the axis to at most 4 labels, which gives
    % -50/0/+50 on a +-50% sweep and -50/0/+50/+100 on a 0.5x-2x sweep.
    step = 200;
    for s = [10 25 50 100 200]
        n = floor(p(2)/s) - ceil(p(1)/s) + 1;
        if n <= 4
            step = s;
            break;
        end
    end
    tp = step * ceil(p(1)/step) : step : step * floor(p(2)/step);
    if ~any(tp == 0)
        tp = sort([tp, 0]);
    end

    % XTick must ascend in DATA coordinates whatever XDir says, so sort there
    % and carry the labels along.
    xt = d * (1 + tp/100);
    [xt, ord] = sort(xt);
    tp = tp(ord);

    labels = arrayfun(@(v) sprintf('%+g%%', v), tp, 'UniformOutput', false);
    labels(tp == 0) = {'0%'};

    set(ax, 'XTick', xt, 'XTickLabel', labels, 'XTickLabelRotation', 0);
    if d < 0
        set(ax, 'XDir', 'reverse');
    end
    xlabel(ax, base_label, ...
        'Interpreter', 'tex', 'FontSize', label_fs);
end
