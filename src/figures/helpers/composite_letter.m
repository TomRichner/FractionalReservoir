function h = composite_letter(fig, block, letter, ax, opts)
% COMPOSITE_LETTER A panel letter, "(A)", in the strip above a block of a composite figure.
%
%   COMPOSITE_LETTER(fig, B{1}(1), 'A')            at the block's left edge
%   COMPOSITE_LETTER(fig, B{2}(1), 'C', ax)        at the left edge of ax (its tick labels and axis name included)
%   COMPOSITE_LETTER(..., 'dx', -0.1)              moved, inches
%
% Parenthesized capital, regular weight, the size of the tick labels
% (docs/Figure_style_guide.md). The letter's bottom sits at the block's top
% edge; composite_figure keeps a strip free there.
%
% See also: composite_figure

arguments
    fig
    block (1,1) struct
    letter (1,:) char
    ax = []
    opts.dx (1,1) double = 0
    opts.dy (1,1) double = 0
end
st = manuscript_style();
x = block.rect(1);
if ~isempty(ax)
    drawnow;
    u = ax.Units; ax.Units = 'inches';
    ti = ax.TightInset; x = ax.Position(1) - ti(1);
    ax.Units = u;
end
sz = getappdata(fig, 'composite_size');
if ~isempty(sz) && any(abs(fig.Position(3:4) - sz) > 1e-6); fig.Position(3:4) = sz; end
fp = fig.Position;                                   % inches
x = max(x + opts.dx, 0.02);
y = block.top + 0.02 + opts.dy;
h = annotation(fig, 'textbox', [x / fp(3), y / fp(4), 0.6 / fp(3), 0.28 / fp(4)], 'String', sprintf('(%s)', letter), ...
    'LineStyle', 'none', 'FontSize', st.tick_fs, 'FontWeight', 'normal', 'Margin', 0, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Tag', ['panel_letter_' letter]);
end
