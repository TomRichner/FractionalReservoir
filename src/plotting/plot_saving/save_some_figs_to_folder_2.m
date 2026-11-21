function [] = save_some_figs_to_folder_2(save_folder, save_name, fig_vec, fig_type)

if not(exist(save_folder,'dir'))
    mkdir(save_folder)
end
% There was a warning('folder already exists') here. Saving into an existing
% folder is the NORMAL case -- every figure that writes more than one file hits
% it, and every rerun into a stable fig_root hits it for all of them -- so it
% fired dozens of times per run and said nothing. Its cost was not the noise
% itself but what the noise buried: a real warning in that output had no chance
% of being noticed.

if isempty(fig_vec)
    figHandles = findobj('Type', 'figure');   % the missing ';' printed the handle array once per save
    for i_f = 1:length(figHandles)
        fig_vec(i_f) = figHandles(i_f).Number;
    end
end

if isempty(fig_type)
    fig_type = {'fig','png'};
end
    
h = get(0,'children');
h = flipud(h);

% Initialize PDF report if requested
if any(strcmpi(fig_type, 'pdf'))
    pdf_file = fullfile(save_folder, [save_name '_report.pdf']);
    if exist(pdf_file, 'file')
        delete(pdf_file);
    end
end

% Each figure is saved inside a try/catch so that ONE failed export cannot
% abort the caller. This matters because run_all_analyses saves figures at the
% end of each stage: an export that throws there discards a completed
% multi-minute sweep and everything downstream of it, to lose a picture that can
% be regenerated from the saved .mat at any time.
%
% The failure this guards against is real and not understood. Roughly one run in
% four, a LIVE figure reaches a state where every rasterizing path --
% exportgraphics at any resolution, print -dpng, saveas png -- throws
% MATLAB:class:InvalidHandle, while the figure handle and all of its descendants
% report valid, saveas to .fig still works, a freshly built figure exports fine,
% and rebuilding the same figure from the same saved data exports fine. So it is
% neither a stale handle in this function nor a broken renderer nor anything
% about the data. Until it is diagnosed, warn and carry on.
% SAVING MUST NOT RAISE OR UNHIDE A FIGURE. exportgraphics needs a graphics
% object, and fig_vec holds figure NUMBERS, so this used to bridge the gap with
% figure(i). That is not a cast: figure(f) brings f above every other window,
% takes keyboard focus, and sets Visible = 'on'. Three consequences, all real:
%
%   - a batch run stole focus once per figure per format, so the machine could
%     not be used while it ran;
%   - it silently DEFEATED the 'visible', false that make_all_paper_figures
%     passes to every figure -- they were built hidden and then shown here;
%   - if i was not an existing figure number, figure(i) CREATED an empty one and
%     exported that, writing a blank file rather than reporting the mistake.
%
% findobj resolves the number to the handle without touching visibility or
% stacking order. An unmatched number now errors into the catch below and is
% reported, which is the behaviour the blank-file case should always have had.
failed = {};
for i=fig_vec
    try
        fh = findobj(groot, 'Type', 'figure', 'Number', i);
        if isempty(fh)
            error('save_some_figs_to_folder_2:NoSuchFigure', ...
                'No open figure numbered %d.', i);
        end
        fh = fh(1);
        set(fh,'PaperPositionMode','auto')
        size_before = fh.Position;      % the size the figure was laid out for
        oversize = larger_than_screen(fh);
        if any(strcmpi(fig_type,'fig'))
            % A figure saved while invisible reopens invisible (openfig keeps the
            % stored Visible), and every pipeline figure is drawn invisible. The
            % CreateFcn runs when the .fig is opened and turns it on (2026-09-14).
            set(fh, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');
            saveas(fh, fullfile(save_folder, [save_name '_f_' num2str(i)]), 'fig');
        end
        if any(strcmpi(fig_type,'png'))
            png_file = fullfile(save_folder, [save_name '_figure_' num2str(i) '.png']);
            if ~oversize
                exportgraphics(fh, png_file, 'Resolution', png_resolution(fh))
                drawnow;                % a resize can arrive after the export returns
                oversize = was_resized(fh, size_before);
            end
            if oversize
                print_whole(fh, size_before, png_file, '-dpng');
            end
        end
        if any(strcmpi(fig_type,'svg'))
            % Was set(gcf, ...), which targeted whatever figure happened to be
            % current rather than fh. It worked only by accident: the png branch
            % above had just raised fh with figure(i), making it current. With
            % the raise gone that accident goes too, so the target is now named.
            set(fh, 'Renderer', 'painters');
            svg_file = fullfile(save_folder, [save_name '_figure_' num2str(i) '.svg']);
            if ~oversize
                exportgraphics(fh, svg_file, 'BackgroundColor', 'none', 'ContentType', 'vector');
                drawnow;
                oversize = was_resized(fh, size_before);
            end
            if oversize
                print_whole(fh, size_before, svg_file, '-dsvg');
            end
        end
        if any(strcmpi(fig_type, 'pdf'))
            % Append to PDF report
            exportgraphics(fh, pdf_file, 'Append', true, 'ContentType', 'vector');
        end
    catch save_err
        failed{end+1} = sprintf('%d (%s)', i, save_err.identifier); %#ok<AGROW>
        warning('save_some_figs_to_folder_2:FigureSaveFailed', ...
            'Could not save figure %d as ''%s'': %s\nContinuing with the rest.', ...
            i, save_name, save_err.message);
    end
end
if ~isempty(failed)
    warning('save_some_figs_to_folder_2:SomeFiguresNotSaved', ...
        '%d of %d figures did not save for ''%s'': %s', ...
        numel(failed), numel(fig_vec), save_name, strjoin(failed, ', '));
end
end

function tf = larger_than_screen(fh)
% True for a figure that does not fit the screen with a tenth to spare.
%
% WHY (2026-11-21, Plan 2 package D). MATLAB shrinks a figure window that
% is larger than the screen to the screen (on the 1710 x 1073 workstation a
% 12 x 15.8 in figure came back 10.97 in tall), also an invisible one, and
% it does so at a moment of its choosing: usually during an export. What
% was measured on a figure 8 x 16 in with marks at its top, middle and
% bottom:
%   exportgraphics  NEVER gave the whole figure. The PNG had either the
%                   bottom 10.97 in only, or all three marks squeezed onto
%                   a canvas 10.97 in tall; the SVG was always 10.97 in.
%                   Whether the window was resized during the call did not
%                   tell the two apart.
%   print           with the paper set to the figure's size gave the whole
%                   figure EXACTLY WHEN the window had its full size before
%                   the call and still had it after; a print during which
%                   the window was shrunk had the top cut off. Setting
%                   Position back takes on one attempt and not on the next.
% On a figure with many small axes laid out in inches the first
% exportgraphics of a session can come out whole, which is why the damage
% looked random: Fig_Ostruct_Plane_rotation and _blockmean (5 to 7 rows)
% lost their column titles and first row on 2026-11-21 while the 3-row
% local group, which fits the screen, was intact. It is also the "resized
% on every other export" that fig_grouped_main works around.
%
% So a figure that does not fit goes to print_whole from the start, and a
% figure that fits takes the path it always took (exportgraphics, cropped
% to its content). The tenth to spare covers the window's title bar and
% the task bar; a figure in that margin is printed uncropped, which is
% harmless.
tf = false;
try
    px  = getpixelposition(fh);
    scr = get(groot, 'ScreenSize');
    tf  = px(3) > 0.9 * scr(3) || px(4) > 0.9 * scr(4);
catch
    % no screen to compare with: keep the usual path
end
end

function tf = was_resized(fh, size_before)
% True when the figure no longer has the size it was laid out for.
tf = any(abs(fh.Position(3:4) - size_before(3:4)) > 1e-6 * max(size_before(3:4)));
end

function print_whole(fh, size_before, file, device)
% Export a figure larger than the screen at its own size (see
% larger_than_screen): print with the paper set to the figure's size,
% repeated until one print has run with the window at its full size before
% AND after. The canvas is not cropped to the content. The size is put
% back at the end as well as it can be: the caller's figure should keep
% the size it was laid out for, but MATLAB may shrink the window again.
units = fh.Units;
fh.Position = size_before;
fh.Units = 'inches';
sz = fh.Position(3:4);
fh.Units = units;
set(fh, 'PaperUnits', 'inches', 'PaperPositionMode', 'manual', 'PaperSize', sz, ...
    'PaperPosition', [0 0 sz], 'InvertHardcopy', 'off');
whole = false;
for attempt = 1:10
    fh.Position = size_before;
    drawnow;
    if was_resized(fh, size_before); continue; end      % the size did not take: set it again
    if strcmp(device, '-dpng')
        print(fh, file, device, sprintf('-r%d', png_resolution(fh)));
    else
        print(fh, file, device, '-vector');
    end
    if ~was_resized(fh, size_before); whole = true; break; end
end
if ~whole
    warning('save_some_figs_to_folder_2:NotWhole', ...
        ['%s: the figure is larger than the screen and no print ran with its window at full size ' ...
         '(10 attempts); the file is probably cut off at the top.'], file);
end
for attempt = 1:3                                       % leave the figure at its own size
    if ~was_resized(fh, size_before); break; end
    fh.Position = size_before;
    drawnow;
end
end

function res = png_resolution(fh)
% PNG resolution for fh: 600 dpi, lowered so the rendered canvas stays within
% MAX_CANVAS pixels on its longer side.
%
% WHY (2026-11-12, FD-006). Above 4096 pixels the raster export draws the
% canvas in tiles that do not line up: text is cut along the seam (a tick
% label '30' came out as '3', '0' as an arch, '0.2' with its top half missing)
% and lines on the far side of the seam are shifted against their axes (a step
% at t = 10 s drawn at 10.3 s). Measured on the single-neuron figure, 1000 x 850
% px: clean at 200 dpi, the seam at 4096 px of the uncropped canvas at 430 and
% 600 dpi, and the same exported from the reopened .fig, so it is the export,
% not the figure. rendererinfo reports a 16384 texture limit, so the limit
% cannot be read from it. The .svg and .fig are unaffected. A figure up to
% 6.67 in (640 px at 96 ppi) still gets 600 dpi; a larger one gets
% proportionally fewer dpi and so about 4000 px, still ample for print.
MAX_CANVAS = 4000;
res = 600;
try
    px  = getpixelposition(fh);
    ppi = get(groot, 'ScreenPixelsPerInch');
    res = min(600, floor(MAX_CANVAS * ppi / max(px(3:4))));
catch
    % no pixel size to be had: keep 600
end
end
