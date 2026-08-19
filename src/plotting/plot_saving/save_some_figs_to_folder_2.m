function [] = save_some_figs_to_folder_2(save_folder, save_name, fig_vec, fig_type)

if not(exist(save_folder,'dir'))
    mkdir(save_folder)
else
    warning('folder already exists')
end

if isempty(fig_vec)
    figHandles = findobj('Type', 'figure')
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
% PNG export resolution, and the raster size above which it must be reduced --
% see the tiling comment in the png branch below.
png_res    = 600;
png_max_px = 3400;   % bisected safe width -- see the png branch below

failed = {};
for i=fig_vec
    try
        set(i,'PaperPositionMode','auto')
        if any(strcmpi(fig_type,'fig'))
            saveas(i, fullfile(save_folder, [save_name '_f_' num2str(i)]), 'fig');
        end
        if any(strcmpi(fig_type,'png'))
            png_file = fullfile(save_folder, [save_name '_figure_' num2str(i) '.png']);
            exportgraphics(figure(i), png_file, 'Resolution', png_res)
            % Wide rasters LOSE GLYPHS. Tick labels came out as "500" -> "50",
            % "1000" -> "10", "700" -> ")0", "+100%" -> ")0%" -- a couple of
            % characters silently eaten out of one label, which reads as a typo
            % rather than as a broken export.
            %
            % What it is NOT, all checked: not a layout collision (the axes'
            % XTickLabel reads back correct and the vector .svg of the same
            % figure is clean); not the intermittent InvalidHandle failure below
            % (this is deterministic and reproduces byte-for-byte); not
            % exportgraphics specifically (`print -dpng -r600` eats the same
            % label); not the renderer (forcing 'painters' changes nothing).
            % Across widths of 2075 / 3520 / 3759 / 4001 / 5999 / 7639 / 8125 px
            % on the same figures, the damage always lands near the end of a
            % rasteriser tile -- so which label is hit moves with the width,
            % which is why resizing a figure appears to "fix" it while only
            % moving the problem to a different panel.
            %
            % The one thing that reliably avoids it is a small enough raster.
            % Bisected on this machine: 3520 px clean, 3759 px damaged, hence the
            % cap below. Re-export at the highest resolution that respects it.
            % Measured rather than predicted, since the pixel size depends on the
            % screen's points-per-inch and so varies by machine and display
            % scaling. Figures small enough to be safe at 600 dpi -- most of
            % them, since this needs a sheet wider than ~6 in -- take this branch
            % never and are bit-identical to before.
            try
                png_info = imfinfo(png_file);
                png_big  = max(png_info(1).Width, png_info(1).Height);
                if png_big > png_max_px
                    safe_res = max(150, floor(png_res * png_max_px / png_big));
                    exportgraphics(figure(i), png_file, 'Resolution', safe_res)
                    warning('save_some_figs_to_folder_2:PngResolutionCapped', ...
                        ['Figure %d (''%s'') was %d px at %d dpi, above the %d px ' ...
                         'tiling threshold where exportgraphics drops glyphs; ' ...
                         're-exported at %d dpi. The .svg is full quality.'], ...
                        i, save_name, png_big, png_res, png_max_px, safe_res);
                end
            catch res_err
                warning('save_some_figs_to_folder_2:PngSizeCheckFailed', ...
                    'Could not check the PNG size for figure %d (%s); left as exported.', ...
                    i, res_err.message);
            end
        end
        if any(strcmpi(fig_type,'svg'))
            set(gcf, 'Renderer', 'painters');
            exportgraphics(figure(i), fullfile(save_folder, [save_name '_figure_' num2str(i) '.svg']), 'BackgroundColor', 'none', 'ContentType', 'vector');
        end
        if any(strcmpi(fig_type, 'pdf'))
            % Append to PDF report
            exportgraphics(figure(i), pdf_file, 'Append', true, 'ContentType', 'vector');
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