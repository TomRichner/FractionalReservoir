%% test_canvas_standard_a.m
% The common canvas of Plan 2 (the manuscript repository's
% docs/Advancing_the_narrative/Plan2_canvas_standard.md, Tom 2026-11-22),
% group A: main Figures 5, 6, 8 and supplementary S6, S7, S8, S10.
%
% Each composite is built from the `fast` run's data in this checkout
% (data/mu7chain_AC_fast; data/IED_model/ied_stim_fast for Figure 8 and
% S10; Figure 8 also needs the clinical panels named by ied_stim_config
% and the model figure figs/IED_model/ied_stim_fast/fig_ied_stim_blocks),
% saved into a temporary folder, and checked:
%   1. the figure is manuscript_style's canvas_w wide and at most
%      canvas_max_h tall;
%   2. every axes, axis label, title, legend, colorbar, text and panel
%      letter is at its role's size (canvas_audit: tick_fs, label_fs,
%      title_fs, letter_fs), except the items named in EXCEPTIONS below,
%      each with its reason; an exception the figure does not use fails
%      too, so the list stays true;
%   3. the exported PNG is within canvas_tol of canvas_w wide (the cropped
%      image fills the width) and no taller than canvas_max_h.
% A figure whose data are not in this checkout is SKIPPED with a line
% saying so; the test fails when every figure is skipped.
%
% Drawing only: nothing is simulated. About two minutes. Assumes
% setup_paths has run and a pool allocation is set (pool_set(0) will do).
% The PNGs are deleted again (Figure 8's hold clinical recordings).
%
% See also: canvas_audit, canvas_type, composite_figure, manuscript_style

st = manuscript_style();
root = fileparts(which('setup_paths'));
run_dir = fullfile(root, 'data', 'mu7chain_AC_fast');
ied_dir = fullfile(root, 'data', 'IED_model', 'ied_stim_fast');
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));

% The exceptions to the type sizes, by figure: {pattern on canvas_audit's
% name, reason}. Every figure not listed has none.
EXCEPTIONS = struct();
EXCEPTIONS.Fig_Main5_Grouped = {'^[EI] \| axes fig_outlier_plane_inset', ...
    'E and I marks inside the response-curve insets (11 pt): at 14 pt they cover the points they label'};

F = {
    'Fig_Main5_Grouped',  @() fig_main5_mechanism('run_dir', run_dir, 'out_dir', fullfile(tmp, 'f5'), 'visible', false), run_dir, fullfile(tmp, 'f5')
    'Fig_Main6_Grouped',  @() fig_main6_amplification('run_dir', run_dir, 'out_dir', fullfile(tmp, 'f6'), 'visible', false), run_dir, fullfile(tmp, 'f6')
    'Fig_Supp6_Grouped',  @() fig_supp6_grouped('run_dir', run_dir, 'out_dir', fullfile(tmp, 's6'), 'visible', false), run_dir, fullfile(tmp, 's6')
    'Fig_Supp7_Grouped',  @() fig_supp7_grouped('run_dir', run_dir, 'out_dir', fullfile(tmp, 's7'), 'visible', false), run_dir, fullfile(tmp, 's7')
    'Fig_Supp8_Grouped',  @() fig_supp8_grouped('run_dir', run_dir, 'out_dir', fullfile(tmp, 's8'), 'visible', false), run_dir, fullfile(tmp, 's8')
    'Fig_Supp10_Grouped', @() fig_supp10_grouped('run_dir', ied_dir, 'ied_run_dir', ied_dir, 'out_dir', fullfile(tmp, 's10'), 'visible', false), ied_dir, fullfile(tmp, 's10')
    'Fig_Main8_Stim',     @() draw_fig8(root, ied_dir, fullfile(tmp, 'f8')), ied_dir, fullfile(tmp, 'f8')
    };

t0 = tic; n_done = 0; skipped = {};
lines = {};
for k = 1:size(F, 1)
    tag = F{k, 1};
    if ~isfolder(F{k, 3})
        skipped{end + 1} = sprintf('%s (no %s)', tag, F{k, 3}); %#ok<SAGROW>
        continue
    end
    before = findall(groot, 'Type', 'figure');
    try
        out = F{k, 2}();
    catch err
        if strcmp(err.identifier, 'test_canvas_standard_a:Skip')
            skipped{end + 1} = sprintf('%s (%s)', tag, err.message); %#ok<SAGROW>
            continue
        end
        rethrow(err);
    end
    fig = out.figs(1);
    ex = cell(0, 2);
    if isfield(EXCEPTIONS, tag); ex = EXCEPTIONS.(tag); end

    % 1. The canvas.
    sz = getappdata(fig, 'composite_size');
    fig.Units = 'inches';
    if isempty(sz) || numel(sz) ~= 2 || sz(1) > 100; sz = fig.Position(3:4); end   % pixels from composite_tools: read the figure
    assert(abs(sz(1) - st.canvas_w) < 1e-6, '%s: %.3f in wide, not canvas_w = %g', tag, sz(1), st.canvas_w);
    assert(sz(2) <= st.canvas_max_h + 1e-6, '%s: %.2f in tall, more than canvas_max_h = %g', tag, sz(2), st.canvas_max_h);

    % 2. The type by role.
    A = canvas_audit(fig, ex);
    bad = A(~A.ok, :);
    assert(isempty(bad), '%s: type off its role''s size: %s', tag, ...
        strjoin(strcat(bad.name, {' ('}, bad.kind, {', '}, string(bad.size), {' pt, want '}, string(bad.want), {')'}), '; '));
    used = ~cellfun(@isempty, A.reason);
    for e = 1:size(ex, 1)
        assert(any(strcmp(A.reason(used), ex{e, 2})), '%s: the exception "%s" is not used any more; take it out', tag, ex{e, 1});
    end

    % 3. The PNG.
    png = fullfile(F{k, 4}, [tag '.png']);
    info = imfinfo(png);
    dpi = info.XResolution;
    if strcmpi(info.ResolutionUnit, 'meter'); dpi = dpi * 0.0254; end
    w_in = info.Width / dpi; h_in = info.Height / dpi;
    assert(w_in >= st.canvas_w - st.canvas_tol && w_in <= st.canvas_w + 0.01, ...
        '%s: the PNG is %.2f in wide; the canvas is %g, the tolerance %g', tag, w_in, st.canvas_w, st.canvas_tol);
    assert(h_in <= st.canvas_max_h + 0.01, '%s: the PNG is %.2f in tall', tag, h_in);

    lines{end + 1} = sprintf('  %-20s %5.2f x %5.2f in, PNG %5.2f x %5.2f in, %3d items of type, %d exception(s)', ...
        tag, sz, w_in, h_in, height(A), nnz(used)); %#ok<SAGROW>
    delete(setdiff(findall(groot, 'Type', 'figure'), before));
    n_done = n_done + 1;
end
assert(n_done > 0, 'Every figure was skipped: %s', strjoin(skipped, '; '));
fprintf('%s\n', lines{:});
note = '';
if ~isempty(skipped); note = sprintf(' (SKIPPED: %s)', strjoin(skipped, '; ')); end
fprintf('PASS test_canvas_standard_a: %d figures on the common canvas%s (%.0f s)\n', n_done, note, toc(t0));

%% ------------------------------------------------------------------------
function out = draw_fig8(root, ied_dir, out_dir)
% Figure 8 through ied_stim_fig_combined, as the ied_stim bundle draws it
% (invisible; it closes its figure), then its saved .fig reopened invisibly
% for the checks.
c = ied_stim_config('fast', 'run_dir', ied_dir, 'fig_root', fullfile(root, 'figs', 'IED_model', 'ied_stim_fast'));
model_fig = fullfile(c.fig_root, 'fig_ied_stim_blocks', 'fig_ied_stim_blocks.fig');
if ~isfolder(fullfile(c.clinical_dir, 'panels')) || ~isfile(model_fig)
    error('test_canvas_standard_a:Skip', 'no clinical panels or no model figure %s', model_fig);
end
o = ied_stim_fig_combined('out_dir', out_dir, 'model_fig', model_fig, 'clinical_dir', c.clinical_dir, ...
    'clinical_files', c.clinical_files, 'run_name', c.run_name, 'sim', c.sim, 'burst_mode_dir', ied_dir, ...
    'tag', 'Fig_Main8_Stim', 'visible', false);
out = struct('figs', openfig(o.files{2}, 'invisible'));
end
