% TEST_IED_STIM_BUNDLE The ied_stim bundle's registry and guards (FD-031 items 4, 5, 14). No simulation.
%
% Fails if a registered figure of the IED stimulation run is removed (the
% model source figure, the standalone rate figure, the supplementary
% state-change figure, Fig 8), if Fig 8 stops reading the run's own model
% figure or the one clinical setting, if fs leaves 500 Hz for new runs, or
% if the run-folder guards stop refusing a mismatched folder.

% ---- registry -------------------------------------------------------------
% Plan 2 (package E) appended Figure 8 in the new order and the appendix
% network-mode figure; the first four entries are as they were.
want = {'fig_ied_stim_blocks', 'fig_ied_stim_blocks_rate', ...
        'fig_ied_stim_blocks_supp_state', 'Fig_Stim_Combined', 'Fig_Main8_Stim', 'fig_ied_network_mode'};
for mode = {'essential', 'fast', 'production', 'stim_blocks_20261117'}
    cfg = ied_stim_config(mode{1});
    F = [cfg.figures{:}]; names = {F.name};
    assert(isequal(names, want), 'registry of %s: %s', mode{1}, strjoin(names, ', '));
    assert(F(3).in_paper && F(4).in_paper, 'supplement and Fig 8 must be in_paper');
    a = struct(F(4).args{:});
    assert(strcmp(a.model_fig, fullfile(cfg.fig_root, 'fig_ied_stim_blocks', 'fig_ied_stim_blocks.fig')), ...
        'Fig 8 must read the run''s own model figure');
    assert(strcmp(a.clinical_dir, cfg.clinical_dir), 'Fig 8 must read cfg.clinical_dir');
    assert(contains(cfg.run_dir, fullfile('data', 'IED_model', cfg.run_name)));
    assert(contains(cfg.fig_root, fullfile('figs', 'IED_model', cfg.run_name)));
    w = cellfun(@(e) e.args{2}, cfg.figures(1:3), 'UniformOutput', false);
    assert(isequal(w, {'main', 'rate', 'supp_state'}));
    a5 = struct(F(5).args{:});
    assert(F(5).in_paper && ~F(6).in_paper, 'new Fig 8 in the paper, the network-mode figure in the appendix');
    assert(isequal(rmfield(a5, {'burst_mode_dir', 'tag'}), a), 'new Fig 8 reads the same model and clinical sources');
    assert(strcmp(a5.burst_mode_dir, cfg.burst_mode_dir) && strcmp(a5.tag, 'Fig_Main8_Stim'));
    ied_network_mode_stage_plan(cfg);                % every mode has a size
end
c = ied_stim_config('fast', 'burst_mode_dir', 'X:/paper_run');
a = struct(c.figures{5}.args{:});
assert(strcmp(a.burst_mode_dir, 'X:/paper_run'), 'burst_mode_dir override must reach the new Fig 8');
assert(endsWith(ied_stim_config('stim_blocks_20261117').burst_mode_dir, 'mu7chain_AC_fastmedium'));
assert(endsWith(ied_stim_config().clinical_dir, 'manuscript_stim/tica_awake_best'), ...
    'clinical default is tica_awake_best (Tom, 2026-11-17)');
c = ied_stim_config('fast', 'clinical_dir', 'X:/elsewhere');
a = struct(c.figures{4}.args{:});
assert(strcmp(a.clinical_dir, 'X:/elsewhere'), 'clinical_dir override must reach Fig 8');

% ---- the experiment ------------------------------------------------------
p = ied_stim_config('production').sim;
assert(isequal(p.seeds, [42 73 101 137 211 307 409 503 601 709 811 907]));
assert(p.n_no == 4 && p.n_stim == 4 && p.block_len == 300 && p.amplitude == 0.2);
for mode = {'essential', 'fast', 'production'}
    assert(ied_stim_config(mode{1}).sim.fs == 500, 'fs 500 Hz for new runs (FD-031 item 14)');
end
assert(ied_stim_config('stim_blocks_20261117').sim.fs == 400, 'the 2026-11-17 run is 400 Hz');

% ---- guards ----------------------------------------------------------------
d = tempname; mkdir(d); cleanup = onCleanup(@() rmdir(d, 's'));
settings = struct('seeds', 42, 'amplitude', 0.2, 'n_no', 1, 'n_stim', 1, ...
    'block_len', 60, 'settle', 10, 'fs', 400); %#ok<NASGU>
save(fullfile(d, 'run_settings.mat'), 'settings');
threw = false;
try
    ied_stim_blocks('data_dir', d, 'seeds', 42, 'n_no', 1, 'n_stim', 1, 'block_len', 60, 'fs', 500);
catch err
    threw = strcmp(err.identifier, 'ied_stim_blocks:SettingsDiffer');
end
assert(threw, 'ied_stim_blocks must refuse a folder written at another fs');

opts = struct('seeds', [42 73], 'amplitude', 0.2, 'n_no', 1, 'n_stim', 1, ...
    'block_len', 60, 'settle', 10, 'fs', 400); rows = table(); %#ok<NASGU>
save(fullfile(d, 'summary.mat'), 'rows', 'opts');
cfg = ied_stim_config('essential', 'run_dir', d, 'fig_root', fullfile(d, 'figs'));
threw = false;
try
    ied_stim_figures(cfg);
catch err
    threw = strcmp(err.identifier, 'ied_stim_figures:WrongRun');
end
assert(threw, 'ied_stim_figures must refuse a run directory of other settings');

% ---- retired scripts stay retired ------------------------------------------
assert(isempty(which('make_fig_stim_combined_v2')) && isempty(which('make_fig_stim_combined')));
r = fileparts(which('ied_stim_blocks'));
assert(endsWith(r, fullfile('scripts', 'paper')), 'ied_stim_blocks lives in scripts/paper');
disp('ied_stim bundle: registry, fs, clinical setting and guards checked.');
