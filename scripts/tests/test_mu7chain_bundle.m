% TEST_MU7CHAIN_BUNDLE The chained-depression bundle's config, presets and name maps (FD-007).
%
% Quick: no stage and no figure runs. Checks that
%   * both figure sets of mu7chain_config resolve, with their own fig_root,
%     conditions and grouped registry, and one shared run_dir;
%   * the conditions simulated and the sets drawn per mode (FD-030): three
%     and {'3cond'} from 'fast' upward, all four and {'4cond'} at
%     'essential', and the three are conditions of every stage's preset;
%   * every preset the config names exists and validates, the network is the
%     manuscript's (mu7revisedMedium) value for value, and the twins differ
%     from it only where they say;
%   * the four conditions are as specified, in order, and the three shared
%     names are the mu7revisedMedium structs exactly;
%   * full_adaptation_condition picks the chain, for the four conditions and
%     for the paper's three, and the plural returns {sfa3_stdc2, sfa3_std2};
%   * every condition name of the bundle has a title, a colour and a short
%     label (two figures index those maps unguarded);
%   * the chained-MTS network builds and runs for one second at reduced size,
%     with finite output and identical routes;
%   * the single-neuron preset builds and runs in all four conditions;
%   * the single-neuron figure's registry entry names its own conditions in
%     each set, and with the entry's args after the bundle's (the order
%     make_all_paper_figures uses) it draws three columns in the 3cond set
%     and four in the 4cond set;
%   * run_all_paper_analyses' run-mode pre-flight accepts the config's mode;
%   * the bundle's defaults of FD-046, in every mode: the four spectral
%     stages (outlier_spectra, frequency_response, eigenvalue_density,
%     outlier_structures) after lyapunov_spectrum_n80, their plans accept
%     the bundle, their figure entries are appended (none in_paper; the
%     density figures also for EEx1p5; fig_ostruct_plane once per group;
%     the conditions the entry's stage simulates), outlier_spectra gives
%     'reference' and 'EEx1p5' the dense pass that three of its figures
%     need, and Figure 3 is
%     the switching version: fig_main3_grouped (in_paper) is the four-row
%     switching composite, the baseline composite is
%     fig_main3_grouped_baseline (not in_paper), K = 100 and 40 states from
%     1 s; the A + C bundle inherits all of it; the regression bundle gets
%     none of it.
%
% Asserts; silent until it throws, then prints one PASS line.
% See also: mu7chain_config, mu7chain_run, test_preset_conditions

main_name = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
ref_name  = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_noStim_noise0p025_dualStd_3cond_mu7revisedMedium';
four_names  = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};
three_names = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
name_of = @(cnd) cellfun(@(c) c.name, cnd, 'UniformOutput', false);

%% The config and its two figure sets
cfg  = mu7chain_config();
cfg3 = mu7chain_config('fast', '3cond');
cfg4 = mu7chain_config('fast', '4cond');
% No set named: the mode's own, the three-condition set from 'fast' upward
% (FD-030; until then '4cond').
assert(strcmp(cfg.run_mode, 'fast') && strcmp(cfg.figure_set, '3cond'));
assert(isequal(cfg.conditions, cfg3.conditions) && strcmp(cfg.fig_root, cfg3.fig_root));
assert(isequal(cfg3.conditions, three_names), '3cond set: %s', strjoin(cfg3.conditions, ', '));
assert(isequal(cfg4.conditions, four_names),  '4cond set: %s', strjoin(cfg4.conditions, ', '));
assert(strcmp(cfg3.run_dir, 'data/mu7chain_fast') && strcmp(cfg4.run_dir, cfg3.run_dir));
assert(strcmp(cfg3.fig_root, 'figs/mu7chain_fast_3cond'));
assert(strcmp(cfg4.fig_root, 'figs/mu7chain_fast_4cond'));
cfg_med = mu7chain_config('medium', '3cond');
assert(strcmp(cfg_med.run_dir, 'data/mu7chain_medium') && ...
       strcmp(cfg_med.fig_root, 'figs/mu7chain_medium_3cond'));
threw = false;
try
    mu7chain_config('not_a_mode');
catch err
    threw = strcmp(err.identifier, 'mu7chain_config:BadRunMode');
end
assert(threw, 'an unknown run mode must be refused');
% These modes run every stage; only 'essential' selects (FD-023; its own
% checks are in test_essential_mode). Every stage is written out, followed
% by the opt-in lyapunov_spectrum_n80 that every mode runs (FD-032) and the
% four spectral stages (FD-046, with_spectral_supplement).
spectral = {'outlier_spectra', 'frequency_response', 'eigenvalue_density', 'outlier_structures'};
n_spec = numel(spectral);
for k = {cfg3, cfg4, cfg_med}
    assert(isequal(k{1}.stages, [paper_stage_names(), {'lyapunov_spectrum_n80'}, spectral]) && ...
        isempty(k{1}.sweeps) && isempty(k{1}.eig_heatmap_examples));
end
% FD-030: the conditions SIMULATED and the sets DRAWN, by mode. Every mode
% from 'fast' upward simulates the paper's three and draws the 3cond set,
% whichever set the config was built for; 'essential' simulates all four
% (run_conditions empty = every condition of the preset, exactly as before)
% and draws the 4cond set.
for m = setdiff(run_mode_names(), {'essential'}, 'stable')
    for fs = {'', '3cond', '4cond'}
        k = mu7chain_config(m{1}, fs{1});
        assert(isequal(k.run_conditions, three_names), '%s: run_conditions %s', m{1}, ...
            strjoin(k.run_conditions, ', '));
        assert(isequal(k.figure_sets, {'3cond'}), '%s: figure_sets', m{1});
    end
end
cfg_ess = mu7chain_config('essential');
assert(isempty(cfg_ess.run_conditions) && isequal(cfg_ess.figure_sets, {'4cond'}));
assert(strcmp(cfg_ess.figure_set, '4cond') && isequal(cfg_ess.conditions, four_names));
assert(isempty(bundle_run_conditions(cfg_ess)));
assert(isequal(bundle_run_conditions(cfg3), three_names));
% The simulated conditions are conditions of every preset a stage runs on
% (run_all_paper_analyses checks the same before any compute), in the
% preset's order.
for p = {cfg3.preset_name, cfg3.mc_preset, cfg3.local_lyapunov_preset}
    [~, ~, sel] = preset_with_conditions(p{1}, cfg3.run_conditions);
    assert(isequal(name_of(sel), three_names), 'preset %s', p{1});
    [~, ~, all_c] = preset_with_conditions(p{1}, {});
    assert(isequal(name_of(all_c), four_names), 'preset %s, no selection', p{1});
end
threw = '';
try
    preset_with_conditions(main_name, {'no_adaptation', 'sfa9_std9'});
catch err
    threw = err.identifier;
end
assert(strcmp(threw, 'preset_with_conditions:NoSuchCondition'), 'unknown condition: %s', threw);
threw = '';
try
    bundle_run_conditions(struct('run_conditions', {{1, 2}}));
catch err
    threw = err.identifier;
end
assert(strcmp(threw, 'bundle_run_conditions:BadConfig'));

% The grouped registry is rebuilt per set: each set's grouped entries read
% their sources from that set's own fig_root (inventory R6).
for c = {cfg3, cfg4}
    k = c{1};
    names = cellfun(@(f) f.name, k.figures, 'UniformOutput', false);
    assert(numel(unique(names)) == numel(names), 'duplicate registry entry');
    grouped = k.figures(startsWith(names, 'fig_main') & endsWith(names, '_grouped'));
    assert(numel(grouped) == 8);
    for g = 1:numel(grouped)
        a = grouped{g}.args;
        src = a{find(strcmp(a(1:2:end), 'source_fig_root')) * 2};
        assert(strcmp(src, k.fig_root), '%s reads %s, not its own %s', ...
            grouped{g}.name, src, k.fig_root);
    end
end
% The same registry as the manuscript bundle's, entry for entry, plus the
% chain's own comparison figure (FD-014), in both sets and in the paper,
% the ranked Lyapunov spectrum of the 80-neuron network (FD-032), and the
% real-axis histogram of the Jacobian's eigenvalues against a uniform disk
% (2026-11-19).
ref_cfg = mu7revisedMedium_config();
% (FD-046) Also the spectral supplement's entries, the switching Figure 3's
% stand-alone versions and the baseline composite kept beside it.
chain_only = {'fig_std_1TS_vs_MTS', 'fig_lyapunov_ranked_spectrum', 'fig_eig_real_histogram', ...
    'fig_eig_heatmap_switch', 'fig_eig_real_histogram_switch', 'fig_main3_grouped_baseline'};
chain_only = [chain_only, spectral_entry_names()];
for c = {cfg3, cfg4}
    nm = cellfun(@(f) f.name, c{1}.figures, 'UniformOutput', false);
    assert(isequal(nm(~ismember(nm, chain_only)), ...
                   cellfun(@(f) f.name, ref_cfg.figures, 'UniformOutput', false)));
    e = c{1}.figures{strcmp(nm, 'fig_std_1TS_vs_MTS')};
    assert(e.in_paper && isequal(e.args, {'preset_name', main_name}));
end

% FD-016: a fixed sample seed, valid for run_all_paper_analyses.
s = cfg.param_space_seed;
assert(isnumeric(s) && isscalar(s) && s >= 0 && s == floor(s) && s <= 2^32 - 1);
assert(isequal(cfg3.param_space_seed, cfg4.param_space_seed, s));

% run_all_paper_analyses' pre-flight is ismember(cfg.run_mode, run_mode_names())
% (it is inline in the entry point, which is not run here).
assert(ismember(cfg.run_mode, run_mode_names()));

%% FD-046: the spectral supplement and the switching Figure 3 are the default
three_spec = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2'};
for m = run_mode_names()
    for fs = {'3cond', '4cond'}
        for ac = [false true]
            if ac; k = mu7chain_AC_config(m{1}, fs{1}); else; k = mu7chain_config(m{1}, fs{1}); end
            tag = sprintf('%s %s%s', m{1}, fs{1}, repmat(' (A + C)', 1, ac));
            % The stages: after lyapunov_spectrum_n80, before block C's; the
            % plans accept the bundle (the reading stages refuse a bundle
            % without outlier_spectra).
            i80 = find(strcmp(k.stages, 'lyapunov_spectrum_n80'));
            assert(isscalar(i80) && isequal(k.stages(i80 + 1:i80 + n_spec), spectral), '%s: the stages', tag);
            assert(isequal(paper_stage_plan(k).stages(end - (n_spec - 1) - 2 * ac:end - 2 * ac), spectral), '%s: run order', tag);
            pl = outlier_spectra_stage_plan(k);
            frequency_response_stage_plan(k); eigenvalue_density_stage_plan(k);
            ps = outlier_structures_stage_plan(k);
            assert(isequal(pl.conditions, three_spec) && isequal(ps.conditions, three_spec), '%s: the stages'' conditions', tag);
            % fig_outlier_plane, _slow_variables and _dwell draw 'reference'
            % and 'EEx1p5' and need their dense pass and full-Jacobian
            % series: the stage's default full_networks, every run of them.
            dense_nets = {'reference', 'EEx1p5'};
            assert(isequal(pl.full_networks, dense_nets) && all(ismember(dense_nets, pl.networks)) && ...
                all([pl.runs(ismember({pl.runs.network}, dense_nets)).full]) && ...
                ~any([pl.runs(~ismember({pl.runs.network}, dense_nets)).full]), '%s: the dense pass of reference and EEx1p5', tag);
            % The entries, none in the paper, each with the conditions it can draw.
            names = cellfun(@(f) f.name, k.figures, 'UniformOutput', false);
            want = spectral_entry_names();
            [present, at] = ismember(want, names);
            assert(all(present) && all(diff(at) == 1), '%s: the spectral entries, in order', tag);
            for e = k.figures(at)
                f = e{1};
                assert(~f.in_paper && isa(f.fn, 'function_handle'), '%s: %s', tag, f.name);
                a = f.args;
                if strcmp(f.name, 'fig_eigdens_theory_check')
                    assert(isempty(a), '%s: %s args', tag, f.name);
                else
                    c = a{find(strcmp(a(1:2:end), 'conditions')) * 2};
                    assert(isequal(c, intersect(k.conditions, three_spec, 'stable')), '%s: %s conditions', tag, f.name);
                    assert(endsWith(f.name, '_EEx1p5') == any(strcmp(a(1:2:end), 'network')), '%s: %s network', tag, f.name);
                    assert(startsWith(f.name, 'fig_ostruct_plane_') == any(strcmp(a(1:2:end), 'group')), '%s: %s group', tag, f.name);
                    if startsWith(f.name, 'fig_ostruct_plane_')
                        assert(strcmp(a{find(strcmp(a(1:2:end), 'group')) * 2}, extractAfter(f.name, 'fig_ostruct_plane_')) && ...
                            strcmp(func2str(f.fn), 'fig_ostruct_plane'), '%s: %s', tag, f.name);
                    elseif endsWith(f.name, '_EEx1p5')
                        assert(strcmp(a{find(strcmp(a(1:2:end), 'network')) * 2}, 'EEx1p5') && ...
                            strcmp(func2str(f.fn), extractBefore(f.name, '_EEx1p5')), '%s: %s', tag, f.name);
                    else
                        assert(strcmp(func2str(f.fn), f.name), '%s: %s', tag, f.name);
                    end
                end
            end
            % Figure 3: the switching composite under the paper's name; the
            % baseline beside it; K = 100 and 40 states at 1, 2, ..., 40 s.
            m3 = k.figures{strcmp(names, 'fig_main3_grouped')};
            b3 = k.figures{strcmp(names, 'fig_main3_grouped_baseline')};
            assert(m3.in_paper && ~b3.in_paper && isequal(m3.fn, b3.fn) && ...
                find(strcmp(names, 'fig_main3_grouped_baseline')) == find(strcmp(names, 'fig_main3_grouped')) + 1, '%s: Figure 3 entries', tag);
            assert(isequal(m3.args, [b3.args, {'spectra_source', 'local_lyapunov', 'ritz_overlay', true, 'ritz_inset', true, ...
                'show_lambda_row', false, 'real_histogram_row', true}]), '%s: the four-row switching composite', tag);
            assert(~any(ismember({'fig_main3_grouped_switch', 'fig_main3_grouped_switch_hist'}, names)) && ...
                all(ismember({'fig_eig_heatmap_switch', 'fig_eig_real_histogram_switch'}, names)), '%s: switching versions', tag);
            assert(k.local_lyapunov_K == 100 && k.local_lyapunov_spectra_n_samples == 40 && ...
                k.local_lyapunov_spectra_start_s == 1 && k.local_lyapunov_accumulation_start_s == 0 && ...
                k.local_lyapunov_warmup_s == 10 && k.local_lyapunov_simulation_start_s == -15, '%s: the timing', tag);
        end
    end
end
% The regression bundle gets none of it.
rg = regression_fast_config();
rgn = cellfun(@(f) f.name, rg.figures, 'UniformOutput', false);
assert(~any(ismember(spectral, paper_stage_plan(rg).stages)) && ~any(ismember(spectral_entry_names(), rgn)) && ...
    ~isfield(rg, 'local_lyapunov_spectra_n_samples') && ~ismember('fig_main3_grouped_baseline', rgn), ...
    'the regression bundle has no spectral stage, entry or switching Figure 3');
% Applying either helper again is refused.
assert(throws_id(@() with_spectral_supplement(cfg3), 'with_spectral_supplement:AlreadyAdded'));
assert(throws_id(@() with_switch_spectra(cfg3), 'with_switch_spectra:AlreadyAdded'));

%% Every preset the config names exists and validates
stage_presets = {cfg.preset_name, cfg.mc_preset, cfg.local_lyapunov_preset, ...
                 cfg.single_neuron_preset, cfg.lyapunov_spectrum_n80_preset};
assert(strcmp(cfg.preset_name, main_name));
for k = 1:numel(stage_presets)
    srnn_param_preset(stage_presets{k});   % validates its conditions, or throws
end

% The 80-neuron twin (FD-032): the network with only its size changed, and
% the same four conditions exactly.
[d80, cls80, cnd80] = srnn_param_preset(cfg.lyapunov_spectrum_n80_preset);
[d500, ~, cnd500] = srnn_param_preset(main_name);
size_fields = {'n', 'indegree', 'F_ref_n', 'F_ref_indegree'};
assert(strcmp(cls80, 'SRNNCellTypePairs') && isequal(cnd80, cnd500));
assert(isequal(cellfun(@(f) d80.(f), size_fields), [80 16 80 16]));
assert(isequal(rmfield(d80, size_fields), rmfield(d500, size_fields)), ...
    'the 80-neuron preset differs from the mu7chain network in more than its size');

[d, cls, cnd] = srnn_param_preset(main_name);
assert(strcmp(cls, 'SRNNCellTypePairs'));
assert(isequal(name_of(cnd), four_names), 'conditions: %s', strjoin(name_of(cnd), ', '));

% The network is the manuscript's, value for value, plus the per-neuron
% depression spread its name states (stdSpread0p25, FD-015); the shared
% names are the manuscript preset's structs exactly.
[d_ref, ~, cnd_ref] = srnn_param_preset(ref_name);
assert(isequal(d.std_spread, [0.25 0.25]), 'the mu7chain network must set std_spread = [0.25 0.25]');
assert(~isfield(d_ref, 'std_spread'), 'the manuscript preset has no depression spread');
assert(isequaln(rmfield(d, 'std_spread'), d_ref), ...
    'the mu7chain network differs from mu7revisedMedium in more than std_spread');
for nm = {'no_adaptation', 'sfa1_std1', 'sfa3_std2'}
    assert(isequaln(cnd{strcmp(name_of(cnd), nm{1})}, cnd_ref{strcmp(name_of(cnd_ref), nm{1})}), ...
        '%s differs from mu7revisedMedium''s', nm{1});
end

% The chained MTS: all four routes, the specified chain, three SFA timescales.
chain = cnd{3};
want = struct('model', 'chain', 'tau_rec', [2 5], 'tau_rel', 0.25, 'n2', 2);
for pre = {'E', 'I'}
    for post = {'E', 'I'}
        assert(isequal(chain.synapse_config.(pre{1}).(post{1}).std, want), ...
            'chain route %s->%s', pre{1}, post{1});
        assert(isequal(cnd{4}.synapse_config.(pre{1}).(post{1}).std, ...
            struct('tau_rec', [2 4], 'tau_rel', [0.25 0.5])));
        assert(isequal(cnd{2}.synapse_config.(pre{1}).(post{1}).std, ...
            struct('tau_rec', 2, 'tau_rel', 0.25)));
    end
end
assert(isequal(chain.tau_a, cnd{4}.tau_a) && numel(chain.tau_a{1}) == 3);

% The twins differ only where they say.
[d_mc, ~, cnd_mc] = srnn_param_preset(cfg.mc_preset);
d_x = d; d_x.sigma_u_noise = 0;
assert(isequaln(d_mc, d_x) && isequaln(cnd_mc, cnd));
[d_st, ~, cnd_st] = srnn_param_preset(cfg.local_lyapunov_preset);
d_x = d; d_x.input_config = d_st.input_config;
assert(isequaln(d_st, d_x) && isequaln(cnd_st, cnd));
assert(d_st.input_config.n_steps == 8 && ~any(d_st.input_config.no_stim_pattern));

%% full_adaptation_condition picks the chain
assert(strcmp(full_adaptation_condition(cnd), 'sfa3_stdc2'));
assert(strcmp(full_adaptation_condition(cnd(1:3)), 'sfa3_stdc2'), ...
    'the paper''s three conditions must give the chain, not sfa1_std1');
assert(isequal(full_adaptation_conditions(cnd), {'sfa3_stdc2', 'sfa3_std2'}));
% By rule, not by order: the product listed first still loses the tie.
assert(strcmp(full_adaptation_condition(cnd([1 2 4 3])), 'sfa3_stdc2'));
assert(isequal(full_adaptation_conditions(cnd([4 3 2 1])), {'sfa3_stdc2', 'sfa3_std2'}));
% Unchanged on a product-only preset.
assert(strcmp(full_adaptation_condition(cnd_ref), 'sfa3_std2'));
assert(isequal(full_adaptation_conditions(cnd_ref), {'sfa3_std2'}));

%% Every condition of the bundle has a title, a colour and a short label
st = manuscript_style();
[~, ~, cnd_one] = srnn_param_preset(cfg.single_neuron_preset);
all_names = unique([name_of(cnd), name_of(cnd_one)]);
for k = 1:numel(all_names)
    nm = all_names{k};
    assert(isKey(st.condition_title, nm), 'no title for %s', nm);
    assert(isKey(st.condition_color, nm), 'no colour for %s', nm);
    assert(isKey(st.condition_short, nm), 'no short label for %s', nm);
end
assert(isequal(st.condition_color('sfa3_std2'), [0 0.7 0]));
assert(isequal(st.condition_color('sfa3_stdc2'), [0.000 0.447 0.698]));
assert(strcmp(st.condition_short('sfa3_stdc2'), 'SFA+STD'));
assert(strcmp(st.condition_short('sfa3_std2'), 'SFA+prodSTD'));
assert(strcmp(st.condition_title('sfa3_stdc2'), 'Multiple-Timescale Adaptation'));
% The four drawn colours are distinct, and so are the titles.
cols = cell2mat(cellfun(@(n) st.condition_color(n), four_names, 'UniformOutput', false)');
assert(size(unique(cols, 'rows'), 1) == 4);
assert(numel(unique(cellfun(@(n) st.condition_title(n), four_names, 'UniformOutput', false))) == 4);

%% The chained MTS builds and runs at reduced size
m = build_from_preset(main_name, 'sfa3_stdc2', 'n', 40, 'indegree', 8, ...
    'lya_method', 'none', 'T_range', [0 1]);
evalc('m.run();');
assert(m.n == 40 && all(isfinite(m.S_out(:))), 'non-finite output on the chain network');
assert(SRNNCellTypePairs.routes_identical(m.get_params()));
p = m.get_params();
assert(all(p.std_chain(:)), 'every route of sfa3_stdc2 must be a chain');

%% The single-neuron preset: four conditions, each builds and runs
assert(isequal(name_of(cnd_one), {'no_adaptation', 'sfa1_std0', 'sfa0_std1', 'sfa0_stdc2'}));
assert(isequal(cnd_one{4}.synapse_config.E.E.std, want));
assert(isequal(cnd_one{3}.synapse_config.E.E.std, struct('tau_rec', 2, 'tau_rel', 0.25)));
assert(isequal(cnd_one{2}.tau_a{1}, cnd{2}.tau_a{1}(1)));
assert(strcmp(full_adaptation_condition(cnd_one), 'sfa0_stdc2'));
for q = 1:4
    m1 = build_from_preset(cfg.single_neuron_preset, cnd_one{q}.name, ...
        'lya_method', 'none', 'T_range', [0 1]);
    evalc('m1.run();');
    assert(m1.n == 1 && all(m1.W(:) == 0) && all(isfinite(m1.S_out(:))), ...
        'single neuron, %s', cnd_one{q}.name);
end

%% The single-neuron figure: its own condition list per set, 3 and 4 columns
% The registry entry names the single-neuron conditions in its own args,
% which follow the bundle-wide 'conditions' in make_all_paper_figures' call
% and so override it. Checked by drawing the figure with exactly that
% argument order, on a synthetic single_neuron stage holding the preset's
% four conditions (the figure only replots; no simulation is needed).
sn_want = struct('c3', {{'no_adaptation', 'sfa1_std0', 'sfa0_std1'}}, ...
                 'c4', {{'no_adaptation', 'sfa1_std0', 'sfa0_std1', 'sfa0_stdc2'}});
tmp_run = tempname; mkdir(fullfile(tmp_run, 'single_neuron'));
cleanup_run = onCleanup(@() rmdir(tmp_run, 's'));
tt = linspace(0, 30, 61);
res = struct('name', name_of(cnd_one), 'title', '', 't', tt, 'u', tt, 'x', tt, 'r', tt, ...
    'syn', tt, 'sfa', tt, 'std', tt, 'std_model', {'none', 'none', 'product', 'chain'}, 'n', 1);
results = res; settings = struct('display_window', [0 30]);
save(fullfile(tmp_run, 'single_neuron', 'single_neuron_data.mat'), 'results', 'settings');
for c = {{cfg3, 'c3', 3}, {cfg4, 'c4', 4}}
    k = c{1}{1}; want = sn_want.(c{1}{2}); n_col = c{1}{3};
    names = cellfun(@(f) f.name, k.figures, 'UniformOutput', false);
    e = k.figures{strcmp(names, 'fig_single_neuron_revised')};
    assert(isequal(e.args, {'conditions', want}), 'single-neuron args, %s set', k.figure_set);
    assert(isequal(k.single_neuron_conditions, want));
    args = [{'run_dir', tmp_run, 'save', false, 'visible', false, 'conditions', k.conditions}, e.args];
    o = fig_single_neuron_revised(args{:});
    n_ax = numel(findall(o.figs, 'Type', 'axes'));
    close(o.figs);
    assert(n_ax == 6 * n_col, '%s set: %d axes, want %d columns', k.figure_set, n_ax, n_col);
end

disp('PASS: test_mu7chain_bundle -- two figure sets, conditions per mode, presets, conditions, name maps, chain and single-neuron runs, the spectral and Figure 3 defaults (FD-046).');

function names = spectral_entry_names()
% The entries with_spectral_supplement appends, in order.
names = {'fig_outlier_plane', 'fig_outlier_dose', 'fig_outlier_slow_variables', 'fig_outlier_dwell', ...
    'fig_freqresp_network', 'fig_freqresp_plane', 'fig_freqresp_phase', 'fig_freqresp_single_neuron'};
for d = {'radial', 'gains', 'measures', 'plane'}
    names = [names, {['fig_eigdens_' d{1}], ['fig_eigdens_' d{1} '_EEx1p5']}]; %#ok<AGROW>
end
names = [names, {'fig_eigdens_theory_check', 'fig_ostruct_summary', 'fig_ostruct_plane_rotation', ...
    'fig_ostruct_plane_local', 'fig_ostruct_plane_blockmean', 'fig_ostruct_rotation', 'fig_ostruct_h1'}];
end

function ok = throws_id(fn, id)
ok = false;
try
    fn();
catch ME
    ok = strcmp(ME.identifier, id);
end
end
