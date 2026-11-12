% TEST_MU7CHAIN_BUNDLE The chained-depression bundle's config, presets and name maps (FD-007).
%
% Quick: no stage and no figure runs. Checks that
%   * both figure sets of mu7chain_config resolve, with their own fig_root,
%     conditions and grouped registry, and one shared run_dir;
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
%   * run_all_paper_analyses' run-mode pre-flight accepts the config's mode.
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
assert(strcmp(cfg.run_mode, 'fast') && strcmp(cfg.figure_set, '4cond'));
assert(isequal(cfg.conditions, cfg4.conditions) && strcmp(cfg.fig_root, cfg4.fig_root));
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
% The same registry as the manuscript bundle's, entry for entry.
ref_cfg = mu7revisedMedium_config();
assert(isequal(cellfun(@(f) f.name, cfg3.figures, 'UniformOutput', false), ...
               cellfun(@(f) f.name, ref_cfg.figures, 'UniformOutput', false)));

% FD-016: a fixed sample seed, valid for run_all_paper_analyses.
s = cfg.param_space_seed;
assert(isnumeric(s) && isscalar(s) && s >= 0 && s == floor(s) && s <= 2^32 - 1);
assert(isequal(cfg3.param_space_seed, cfg4.param_space_seed, s));

% run_all_paper_analyses' pre-flight is ismember(cfg.run_mode, run_mode_names())
% (it is inline in the entry point, which is not run here).
assert(ismember(cfg.run_mode, run_mode_names()));

%% Every preset the config names exists and validates
stage_presets = {cfg.preset_name, cfg.mc_preset, cfg.local_lyapunov_preset, ...
                 cfg.single_neuron_preset};
assert(strcmp(cfg.preset_name, main_name));
for k = 1:numel(stage_presets)
    srnn_param_preset(stage_presets{k});   % validates its conditions, or throws
end

[d, cls, cnd] = srnn_param_preset(main_name);
assert(strcmp(cls, 'SRNNCellTypePairs'));
assert(isequal(name_of(cnd), four_names), 'conditions: %s', strjoin(name_of(cnd), ', '));

% The network is the manuscript's, value for value; the shared names are the
% manuscript preset's structs exactly.
[d_ref, ~, cnd_ref] = srnn_param_preset(ref_name);
assert(isequaln(d, d_ref), 'the mu7chain network differs from mu7revisedMedium');
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

disp('PASS: test_mu7chain_bundle -- two figure sets, presets, conditions, name maps, chain and single-neuron runs.');
