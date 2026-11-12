% TEST_FIG_STD_1TS_VS_MTS The 1TS-vs-chain depression figure against an independent integration (FD-014).
%
% About 30 s with no pool (most of it the two exports). Checks that
%   * the figure builds in both figure sets of mu7chain_config (invisible,
%     closed afterwards), skips no_adaptation, and writes its three files and
%     its table;
%   * the traces it returns AND the lines it draws agree with an independent
%     ode45 integration of (C1)-(C2) (the chain), of the single pool, and of
%     the product of two independent pools, all written out below with their
%     parameters read from the preset's synapse_config, not from a route model;
%   * the first- and last-cycle peak outputs printed in the table are those of
%     that independent integration;
%   * at the onset of the first rate step from rest the 1TS and chain factors
%     coincide: |F_chain - F_1TS| <= 2e-4 over the first 0.25 s (4.4e-5
%     measured), while F_1TS has fallen to 0.79;
%   * the robustness grid's (0.1 Hz, 0.5) entry is the main sinusoid's ratio.
%
% Asserts; silent until it throws, then prints one PASS line.
% See also: fig_std_1TS_vs_MTS, route_pool, test_std_chain_model

preset = 'celltype_pairs_sfaEI_Sc0p2sig0p1_tauSpread0p25_stdSpread0p25_noStim_noise0p025_chainStd_4cond_mu7chain';
[~, ~, cnd] = srnn_param_preset(preset);
cnames = cellfun(@(c) c.name, cnd, 'UniformOutput', false);
par = @(nm) cnd{strcmp(cnames, nm)}.synapse_config.E.E.std;
p1 = par('sfa1_std1');  p2 = par('sfa3_stdc2');  pp = par('sfa3_std2');

% ---- The equations, written out (independent of src/model/routes) --------
% state y: 1TS b; chain b1, b2; product b_a, b_b
rhs = @(y, r) [ ...
    (1 - y(1)) / p1.tau_rec - y(1) * r / p1.tau_rel; ...
    y(3) * (1 - y(2)) / p2.tau_rec(1) - y(2) * r / p2.tau_rel; ...
    (1 - y(3)) / p2.tau_rec(2) - y(3) * (1 - y(2)) / (p2.n2 * p2.tau_rec(1)); ...
    (1 - y(4)) / pp.tau_rec(1) - y(4) * r / pp.tau_rel(1); ...
    (1 - y(5)) / pp.tau_rec(2) - y(5) * r / pp.tau_rel(2)];
ode = odeset('RelTol', 1e-11, 'AbsTol', 1e-13);
indep = @(Y, r) struct( ...
    'sfa1_std1',  struct('D', Y(:, 1)', 'b2', [], 'theta', r .* Y(:, 1)'), ...
    'sfa3_stdc2', struct('D', Y(:, 2)', 'b2', Y(:, 3)', 'theta', r .* Y(:, 2)'), ...
    'sfa3_std2',  struct('D', (Y(:, 4) .* Y(:, 5))', 'b2', [], 'theta', r .* (Y(:, 4) .* Y(:, 5))'));

sets = {mu7chain_config('fast', '3cond'), mu7chain_config('fast', '4cond')};
for s = 1:2
    k_cfg = sets{s};
    tmp = tempname; mkdir(tmp);
    cleanup = onCleanup(@() rmdir(tmp, 's'));
    names = cellfun(@(f) f.name, k_cfg.figures, 'UniformOutput', false);
    e = k_cfg.figures{strcmp(names, 'fig_std_1TS_vs_MTS')};
    assert(e.in_paper, 'fig_std_1TS_vs_MTS must be in_paper in the %s set', k_cfg.figure_set);
    % The call make_all_paper_figures makes: standard fields, the bundle's
    % conditions, then the entry's args. The robustness grid is cut to two
    % frequencies at one amplitude (it is most of the figure's run time, and
    % only its 0.1 Hz / 0.5 cell is checked here); everything drawn is the
    % registry's call.
    args = [{'out_dir', tmp, 'save', true, 'visible', false, 'conditions', k_cfg.conditions}, ...
        e.args, {'grid_freq', [0.1 0.2], 'grid_rmax', 0.5}];
    o = fig_std_1TS_vs_MTS(args{:});
    fig_closer = onCleanup(@() close(o.figs(isvalid(o.figs))));
    want = k_cfg.conditions(~strcmp(k_cfg.conditions, 'no_adaptation'));
    assert(isequal(o.conditions, want), '%s set draws {%s}', k_cfg.figure_set, strjoin(o.conditions, ', '));
    assert(isequal(o.skipped, {'no_adaptation'}));
    assert(strcmp(o.reference, 'sfa1_std1'));
    assert(isequal(sort(o.files), sort({'Fig_STD_1TS_vs_MTS.png', 'Fig_STD_1TS_vs_MTS.svg', ...
        'Fig_STD_1TS_vs_MTS.fig'})), 'files: %s', strjoin(o.files, ', '));
    assert(isfile(o.table));

    %% Independent integration: sinusoid (one ode45 call, smooth rate)
    sn = o.sine;
    rfun = @(t) (sn.rmax / 2) * (1 - cos(2 * pi * sn.f * t));
    [~, Ys] = ode45(@(t, y) rhs(y, rfun(t)), sn.t, ones(5, 1), ode);
    Is = indep(Ys, sn.r);
    %% Rate steps (one ode45 call per constant-rate segment)
    stp = o.steps;
    Yt = ones(numel(stp.t), 5);
    % stp.r(i + 1) is the rate over (t(i), t(i + 1)], so the rate changes at
    % sample i wherever stp.r(i + 1) ~= stp.r(i).
    edges = unique([1, find(diff(stp.r) ~= 0), numel(stp.t)]);
    for q = 1:numel(edges) - 1
        i0 = edges(q); i1 = edges(q + 1);
        rq = stp.r(i0 + 1);
        [~, Yq] = ode45(@(~, y) rhs(y, rq), stp.t(i0:i1), Yt(i0, :)', ode);
        if i1 - i0 == 1; Yq = Yq([1 end], :); end
        Yt(i0 + 1:i1, :) = Yq(2:end, :);
    end
    It = indep(Yt, stp.r);

    tol = 1e-7;   % RK4 at 1 kHz against ode45 at 1e-11
    for k = 1:numel(o.conditions)
        nm = o.conditions{k};
        for f = {'D', 'theta', 'b2'}
            pairs = {sn.traces(k).(f{1}), Is.(nm).(f{1}), 'sine'; ...
                     stp.traces(k).(f{1}), It.(nm).(f{1}), 'steps'};
            for q = 1:2
                a = pairs{q, 1}; b = pairs{q, 2};
                assert(isequal(size(a), size(b)), '%s %s %s: size', nm, pairs{q, 3}, f{1});
                if isempty(a); continue; end              % no reserve: none on either side
                err = max(abs(a - b));
                assert(err <= tol, '%s %s %s: max diff %.2e', nm, pairs{q, 3}, f{1}, err);
            end
        end
    end

    %% The lines the figure DRAWS are those traces
    st = manuscript_style();
    lines_ = findall(o.figs, 'Type', 'line');
    n_checked = 0;
    for h = lines_'
        x = h.XData; y = h.YData;
        if numel(x) < 100; continue; end                 % markers, steady-state segments
        cname = '';
        for k = 1:numel(o.conditions)
            if isequal(h.Color, st.condition_color(o.conditions{k})); cname = o.conditions{k}; end
        end
        if isempty(cname); continue; end                 % the rate, legend dummies
        if abs(x(end) - sn.t(end)) < 1e-9
            I = Is.(cname); tt = sn.t;
        elseif abs(x(end) - stp.t(end)) < 1e-9
            I = It.(cname); tt = stp.t;
        else
            continue                                     % steady-state curves (x = rate)
        end
        idx = round(x / (tt(2) - tt(1))) + 1;
        cand = {I.D, I.theta, I.b2};
        err = inf;
        for c = 1:numel(cand)
            if ~isempty(cand{c}); err = min(err, max(abs(y - cand{c}(idx)))); end
        end
        assert(err <= tol, 'a drawn %s line differs from the independent integration by %.2e', cname, err);
        n_checked = n_checked + 1;
    end
    % sine: D + theta per condition + b2 of the chain; steps the same
    assert(n_checked == 2 * (2 * numel(o.conditions) + 1), 'checked %d drawn traces', n_checked);

    %% The table's first- and last-cycle peaks are the independent ones
    txt = fileread(o.table);
    tok = regexp(txt, '\| `(\w+)` \| ([\d.]+) \| [+-][\d.]+ \| ([\d.]+) \| [+-]', 'tokens');
    assert(numel(tok) == numel(o.conditions), 'sinusoid table: %d rows', numel(tok));
    n_per = round(1 / (sn.f * (sn.t(2) - sn.t(1))));
    for q = 1:numel(tok)
        th = Is.(tok{q}{1}).theta;
        first = max(th(1:n_per + 1)); last = max(th(end - n_per:end));
        assert(abs(str2double(tok{q}{2}) - first) <= 5.1e-5 && ...
               abs(str2double(tok{q}{3}) - last) <= 5.1e-5, ...
            '%s: table %s / %s, independent %.5f / %.5f', tok{q}{1}, tok{q}{2}, tok{q}{3}, first, last);
    end

    %% At the first onset from rest, 1TS and the chain coincide
    i0 = stp.onset_idx(1);
    w = i0:i0 + round(0.25 / (stp.t(2) - stp.t(1)));
    d = max(abs(It.sfa3_stdc2.D(w) - It.sfa1_std1.D(w)));
    assert(d <= 2e-4, 'first onset: 1TS and chain factors differ by %.2e in 0.25 s', d);
    assert(It.sfa1_std1.D(w(end)) < 0.95, 'the window must include real depression');

    %% The grid's (0.1 Hz, 0.5) cell is the main sinusoid
    gi = find(o.grid.freq == sn.f); ga = find(o.grid.rmax == sn.rmax);
    kc = find(strcmp(o.conditions, 'sfa3_stdc2'));
    assert(abs(o.grid.ratio_first(gi, ga, kc) - o.sine.first(kc) / o.sine.first(1)) < 1e-12);
    assert(abs(o.grid.ratio_last(gi, ga, kc) - o.sine.last(kc) / o.sine.last(1)) < 1e-12);

    clear fig_closer cleanup
end

disp('PASS: test_fig_std_1TS_vs_MTS -- both sets build; traces, drawn lines and table match an independent integration; 1TS and chain coincide at onset.');
