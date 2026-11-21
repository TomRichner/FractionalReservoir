function Mz = soied_measures(R, thr)
% SOIED_MEASURES Fractions beyond the boundary, crossings and excursions against discharges, subnet attribution, full-Jacobian network mode.
%
%   Mz = soied_measures(R)          % R from soied_run (one network)
%   Mz = soied_measures(R, thr)     % coherence threshold (default 0.3)
%
% THE COHERENT MODE at a sample is the rightmost eigenvalue of M (Im >= 0,
% among the n_top rightmost) whose unit eigenvector has at least thr of
% its norm^2 in the span of the ten population indicators (coh10; random
% vector 10/350 = 0.029). rc is its real part.
%
% Per block (1 baseline, 2 stimulation), on the uniform coarse grid:
%   frac_l1    fraction of time Re lambda_1(M) > 1 (any mode, bulk included)
%   frac_coh   fraction of time rc > 1
%   frac_mc    the same for the most coherent of the n_top rightmost modes
%   R_eff, R10 mean bulk radius predicted from the gains (2-type, 10-pop)
%   dbar       mean gains (post x pre, E/I), rate (E, I)
% Full Jacobian (1-s grid): frac_net_pos (rightmost network mode Re > 0),
%   net_re (median, 1/s), net_im (median |Im|), net_coh (median coh10),
%   frac_netc_pos (rightmost network mode with coh10 >= 0.2, Re > 0).
%
% Against discharges (merged coarse + dense grid):
%   STRICT CROSSINGS: upward passages of rc through 1. ev_above: the
%   discharge has rc > 1 somewhere in [t - 0.6, t + 0.2]; ev_cross: an
%   upward passage in [t - 0.6, t]; lag_cross: its time minus t (the last
%   one); n_cross, n_cross_free: passages, and those with no discharge in
%   [passage - 0.2, passage + 0.6].
%   EXCURSIONS: theta = 95th percentile of rc over baseline coarse samples
%   at least 1 s from every discharge. ev_exc: rc > theta in
%   [t - 0.6, t + 0.05]; lag_exc: first such time minus t; n_exc:
%   excursions (upward passages of theta on the coarse grid, passages less
%   than 0.3 s apart merged); exc_rate_free: those with no discharge in
%   [onset - 0.2, onset + 0.6] s, per minute.
%   ATTRIBUTION: at the sample of max rc in [t - 0.3, t + 0.05] the
%   coherent mode's subnet shares; dom_in_event = its largest-share subnet
%   is one of the subnets the detector marked; share_event = its share on
%   those subnets.
% event_full (soied_event_full): median profiles of the full network
%   mode around discharges and control times, and the fraction with Re > 0
%   somewhere in [-0.3, 0] s.
%
% Event time t is the detector's: the dendritic peak of the discharge.
%
% See also: soied_run, soied_main_figure, soied_write_sidecar

if nargin < 2, thr = 0.3; end
Mz = struct('seed', R.seed, 'thr', thr, 'check', R.check);
for b = 1:2
    B = R.block(b); Q = B.spec; t = Q.t(:); c = Q.is_coarse(:);
    rc = coherent_re(Q, thr);
    [cm, k] = max(Q.coh10, [], 2);
    re_mc = real(Q.lam(sub2ind(size(Q.lam), (1:numel(t))', k)));
    w = struct();
    w.n_events = numel(B.events.time);
    w.events_per_min = w.n_events / (R.settings.block_len / 60);
    w.frac_l1 = mean(real(Q.lam(c, 1)) > 1);
    w.frac_coh = mean(rc(c) > 1);
    w.frac_mc = mean(re_mc(c) > 1);
    w.coh_l1 = median(Q.coh10(c, 1));
    w.re_l1 = mean(real(Q.lam(c, 1)));
    w.re_coh = mean(rc(c), 'omitnan');
    w.re_coh_p = prctile(rc(c), [5 50 95]);
    w.R_eff = mean(Q.R_eff(c)); w.R10 = mean(Q.R10(c));
    w.K10_1 = mean(real(Q.K10_eig(c, 1)));
    w.dbar = squeeze(mean(Q.dbar(c, :, :), 1));
    w.rate = mean(Q.rate(c, :), 1);
    F = B.full;
    w.frac_net_pos = mean(real(F.net) > 0);
    w.net_re = median(real(F.net)); w.net_re_mean = mean(real(F.net));
    w.net_im = median(imag(F.net)); w.net_coh = median(F.net_coh10);
    w.frac_netc_pos = mean(real(F.netc) > 0);
    w.kin_edge = mean(F.kin_edge);

    % strict crossings of 1
    te = B.events.time(:);
    up = find(rc(1:end - 1) <= 1 & rc(2:end) > 1 & diff(t) < 0.06) + 1;
    tu = t(up);
    w.n_cross = numel(tu);
    w.n_cross_free = sum(arrayfun(@(x) ~any(te >= x - 0.2 & te <= x + 0.6), tu));
    w.ev_above = arrayfun(@(x) any(rc(t >= x - 0.6 & t <= x + 0.2) > 1), te);
    w.lag_cross = nan(numel(te), 1);
    for e = 1:numel(te)
        q = tu(tu >= te(e) - 0.6 & tu <= te(e));
        if ~isempty(q), w.lag_cross(e) = q(end) - te(e); end
    end
    w.ev_cross = isfinite(w.lag_cross);
    % excursions above the baseline quiet 95th percentile
    if b == 1
        quiet = c & arrayfun(@(x) all(abs(te - x) >= 1), t);
        theta = prctile(rc(quiet), 95);
        Mz.theta = theta; Mz.rc_quiet_p = prctile(rc(quiet), [5 50 95]);
    end
    w.ev_exc = false(numel(te), 1); w.lag_exc = nan(numel(te), 1);
    for e = 1:numel(te)
        s = t >= te(e) - 0.6 & t <= te(e) + 0.05 & rc > Mz.theta;
        if any(s), w.ev_exc(e) = true; w.lag_exc(e) = t(find(s, 1)) - te(e); end
    end
    tcg = t(c); rcg = rc(c);
    upc = tcg(find(rcg(1:end - 1) <= Mz.theta & rcg(2:end) > Mz.theta) + 1);
    if numel(upc) > 1                      % one excursion: passages less than 0.3 s apart
        upc = upc([true; diff(upc(:)) > 0.3]);
    end
    free = arrayfun(@(x) ~any(te >= x - 0.2 & te <= x + 0.6), upc);
    w.n_exc = numel(upc); w.n_exc_free = sum(free);
    w.exc_rate_free = sum(free) / (R.settings.block_len / 60);
    % attribution
    w.dom_in_event = nan(numel(te), 1); w.share_event = nan(numel(te), 1);
    w.dom = nan(numel(te), 1); w.n_groups = nan(numel(te), 1);
    for e = 1:numel(te)
        s = find(t >= te(e) - 0.3 & t <= te(e) + 0.05);
        if isempty(s), continue; end
        [~, i] = max(rc(s)); i = s(i);
        kk = coherent_index(Q, thr, i);
        if isnan(kk), continue; end
        sh = squeeze(double(Q.share(i, kk, :)));
        g = B.events.groups{e};
        [~, w.dom(e)] = max(sh);
        w.dom_in_event(e) = ismember(w.dom(e), g);
        w.share_event(e) = sum(sh(g));
        w.n_groups(e) = numel(g);
    end
    % event-triggered median of rc on the dense grid
    lags = R.settings.dense_win(1):R.settings.dense_dt:R.settings.dense_win(2);
    P = nan(numel(te), numel(lags)); P1 = P;
    for e = 1:numel(te)
        [dd, ii] = min(abs(t - (te(e) + lags)), [], 1);
        ok = dd < R.settings.dense_dt / 2;
        P(e, ok) = rc(ii(ok)); P1(e, ok) = real(Q.lam(ii(ok), 1));
    end
    w.trig_lags = lags; w.trig_rc = median(P, 1, 'omitnan'); w.trig_l1 = median(P1, 1, 'omitnan');
    Mz.block(b) = w;
end
if isfield(R, 'event_full')
    E = R.event_full;
    pre = E.lags >= -0.3 & E.lags <= 0;
    Mz.ef = struct('lags', E.lags, ...
        'ev_net_med', median(real(E.event.net), 1), 'ctrl_net_med', median(real(E.ctrl.net), 1), ...
        'ev_Mcoh_med', median(real(E.event.M_coh), 1, 'omitnan'), ...
        'ctrl_Mcoh_med', median(real(E.ctrl.M_coh), 1, 'omitnan'), ...
        'ev_pos_pre', mean(any(real(E.event.net(:, pre)) > 0, 2)), ...
        'ctrl_pos_pre', mean(any(real(E.ctrl.net(:, pre)) > 0, 2)), ...
        'ev_max_pre', max(real(E.event.net(:, pre)), [], 2)', ...
        'ctrl_max_pre', max(real(E.ctrl.net(:, pre)), [], 2)', ...
        'ev_net_coh_pre', median(E.event.net_coh10(:, pre), 'all'), ...
        'repro_ok', E.repro_ok, 'n', numel(E.t_event));
    [~, lp] = max(Mz.ef.ev_net_med); Mz.ef.lag_peak = E.lags(lp);
    fp = find(Mz.ef.ev_net_med > 0 & E.lags >= -0.4, 1);
    if isempty(fp), Mz.ef.lag_first_pos = NaN; else, Mz.ef.lag_first_pos = E.lags(fp); end
end
Mz.ospec_full_probe = R.ospec_full_probe;
end

function rc = coherent_re(Q, thr)
re = real(double(Q.lam)); re(Q.coh10 < thr | isnan(re)) = -Inf;
rc = max(re, [], 2); rc(~isfinite(rc)) = NaN;
end

function kk = coherent_index(Q, thr, i)
re = real(double(Q.lam(i, :))); re(Q.coh10(i, :) < thr) = -Inf;
[v, kk] = max(re); if ~isfinite(v), kk = NaN; end
end
