function [T, settings, truth] = make_synthetic_mc_gain_sweep(conditions)
% MAKE_SYNTHETIC_MC_GAIN_SWEEP A fake block C result of the planned shape, for the figures' tests.
%
%   [T, settings, truth] = MAKE_SYNTHETIC_MC_GAIN_SWEEP()
%   [T, settings, truth] = MAKE_SYNTHETIC_MC_GAIN_SWEEP(conditions)
%
% NOT DATA. A table with the columns Plan C section 2 fixes for
% run_mc_gain_sweep (one row per gain, bias, seed and condition), filled
% with a pattern chosen so that the statistics have a known answer:
%   * gains logspace(-1, 1, 9) (gain 1 is a level), biases -0.2:0.2:0.8
%     (0.2 is a level), seeds 1:3;
%   * lambda_1 = 1.5 log10(gain / g_c(bias)), g_c = 10^(0.5 (bias - 0.2)),
%     zero along a ridge through (gain 1, bias 0.2); adapting conditions
%     sit 0.4 1/s lower (more stable);
%   * no_adaptation: a NARROW ridge of memory capacity along lambda_1 = 0,
%     peak 11 at (1, 0.2);
%   * sfa1_std1: a broad bump, peak 8 at (1, 0.2);
%   * sfa3_stdc2: a broader plateau, peak 9 at (gain 1.78, bias 0.4);
%   * sfa3_std2 (when asked for): the chain's shape at 0.9 of its height;
%   * a saturated corner (gain >= 5.6, bias <= 0: 80 % of neurons pinned,
%     memory capacity 0.5) and a silent corner (gain <= 0.18, bias >= 0.6:
%     90 % silent, memory capacity 0.3), in every condition;
%   * one failed job (no_adaptation, gain 0.316, bias 0.4, seed 2);
%   * +-4 % seed-to-seed noise (fixed generator), so medians recover the
%     shapes above.
% truth holds the peak locations and the expected order of the shares.
%
% See also: test_mc_gain_bias_figures, read_mc_gain_sweep

if nargin < 1 || isempty(conditions)
    conditions = {'no_adaptation', 'sfa1_std1', 'sfa3_stdc2', 'sfa3_std2'};
end
gains  = logspace(-1, 1, 9);
biases = -0.2:0.2:0.8;
seeds  = 1:3;
delay_s = 0.3:0.3:15;

s = RandStream('mt19937ar', 'Seed', 7);
rows = {};
for c = 1:numel(conditions)
    cn = conditions{c};
    for g = gains
        for b = biases
            for sd = seeds
                gc  = 10^(0.5 * (b - 0.2));
                lam = 1.5 * log10(g / gc);
                lg = log10(g);
                switch cn
                    case 'no_adaptation'
                        mc = 0.5 + 10.5 * exp(-(lam / 0.12)^2);
                    case 'sfa1_std1'
                        lam = lam - 0.4;
                        mc = 0.5 + 7.5 * exp(-lg^2 / (2 * 0.5^2)) * exp(-(b - 0.2)^2 / (2 * 0.45^2));
                    case 'sfa3_stdc2'
                        lam = lam - 0.4;
                        mc = 0.5 + 8.5 * exp(-(lg - 0.25)^2 / (2 * 0.6^2)) * exp(-(b - 0.4)^2 / (2 * 0.6^2));
                    otherwise   % sfa3_std2 and anything else: the chain at 0.9
                        lam = lam - 0.4;
                        mc = 0.9 * (0.5 + 8.5 * exp(-(lg - 0.25)^2 / (2 * 0.6^2)) * exp(-(b - 0.4)^2 / (2 * 0.6^2)));
                end
                pinned = 0.02; silent = 0.02;
                if g >= 5.6 && b <= 0
                    pinned = 0.8; mc = 0.5;
                elseif g <= 0.18 && b >= 0.6
                    silent = 0.9; mc = 0.3;
                end
                mc  = mc * (1 + 0.04 * (2 * rand(s) - 1));
                lam = lam + 0.02 * randn(s);
                failed = strcmp(cn, 'no_adaptation') && abs(g - gains(3)) < 1e-12 && ...
                         abs(b - 0.4) < 1e-12 && sd == 2;
                tau = 0.3 / log(1 + 1 / mc);              % sum of exp(-d/tau) over delays ~ mc
                R2 = exp(-delay_s / tau);
                H = tau * log(2);
                if failed; mc = NaN; H = NaN; R2(:) = NaN; lam = NaN; end
                rows(end+1, :) = {g, b, sd, cn, mc, H, R2, lam, lam - 0.05, ...
                    0.2 + 0.1 * lg, 0.25 + 0.1 * lg, pinned, silent, failed}; %#ok<AGROW>
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', {'gain', 'bias', 'seed', 'condition', 'mc_total', ...
    'horizon', 'R2_delay', 'lle_input', 'lle_noinput', 'rate_E', 'rate_I', ...
    'frac_pinned', 'frac_silent', 'failed'});
settings = struct('conditions', {conditions}, 'delay_s', delay_s, ...
    'preset_name', '', 'preset_gain', 1, 'preset_bias', 0.2, 'synthetic', true);

truth = struct();
truth.peak = struct('no_adaptation', [1 0.2], 'sfa1_std1', [1 0.2], ...
                    'sfa3_stdc2', [gains(6) 0.4], 'sfa3_std2', [gains(6) 0.4]);
% Share of the map at >= half of the best of ANY condition: broad beats narrow.
truth.share_any_order = {'sfa3_stdc2', 'sfa1_std1', 'no_adaptation'};
truth.n_failed = 1;
end
