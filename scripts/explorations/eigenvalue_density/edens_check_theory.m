function C = edens_check_theory(opts)
%EDENS_CHECK_THEORY The radial-density formulas against sampled random matrices, and against each other.
%
%   C = edens_check_theory()
%   C = edens_check_theory('n', 500, 'n_draws', 20, 'out_dir', d, 'save', true)
%
% Five cases, each a matrix X diag(d) with X i.i.d. Gaussian of variance
% 1 / n (cases 1 to 4) or a block variance profile (case 5):
%   1 uniform       every gain 1: the circular law
%   2 broad         gains log-normal, coefficient of variation 1
%   3 zeros         half the gains 0, the rest 1: an atom of mass 1/2 at
%                   the origin and a uniform disk of radius sqrt(1/2)
%   4 two types     Harris et al. (2023), Fig. 2(c): f = 0.8, mu_e =
%                   sigma_e = 1, mu_i = -4, sigma_i = 4, alpha = 0.3
%                   (variances only; the mean is left out)
%   5 blocks        variances that depend on (post, pre) type, with broad
%                   gains and a fifth of them zero: the non-separable case
%                   that edens_radial_theory_blocks is for
% For each: the sampled cumulative radial distribution (n_draws matrices
% of size n) and the theory's, the largest difference between them
% (max_diff, to be read against the sampling error 1 / sqrt(n n_draws)),
% and for cases 1 to 4 the agreement of edens_radial_theory with
% edens_radial_theory_blocks (max_diff_solvers) and, for case 4, with
% Harris's closed form (max_diff_harris, in density relative to uniform).
%
% The figure (eigenvalue_density_theory_check): one panel per case; x is
% the radius over the bulk radius, y the proportion of eigenvalues within
% it; solid: sampled; dashed black: theory; dotted gray: a uniform disk.
%
% See also: edens_radial_theory, edens_radial_theory_blocks, edens_harris_density

arguments
    opts.n       (1,1) double = 500
    opts.n_draws (1,1) double = 20
    opts.seed    (1,1) double = 1
    opts.out_dir (1,:) char = ''
    opts.save    (1,1) logical = true
    opts.visible (1,1) logical = false
end
n = opts.n;
stream = RandStream('mt19937ar', 'Seed', opts.seed);
u = linspace(0, 1, 201);
type_of = [ones(n / 2, 1); 2 * ones(n / 2, 1)];

H = edens_harris_density(0.3, 0.8, 1, -4, 1, 4);
n_e = round(0.8 * n);
d_broad = exp(randn(stream, n, 1) * sqrt(log(2)));       % log-normal, cv 1
d_zero = [zeros(n / 2, 1); ones(n / 2, 1)];
d_mixed = rand(stream, n, 1).^2; d_mixed(rand(stream, n, 1) < 0.2) = 0;
cases = { ...
    'uniform',   ones(n, 1),                                         ones(2) / n; ...
    'broad',     d_broad,                                            ones(2) / n; ...
    'zeros',     d_zero(randperm(stream, n)),                        ones(2) / n; ...
    'two types', [sqrt(H.v_e) * ones(n_e, 1); sqrt(H.v_i) * ones(n - n_e, 1)], ones(2) / n; ...
    'blocks',    d_mixed,                                            [0.10 0.046; 0.046 0.03]};

C = struct('name', cases(:, 1)', 'u', u, 'sampled', [], 'theory', [], 'R', [], ...
    'max_diff', [], 'max_diff_solvers', NaN, 'max_diff_harris', NaN, 'edge_mean_max', []);
for c = 1:size(cases, 1)
    d = cases{c, 2}; s2 = cases{c, 3};
    T = edens_radial_theory_blocks(s2, type_of, d, []);
    T = edens_radial_theory_blocks(s2, type_of, d, u * T.R);
    sd = sqrt(s2(type_of, type_of)) .* d';
    acc = zeros(size(u)); mx = 0;
    for k = 1:opts.n_draws
        a = sort(abs(eig(randn(stream, n) .* sd)));
        acc = acc + sum(a <= u * T.R, 1) / n / opts.n_draws;
        mx = mx + a(end) / opts.n_draws;
    end
    C(c).sampled = acc; C(c).theory = T.n_less; C(c).R = T.R;
    C(c).max_diff = max(abs(acc - T.n_less));
    C(c).edge_mean_max = mx / T.R;
    if c <= 4
        T1 = edens_radial_theory(d.^2, u * T.R);
        C(c).max_diff_solvers = max(abs(T1.n_less - T.n_less));
    end
    if c == 4
        T1 = edens_radial_theory(d.^2, H.r);
        C(c).max_diff_harris = max(abs(T1.rho_rel(2:end - 1) - H.rho_rel(2:end - 1)));
    end
    fprintf('%-10s  sampled vs theory %.4f (sampling error about %.4f)   solvers %.1e   Harris %.1e   largest modulus / R %.3f\n', ...
        C(c).name, C(c).max_diff, 1 / sqrt(n * opts.n_draws), C(c).max_diff_solvers, C(c).max_diff_harris, C(c).edge_mean_max);
end

if ~opts.save; return; end
out_dir = opts.out_dir;
if isempty(out_dir)
    out_dir = fullfile(fileparts(which('setup_paths')), 'figs', 'explorations', 'eigenvalue_density');
end
if ~isfolder(out_dir); mkdir(out_dir); end
st = manuscript_style();
fig = figure('Visible', matlab.lang.OnOffSwitchState(opts.visible), 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 34 8.5]);
cleanup = onCleanup(@() close(fig));
tl = tiledlayout(fig, 1, numel(C), 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:numel(C)
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, u, u.^2, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    plot(ax, u, C(c).sampled, '-', 'Color', [0.000 0.447 0.698], 'LineWidth', 3);
    plot(ax, u, C(c).theory, '--', 'Color', 'k', 'LineWidth', 1.2);
    set(ax, 'FontSize', st.tick_fs, 'LineWidth', st.axis_lw, 'Box', 'off', 'TickDir', 'out', ...
        'XLim', [0 1], 'YLim', [0 1], 'XTick', 0:0.5:1, 'YTick', 0:0.5:1);
    title(ax, C(c).name, 'FontSize', st.title_fs, 'FontWeight', 'normal');
    if c == 1; ylabel(ax, 'fraction within r', 'FontSize', st.label_fs); end
    xlabel(ax, 'r / R', 'FontSize', st.label_fs);
end
tag = 'eigenvalue_density_theory_check';
exportgraphics(fig, fullfile(out_dir, [tag '.png']), 'Resolution', 200);

fid = fopen(fullfile(out_dir, [tag '.md']), 'w');
fprintf(fid, '# %s\n\n', tag);
fprintf(fid, ['The radial-density theory against sampled random matrices (`edens_check_theory`; %d matrices of ' ...
    'size %d per case, Gaussian entries, no mean). Each panel: the fraction of eigenvalues within radius r of the ' ...
    'origin, r in units of the predicted bulk radius R. Thick blue: sampled. Dashed black: theory ' ...
    '(Ahmadian et al. 2015, Eqs. 2.13 to 2.16, solved by `edens_radial_theory_blocks`). Dotted gray: a uniform disk, ' ...
    '(r/R)^2. A curve ABOVE the dotted one is a spectrum concentrated toward the center.\n\n'], opts.n_draws, n);
fprintf(fid, ['Cases: `uniform` every gain 1; `broad` log-normal gains with coefficient of variation 1; `zeros` half ' ...
    'the gains 0 and half 1 (the curve starts at 0.5: the zero-gain neurons'' eigenvalues are exactly at the origin); ' ...
    '`two types` the E and I column variances of Harris et al. (2023) Fig. 2(c) at alpha = 0.3; `blocks` variances by ' ...
    '(post, pre) type [0.10 0.046; 0.046 0.03] with broad gains, a fifth of them zero (the non-separable case).\n\n']);
fprintf(fid, '| case | R | largest difference, sampled vs theory | column solver vs block solver | vs Harris closed form (density / uniform) | mean largest modulus / R |\n|---|---|---|---|---|---|\n');
for c = 1:numel(C)
    fprintf(fid, '| %s | %.4g | %.4f | %.1e | %.1e | %.3f |\n', C(c).name, C(c).R, C(c).max_diff, ...
        C(c).max_diff_solvers, C(c).max_diff_harris, C(c).edge_mean_max);
end
fprintf(fid, '\nSampling error of a proportion here is about %.4f. NaN: not applicable.\n', 1 / sqrt(n * opts.n_draws));
fclose(fid);
end
