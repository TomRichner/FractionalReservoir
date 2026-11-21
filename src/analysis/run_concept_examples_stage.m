function detail = run_concept_examples_stage(cfg, run_dir, plan)
% RUN_CONCEPT_EXAMPLES_STAGE A bundle's concept_examples stage: the example networks of Figure 1A and the gallery of S3, into <run_dir>/concept_examples.
%
%   detail = RUN_CONCEPT_EXAMPLES_STAGE(cfg, run_dir, plan)
%
%   cfg      the bundle's config (unused; the plan holds what is used)
%   run_dir  the bundle's run directory
%   plan     from concept_examples_stage_plan(cfg)
%   detail   the file written and what is in it, on one line
%
% For every network of plan.examples and plan.gallery (cx_network builds
% it; nothing in the model is changed):
%   - the weight statistics: realized block means and variances (zeros
%     included), the fraction of nonzero entries and the mean and SD of the
%     nonzero entries per (post, pre) block, the structural mean matrix and
%     its eigenvalues, the structural bulk radius R (ospec_weight_stats)
%   - the eigenvalues of the Jacobian with every gain equal to 1 (the
%     slope of the piecewise nonlinearity on its linear part), as rates:
%     ev = (eig(W) - 1) / tau_d in 1/s, with the bulk circle's center
%     -1/tau_d and radius R / tau_d
%   - the representative weight image (weight_image_sample: plan.n_image
%     neurons, half E and half I, from the private stream plan.image_seed)
%   - for a network with added low-rank structure, where the addition was
%     built to put its eigenvalue (as a rate)
% and for each EXAMPLE only, a run of the unadapted network over [0, T]
% from the class's initial state with the preset's noise (OspecNet's
% run_sampled, the class's integrator, every plan.dt_sample seconds):
%   - the dendritic potential x of the first plan.n_trace E and the first
%     plan.n_trace I neurons, their setpoints and rates
%   - the population mean rates of E and of I
%   - lambda_1 over [lya_from, T] (Benettin), when plan.lya
%   - summaries over the second half of the run: mean rates, the fraction
%     of neuron-samples at the floor or the ceiling, the SD of x
%
% WRITES <run_dir>/concept_examples/concept_examples.mat with
%   examples   1 x m struct, the fields above (cx_load documents them)
%   gallery    1 x k struct, the same without the run's fields
%   settings   the plan without its lists
% A rerun into the same folder overwrites the file. The runs go one after
% another in the client: four networks of 500 neurons for 12 s are about a
% minute.
%
% See also: concept_examples_stage_plan, cx_load, cx_network, weight_image_sample

arguments
    cfg     (1,1) struct %#ok<INUSA>
    run_dir (1,:) char
    plan    (1,1) struct
end

verbose = plan.verbose;
out_dir = fullfile(run_dir, plan.folder);
if ~isfolder(out_dir); mkdir(out_dir); end
out_file = fullfile(out_dir, plan.file);
if isfile(out_file); delete(out_file); end

t0 = tic;
examples = cell(1, numel(plan.examples));
for k = 1:numel(plan.examples)
    e = plan.examples(k);
    [model, params, ws, net] = cx_network(plan, e, plan.T);
    r = structure_record(plan, e, model, params, ws, net);
    r = add_run(r, plan, model, params);
    examples{k} = r;
    vprintf(verbose, 'minimal', ['  concept_examples: %-12s %-26s level %.2f  rightmost %+.1f%+.1fi 1/s  R %.1f 1/s  ' ...
        'rates E %.2f I %.2f  lambda_1 %+.2f 1/s\n'], e.key, e.name, e.level_of_chaos, real(r.rightmost), ...
        abs(imag(r.rightmost)), r.R_rate, r.mean_rate_E, r.mean_rate_I, r.LLE);
end
gallery = cell(1, numel(plan.gallery));
for k = 1:numel(plan.gallery)
    e = plan.gallery(k);
    [model, params, ws, net] = cx_network(plan, e);
    gallery{k} = structure_record(plan, e, model, params, ws, net);
end
examples = [examples{:}];
gallery = [gallery{:}];
if isempty(gallery); gallery = struct([]); end
settings = rmfield(plan, {'examples', 'gallery'});
save(out_file, 'examples', 'gallery', 'settings', '-v7');
detail = sprintf('%s (%d examples run for %g s, %d gallery structures; %.0f s)', out_file, numel(examples), ...
    plan.T, numel(gallery), toc(t0));
vprintf(verbose, 'minimal', '  concept_examples: %s\n', detail);
end

%% ------------------------------------------------------------------------
function r = structure_record(plan, e, model, params, ws, net)
% What both lists keep: the weights' statistics and the equal-gain spectrum.
tau = model.tau_d;
lam = eig(ws.Wf);
ev = (lam - 1) / tau;
[~, q] = max(real(ev));
ti = ws.ti;
dens = zeros(2); nz_mean = zeros(2); nz_sd = zeros(2);
for a = 1:2
    for b = 1:2
        blk = ws.Wf(ti{a}, ti{b});
        nz = blk(blk ~= 0);
        dens(a, b) = numel(nz) / numel(blk);
        nz_mean(a, b) = mean(nz);
        nz_sd(a, b) = std(nz, 1);
    end
end
lowrank_rate = NaN;
if ~isempty(model.added); lowrank_rate = (model.added.target - 1) / tau; end
r = struct();
r.key = e.key;
r.name = e.name;
r.label = net.label;
r.target = net.target;
r.level_of_chaos = e.level_of_chaos;
r.seed = e.seed;
r.n = model.n;
r.n_E = numel(ti{1});
r.n_I = numel(ti{2});
r.alpha = params.alpha;
r.tau_d = tau;
r.ev = ev;                                  % (eig(W) - 1) / tau_d, 1/s
r.center = -1 / tau;
r.R_struct = ws.R_struct;                   % units of W
r.R_rate = ws.R_struct / tau;               % 1/s: the bulk circle's radius
r.R_theory_rate = model.R / tau;
r.rightmost = ev(q);
r.mu = ws.mu;                               % (post, pre), zeros included
r.s2 = ws.s2;
r.density = dens;                           % fraction of nonzero entries per block
r.nz_mean = nz_mean;
r.nz_sd = nz_sd;
r.K_struct = ws.K_struct;
r.K_eig_rate = (ws.K_struct_eig - 1) / tau;
r.lowrank_rate = lowrank_rate;              % NaN: nothing added
r.sample = weight_image_sample(ws.Wf, ws.type_of, plan.n_image, plan.image_seed);
end

function r = add_run(r, plan, model, params)
% The unadapted run: kept neurons, population rates, lambda_1.
step = round(plan.dt_sample * plan.fs);
nt = numel(model.t_ex);
keep = 1:step:nt;
[S, lya] = model.run_sampled(keep);
n = model.n;
x = S(1:n, :);
c = model.S_c;
if ~isempty(model.S_c_vec); c = model.S_c_vec(:); end
rate = model.activation_function(x);        % the class's own nonlinearity, each neuron at its setpoint
ti = params.type_indices;
tr = [ti{1}(1:plan.n_trace), ti{2}(1:plan.n_trace)];
r.t = model.t_ex(keep);
r.t = r.t(:)';
r.trace_idx = tr(:);
r.trace_type = [ones(plan.n_trace, 1); 2 * ones(plan.n_trace, 1)];
r.x_trace = x(tr, :);
r.rate_trace = rate(tr, :);
if isscalar(c); r.S_c_trace = c * ones(numel(tr), 1); else; r.S_c_trace = c(tr); end
r.S_a = model.S_a;
r.activation = model.activation;
r.rate_E = mean(rate(ti{1}, :), 1);
r.rate_I = mean(rate(ti{2}, :), 1);
late = r.t >= r.t(end) / 2;
r.mean_rate_E = mean(r.rate_E(late));
r.mean_rate_I = mean(r.rate_I(late));
rl = rate(:, late);
r.frac_floor = mean(rl(:) < 1e-6);
r.frac_ceiling = mean(rl(:) > 1 - 1e-6);
r.x_sd = mean(std(x(:, late), 0, 2));
r.LLE = NaN;
r.lya_interval = [NaN NaN];
if plan.lya
    r.LLE = lya.LLE;
    r.lya_interval = [plan.lya_from, plan.T];
end
end
