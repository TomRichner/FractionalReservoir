function [t_out, S_out, params, lya_results, plot_data] = full_SRNN_run(u_ex_scale, n_a_E, n_b_E, level_of_chaos, rng_seeds, save_dir, save_figs, save_workspace, note, time_config)
% FULL_SRNN_RUN Runs the SRNN simulation with specified parameters.
%
% Inputs:
%   u_ex_scale    - Scaling factor for external input
%   n_a_E         - Number of adaptation time constants for E neurons
%   n_b_E         - Number of STD timescales for E neurons
%   level_of_chaos- Level of chaos parameter
%   rng_seeds     - Seeds for random number generation
%   save_dir      - Directory path to save outputs
%   save_figs     - Boolean, whether to save figures
%   save_workspace- Boolean, whether to save the workspace
%   note          - Optional string for naming figures (default: 'SRNN')
%   time_config   - Optional struct with fields:
%                   .T_range - [start, end] simulation time (default: [0, 50])
%                   .T_plot  - [start, end] plotting time (default: T_range)

if nargin < 9 || isempty(note)
    note = 'SRNN';
end

if nargin < 10 || isempty(time_config)
    time_config = struct();
end
if ~isfield(time_config, 'T_range')
    time_config.T_range = [0, 50];
end
if ~isfield(time_config, 'T_plot')
    time_config.T_plot = time_config.T_range;
end

% params.rng_seeds = [5 3 3 4 5];
params.rng_seeds = rng_seeds;
rng(params.rng_seeds(1))

% Lyapunov method selection
Lya_method = 'benettin'; % 'benettin', 'qr', or 'none'

%% time 
fs = 200;  % 100 Hz is good enough for 0.01, 200 Hz is good for 0.001 resolution of LLE Sampling frequency (Hz)
dt = 1/fs;
T_range = time_config.T_range;

%% input
% u_ex_scale passed as argument

% set simulation parameters
params.n = 100;
params.indegree = 20; % expected indegree
params.d = params.indegree/params.n; % expected density
params.f = 0.5; % fraction of neurons that are E
params.G_stdev = 1/sqrt(params.indegree); % fix this equation such that the expected spectral radius of sparse W is 1. 
params.mu_E = 1;
params.level_of_chaos = level_of_chaos; % Passed as argument

%% set adaptation params
params.n_a_E = n_a_E;  % Passed as argument
params.n_a_I = 0;  % 0 to 3 adaptation for I neurons
params.tau_a_E = logspace(log10(0.25), log10(10), params.n_a_E);  % Logarithmically spaced from 0.1 to 10
params.tau_a_I = logspace(log10(0.25), log10(10), params.n_a_I);  % Logarithmically spaced from 0.1 to 10
params.c_E = 0.1/3;  % Adaptation scaling for E neurons (scalar, typically 0-3)
params.c_I = .1;  % Adaptation scaling for I neurons (scalar, typically 0-3)

params.n_b_E = n_b_E;  % Passed as argument
params.n_b_I = 0;  % Number of STD timescales for I neurons (0 or 1)
params.tau_b_E_rel = 0.25;  % STD release time constant for E neurons (s)
params.tau_b_I_rel = 0.25;  % STD release time constant for I neurons (s)
params.tau_b_E_rec = 1;  % STD recovery time constant for E neurons (s)
params.tau_b_I_rec = 1;  % STD recovery time constant for I neurons (s)

%% set activaiton function
S_a = 0.9;
S_c = 0.35;
params.activation_function = @(x) piecewiseSigmoid(x, S_a, S_c);
params.activation_function_derivative = @(x) piecewiseSigmoidDerivative(x, S_a, S_c);

%% dependent E vs I params
params.n_E = round(params.f*params.n);
params.n_I = params.n-params.n_E;
params.N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + params.n_E * params.n_b_E + params.n_I * params.n_b_I + params.n;
params.E_indices = 1:params.n_E;
params.I_indices = params.n_E+1:params.n;
imbalance = 1;
params.mu_I = -imbalance*params.f*params.mu_E/(1-params.f);


%% make W connection matrix
[W, M, G, Z] = create_W_matrix(params); % should use params.d not recompute it
% Scale W based on abscissa to control level of chaos
W_eigs = eig(W);
radius_0 = max(abs(W_eigs))
abscissa_0 = max(real(W_eigs))

% % add a plot of eigs of W here
% figure;
% plot_eigenvalues(W_eigs, gca, 0);
% title('Eigenvalues of unscaled W');

gamma = 1 / abscissa_0;  % Scaling to reach edge of chaos
params.W = params.level_of_chaos * gamma * W;

% verify new abscissa and spectral radius
W_eigs_1 = eig(params.W);
radius_1 = max(abs(W_eigs_1))
abscissa_1 = max(real(W_eigs_1))

params.tau_d = 0.1;  % dentdritic time constant, seconds

%% ICs
S0 = initialize_state(params);

%% External input
% Configure external input parameters
input_config.n_steps = 3;
input_config.step_density = 0.2;  % Fraction of neurons receiving input at each step
input_config.amp = 0.5;  % Amplitude scaling factor
input_config.no_stim_pattern = false(1, input_config.n_steps);
input_config.no_stim_pattern(1:2:end) = true;  % No stimulation for odd steps
input_config.intrinsic_drive = -0 + 0*randn(params.n, 1);  % Intrinsic drive

% Generate external input
T_stim = T_range(2);
[u_stim, t_stim] = generate_external_input(params, T_stim, fs, params.rng_seeds(2), input_config);

if T_range(1) < 0
    % Prepend zeros for settling time
    t_pre = (T_range(1):dt:-dt)'; 
    u_pre = zeros(params.n, length(t_pre));
    
    t_ex = [t_pre; t_stim];
    u_ex = [u_pre, u_stim];
else
    % Slice if start time is positive
    indices = t_stim >= T_range(1);
    t_ex = t_stim(indices);
    u_ex = u_stim(:, indices);
end

u_ex = u_ex.*u_ex_scale; % way to turn it off easily
nt = length(t_ex);
step_period = fix(T_stim/input_config.n_steps);
n_steps = input_config.n_steps;

% Integrate
ode_solver = @ode45;
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
jac_wrapper = @(t, S) compute_Jacobian_fast(S, params);
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep', 1*dt, 'Jacobian', jac_wrapper);
% opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'MaxStep', 0.1*dt, 'Jacobian', jac_wrapper); % high precision comparison
% [t_out, S_out] = ode45(rhs, t_ex, S0, opts);
disp('Integrating Equations...')
tic
[t_out, S_out] = ode_solver(rhs, t_ex, S0, opts);
toc

%% Compute Lyapunov exponent(s)
T_interval = [T_range(1)+1, T_range(2)];
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
lya_results = compute_lyapunov_exponents(Lya_method, S_out, t_out, dt, fs, T_interval, params, opts, ode_solver, rhs, t_ex, u_ex);

%% Decimate state vector for plotting
plot_deci = round(fs/20);
[t_plot, S_plot, plot_indices] = decimate_states(t_out, S_out, plot_deci);

%% Unpack state vector and compute firing rates
[x_plot, a_plot, b_plot, r_plot, br_plot] = unpack_and_compute_states(S_plot, params);

% Split external input into E and I
u_ex_plot = u_ex(:, plot_indices);
u_plot.E = u_ex_plot(params.E_indices, :);
u_plot.I = u_ex_plot(params.I_indices, :);

%% Trim data to T_plot range
if isfield(time_config, 'T_plot') && ~isempty(time_config.T_plot)
    keep_mask = t_plot >= time_config.T_plot(1) & t_plot <= time_config.T_plot(2);
    
    t_plot = t_plot(keep_mask);
    
    % Trim u_plot (manually as it's not generated by unpack)
    if ~isempty(u_plot.E), u_plot.E = u_plot.E(:, keep_mask); end
    if ~isempty(u_plot.I), u_plot.I = u_plot.I(:, keep_mask); end
    
    % Trim other structs
    x_plot = trim_struct_data(x_plot, 2, keep_mask);
    r_plot = trim_struct_data(r_plot, 2, keep_mask);
    b_plot = trim_struct_data(b_plot, 2, keep_mask);
    br_plot = trim_struct_data(br_plot, 2, keep_mask);
    a_plot = trim_struct_data(a_plot, 3, keep_mask);
    
    % Trim Lyapunov results
    if ~isempty(lya_results) && isfield(lya_results, 't_lya')
        keep_mask_lya = lya_results.t_lya >= time_config.T_plot(1) & lya_results.t_lya <= time_config.T_plot(2);
        
        lya_results.t_lya = lya_results.t_lya(keep_mask_lya);
        
        if isfield(lya_results, 'local_lya')
             lya_results.local_lya = lya_results.local_lya(keep_mask_lya);
        end
        if isfield(lya_results, 'finite_lya')
             lya_results.finite_lya = lya_results.finite_lya(keep_mask_lya);
        end
        
        % Handle QR fields if present
        if isfield(lya_results, 'local_LE_spectrum_t')
            % Spectrum is [time x exponents]
            lya_results.local_LE_spectrum_t = lya_results.local_LE_spectrum_t(keep_mask_lya, :);
        end
        
        if isfield(lya_results, 'finite_LE_spectrum_t')
            lya_results.finite_LE_spectrum_t = lya_results.finite_LE_spectrum_t(keep_mask_lya, :);
        end
    end
end

%% Plotting
plot_data.t = t_plot;
plot_data.u = u_plot;
plot_data.x = x_plot;
plot_data.r = r_plot;
plot_data.a = a_plot;
plot_data.b = b_plot;
plot_data.br = br_plot;

[~, ~] = plot_SRNN_tseries(t_plot, u_plot, x_plot, r_plot, a_plot, b_plot, br_plot, params, lya_results, Lya_method, time_config.T_plot);

%% Compute Jacobian eigenvalues at multiple time points
% Sample at the center of each stimulus ON period
J_times_sec = [];
for k = 1:n_steps
    if ~input_config.no_stim_pattern(k)
        t_center = (k-1)*step_period + step_period/2;
        J_times_sec = [J_times_sec, t_center];
    end
end

% Convert times in seconds to indices (MATLAB is 1-indexed)
% Adjust for start time of simulation (t_out(1) corresponds to index 1)
J_times = round((J_times_sec - t_out(1)) * fs) + 1;
J_times = unique(J_times);  % Ensure unique indices and sort

fprintf('Computing Jacobian at %d time points...\n', length(J_times));
J_array = compute_Jacobian_at_indices(S_out, J_times, params);

% Compute eigenvalues for each Jacobian
eigenvalues_all = cell(length(J_times), 1);
for i = 1:length(J_times)
    eigenvalues_all{i} = eig(J_array(:,:,i));
end

%% Plot eigenvalue distributions on complex plane
n_plots = length(J_times);
if n_plots <= 4
    n_rows = 1;
    n_cols = n_plots;
else
    n_cols = ceil(sqrt(n_plots));
    n_rows = ceil(n_plots / n_cols);
end

figure('Position', [  1312         526         600         360]);
ax_handles = zeros(n_plots, 1);

% Compute global axis limits across all eigenvalue sets
all_real = [];
all_imag = [];
for i = 1:n_plots
    evals = eigenvalues_all{i};
    all_real = [all_real; real(evals)];
    all_imag = [all_imag; imag(evals)];
end
global_xlim = [min(all_real), max(all_real)];
global_ylim = [min(all_imag), max(all_imag)];

% Add some padding (10% of range on each side)
x_range = diff(global_xlim);
y_range = diff(global_ylim);
global_xlim = global_xlim + [-0.1, 0.1] * x_range;
global_ylim = global_ylim + [-0.1, 0.1] * y_range;

for i = 1:n_plots
    ax_handles(i) = subplot(n_rows, n_cols, i);
    evals = eigenvalues_all{i};
    time_val = t_out(J_times(i));
    ax_handles(i) = plot_eigenvalues(evals, ax_handles(i), time_val, global_xlim, global_ylim);
    set(ax_handles(i), 'Color', 'none');
    
    % Add time annotation in lower left corner
    % x_range = diff(global_xlim);
    % y_range = diff(global_ylim);
    % text_x = global_xlim(1) + 0.05 * x_range;
    % text_y = global_ylim(1) + 0.05 * y_range;
    % text(text_x, text_y, sprintf('t = %.2f s', time_val), ...
    %     'FontSize', 14, 'Color', 'k', 'BackgroundColor', 'white', ...
    %     'EdgeColor', 'none');
end

% Link axes of all subplots
linkaxes(ax_handles, 'xy');


%% new figure with imagesc plot of J_eff at J_times

fprintf('Computing J_eff at %d time points...\n', length(J_times));

omit_diagonal_in_J_eff = true;

% Compute J_eff for each time point
J_eff_array = zeros(params.n, params.n, length(J_times));
for i = 1:length(J_times)
    J_eff_array(:,:,i) = full(compute_J_eff(S_out(J_times(i),:)', params));
    
    % Optionally zero out diagonal to better visualize off-diagonal terms
    if omit_diagonal_in_J_eff
        J_eff_array(:,:,i) = J_eff_array(:,:,i) - diag(diag(J_eff_array(:,:,i)));
    end
end

% Compute global color limits across all J_eff matrices
if omit_diagonal_in_J_eff
    % Exclude diagonal when computing color limits
    all_J_eff_values = [];
    for i = 1:length(J_times)
        J_temp = J_eff_array(:,:,i);
        % Get off-diagonal values only
        all_J_eff_values = [all_J_eff_values; J_temp(~eye(size(J_temp)))];
    end
else
    all_J_eff_values = J_eff_array(:);
end

% Use 90th percentile of absolute values of non-zero elements
nonzero_abs_values = abs(all_J_eff_values(all_J_eff_values ~= 0));
max_abs = prctile(nonzero_abs_values, 90);
global_clim = [-max_abs, max_abs];

% Create figure with same layout as eigenvalue figure
figure('Position', [1312         940         600         310]);

for i = 1:n_plots
    subplot(n_rows, n_cols, i);
    
    % Plot J_eff matrix
    imagesc(J_eff_array(:,:,i));
    
    % Set colormap and color limits
    colormap(bluewhitered_colormap(256));
    clim(global_clim);
    
    % Make axes square
    axis square;
    
    % Remove ticks, labels, and box
    set(gca, 'XTick', [], 'YTick', []);
    box off;
    set(gca, 'Color', 'none');
    
    % Make axis lines white and send to back
    set(gca, 'XColor', 'white', 'YColor', 'white', 'Layer', 'bottom');
    
    % Add time annotation below the plot
    % For imagesc, y-axis is inverted (top=1, bottom=n), so below means y > max
    time_val = t_out(J_times(i));
    % ax_lim = axis;  % Get current axis limits [xmin xmax ymin ymax]
    % text_x = ax_lim(1);
    % text_y = ax_lim(4) + 0 * (ax_lim(4) - ax_lim(3));  % Below the plot
    % text(text_x, text_y, sprintf('t = %.2f s', time_val), ...
    %     'FontSize', 14, 'Color', 'k', ...
    %     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

%% Create separate colorbar figure for J_eff
figure('Position', [100   346   285   154], 'Color', 'white');
ax = axes('Position', [0.3, 0.1, 0.3, 0.8]);
colormap(bluewhitered_colormap(256));
cb = colorbar('Location', 'east');
clim(global_clim);
set(gca, 'Visible', 'off', 'Color', 'none');
set(cb, 'AxisLocation', 'out', 'Ticks', []);
ylabel(cb, 'J_{eff}', 'Interpreter', 'tex', 'FontSize', 14);

%% Plot J_eff as directed graph
figure('Position', [300, 400,    900,   400]);

for i = 1:n_plots
    subplot(n_rows, n_cols, i);
    
    % Plot J_eff graph
    % Reuse max_abs and global_clim from imagesc section
    plot_J_eff_graph(J_eff_array(:,:,i), max_abs, global_clim);
    set(gca, 'Color', 'none');
    
    % Add time annotation
    time_val = t_out(J_times(i));
    % text(-1.3, 1.3, sprintf('t = %.2f s', time_val), ...
    %     'FontSize', 14, 'Color', 'k', 'HorizontalAlignment', 'left');
end

%% Save results
if save_figs
    if ~isempty(save_dir)
        figs_folder = fullfile(save_dir, 'figs');
        if ~exist(figs_folder, 'dir')
            mkdir(figs_folder);
        end
        
        % Find all open figures
        fig_handles = findobj('Type', 'figure');
        % Sort by figure number
        [~, idx] = sort([fig_handles.Number]);
        fig_handles = fig_handles(idx);
        fig_vec = [fig_handles.Number];
        
        save_some_figs_to_folder_2(figs_folder, note, fig_vec, {'fig', 'svg', 'png', 'pdf'});
    else
        warning('save_figs is true but save_dir is empty. Figures not saved.');
    end
end

if save_workspace
    if ~isempty(save_dir)
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        save(fullfile(save_dir, 'SRNN_run.mat'), '-v7.3', '-nocompression');
    else
        warning('save_workspace is true but save_dir is empty. Workspace not saved.');
    end
end

end

function s_out = trim_struct_data(s_in, dim, mask)
    % Helper to trim fields of a struct along a dimension
    s_out = s_in;
    fields = fieldnames(s_in);
    for i = 1:length(fields)
        val = s_in.(fields{i});
        if ~isempty(val)
            if dim == 2
                s_out.(fields{i}) = val(:, mask);
            elseif dim == 3
                s_out.(fields{i}) = val(:, :, mask);
            end
        end
    end
end
