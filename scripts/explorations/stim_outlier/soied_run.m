function soied_run(seeds, varargin)
% SOIED_RUN Discharges, stimulation and the frozen-gain coherent mode in Fig 8's IED network (pilot, 2026-11-21).
%
%   soied_run(42)            % one network, the bundle's first seed
%   soied_run([73 101])      % more networks (each saved on its own)
%   soied_run(seeds, 'block_len', 60, ...)   override a setting
%
% The question (extends scripts/explorations/stim_outlier, the
% bursting_pairs pilot, to the modular IED-like network of the paper's
% Figure 8, the ied_stim bundle): is each discharge preceded by a coherent
% mode of the frozen-gain connectivity M = W diag(d) crossing Re = 1, and
% what do the frozen-gain mode and the full Jacobian's network mode do
% under the constant uniform input 0.2?
%
% Each network: a baseline block and a stimulation block of block_len s
% after the bundle's 10-s warmup and 10-s settle (soied_run_network). The
% bundle analyses 300 s per block; this pilot analyses block_len (80 s) to
% fit its compute budget. Output: data/explorations/stim_outlier_ied/
% soied_seed<seed>.mat (a seed whose file exists is not re-simulated; the
% dense full-Jacobian pass around discharges, soied_event_full, is added
% when missing; it uses the block's own settings, with full_event_dt,
% n_event_full and coh_thr from this call). Measures: soied_measures; figure: soied_main_figure;
% sidecar: soied_write_sidecar.
%
% See also: soied_run_network, soied_measures, soied_main_figure, soied_write_sidecar

root = fileparts(which('setup_paths'));
S = struct('block_len', 80, 'settle', 10, 'fs', 500, 'amplitude', 0.20, ...
    'coarse_dt', 0.05, 'dense_dt', 0.01, 'dense_win', [-0.6 0.2], ...
    'full_dt', 1.0, 'n_top', 40, 'full_event_dt', 0.05, 'n_event_full', 6, 'coh_thr', 0.3, ...
    'data_dir', fullfile(root, 'data', 'explorations', 'stim_outlier_ied'), 'log', '');
for k = 1:2:numel(varargin), S.(varargin{k}) = varargin{k + 1}; end
if ~isfolder(S.data_dir), mkdir(S.data_dir); end
for s = seeds
    f = fullfile(S.data_dir, sprintf('soied_seed%03d.mat', s));
    if isfile(f)                    % resume: only what is missing
        R = load(f, 'R'); R = R.R;
    else
        R = soied_run_network(s, S);
        save(f, 'R', '-v7.3');
    end
    if ~isfield(R, 'event_full')
        t0 = tic;
        R.event_full = soied_event_full(s, R.block(1).events, R.settings, S);
        R.event_full.repro_ok = abs(R.event_full.repro_mean_x_E - ...
            mean(double(R.block(1).signals.mean_x_E))) < 1e-5;
        save(f, 'R', '-v7.3');
        if ~isempty(S.log)
            fid = fopen(S.log, 'a'); fprintf(fid, 'seed %d event_full in %.0f s\n', s, toc(t0)); fclose(fid);
        end
    end
end
end
