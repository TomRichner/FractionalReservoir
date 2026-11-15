classdef wiener_segments < handle
    %WIENER_SEGMENTS A run's Wiener increments, handed out a window of steps at a time.
    %
    %   w = WIENER_SEGMENTS(n_rows, n_steps, seed)
    %   w = WIENER_SEGMENTS(n_rows, n_steps, seed, schedule, t0, fs, err_id_prefix)
    %   [xi1, xi2] = w.next(n_cols)     % the next n_cols steps, in order
    %
    % THE SAME NUMBERS AS WIENER_INCREMENTS, bit for bit, without the two
    % n_rows x n_steps tensors: successive calls of next return successive
    % blocks of columns of the xi1 and xi2 that wiener_increments would
    % return for the same arguments (test_release_noise_reservoir asserts
    % isequal, with and without a schedule). For a run that is integrated
    % in segments and never looks back (the lean memory-capacity path of
    % SRNN_ESN_reservoir), where the whole tensors are 1.2 GB at n = 500
    % over 760 s at 200 Hz and a window is a few MB.
    %
    % HOW. wiener_increments draws a seed's path from one twister stream,
    % all of xi1 and then all of xi2. A twister stream cannot jump to a
    % column (a normal consumes a varying number of its uniforms), but it
    % can be read in order, and a matrix of normals is filled column by
    % column: so xi1 is read window by window from one stream of the seed,
    % and xi2 from a second stream of the same seed that has first been run
    % past the whole of xi1 (drawn in blocks and thrown away: about 1 s per
    % 1e8 normals, once). A noise_seed_schedule adds one such pair of
    % streams per row; every pair is read at every window, so each stays at
    % the run's current column, and a window takes its columns from the
    % seed that owns them (noise_schedule_columns, as wiener_increments
    % splices).
    %
    % ORDER ONLY. There is no going back and no skipping: whoever
    % re-integrates a stretch of the run against the same increments
    % (Benettin) needs wiener_increments.
    %
    % The global random stream is not touched.
    %
    % See also: wiener_increments, release_increments, sde_fixed_step,
    %           SRNN_ESN_reservoir.run_memory_capacity_lean

    properties (SetAccess = private)
        n_rows
        n_steps
        position = 0        % steps handed out so far
    end

    properties (Access = private)
        streams_xi1 = {}    % one stream per seed (the base seed, then each schedule row)
        streams_xi2 = {}
        first_col           % the first column each schedule row owns
    end

    methods
        function obj = wiener_segments(n_rows, n_steps, seed, schedule, t0, fs, err_id_prefix)
            if nargin < 4, schedule = []; end
            if nargin < 7 || isempty(err_id_prefix), err_id_prefix = 'wiener_segments'; end
            obj.n_rows = n_rows;
            obj.n_steps = n_steps;
            seeds = seed;
            obj.first_col = zeros(1, 0);
            if ~isempty(schedule)
                check_noise_seed_schedule(schedule, err_id_prefix);
                obj.first_col = noise_schedule_columns(schedule(:, 1), n_steps, t0, fs, err_id_prefix);
                seeds = [seed; schedule(:, 2)];
            end
            BLOCK = 4096;       % columns thrown away per draw while passing xi1
            for k = 1:numel(seeds)
                % The stream rng(seed, 'twister') makes: wiener_increments' draw.
                xi1_stream = RandStream('mt19937ar', 'Seed', seeds(k));
                xi2_stream = RandStream('mt19937ar', 'Seed', seeds(k));
                left = n_steps;
                while left > 0
                    take = min(left, BLOCK);
                    randn(xi2_stream, n_rows, take);
                    left = left - take;
                end
                obj.streams_xi1{k} = xi1_stream;
                obj.streams_xi2{k} = xi2_stream;
            end
        end

        function [xi1, xi2] = next(obj, n_cols)
            %NEXT The increments of the next n_cols steps of the run.
            cols = obj.position + (1:n_cols);
            if n_cols < 0 || obj.position + n_cols > obj.n_steps
                error('wiener_segments:OutOfRange', ...
                    'Asked for steps %d..%d of a run of %d steps.', ...
                    obj.position + 1, obj.position + n_cols, obj.n_steps);
            end
            xi1 = randn(obj.streams_xi1{1}, obj.n_rows, n_cols);
            xi2 = randn(obj.streams_xi2{1}, obj.n_rows, n_cols);
            for k = 1:numel(obj.first_col)
                % Read every scheduled path, so it stays at this column;
                % keep the columns its seed owns.
                s1 = randn(obj.streams_xi1{k + 1}, obj.n_rows, n_cols);
                s2 = randn(obj.streams_xi2{k + 1}, obj.n_rows, n_cols);
                last_col = obj.n_steps;
                if k < numel(obj.first_col), last_col = obj.first_col(k + 1) - 1; end
                owned = cols >= obj.first_col(k) & cols <= last_col;
                xi1(:, owned) = s1(:, owned);
                xi2(:, owned) = s2(:, owned);
            end
            obj.position = obj.position + n_cols;
        end
    end
end
