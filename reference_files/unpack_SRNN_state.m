function [a_E, a_I, b_E, b_I, u_d] = unpack_SRNN_state(X, params)
    % Unpack SRNN state vector/matrix X into named variables for E and I neurons.
    %
    % If X is N_sys_eqs x 1 (single time point, e.g., from within SRNN.m):
    %   a_E will be n_E x n_a_E (if n_a_E > 0, else [])
    %   a_I will be n_I x n_a_I (if n_a_I > 0, else [])
    %   b_E will be n_E x n_b_E (if n_b_E > 0, else [])
    %   b_I will be n_I x n_b_I (if n_b_I > 0, else [])
    %   u_d will be n x 1
    %
    % If X is nt x N_sys_eqs (multiple time points, e.g., output from ode solver):
    %   a_E will be n_E x n_a_E x nt (if n_a_E > 0, else [])
    %   a_I will be n_I x n_a_I x nt (if n_a_I > 0, else [])
    %   b_E will be n_E x n_b_E x nt (if n_b_E > 0, else [])
    %   b_I will be n_I x n_b_I x nt (if n_b_I > 0, else [])
    %   u_d will be n x nt

    n_E   = params.n_E;
    n_I   = params.n_I;
    n_a_E = params.n_a_E;
    n_a_I = params.n_a_I;
    n_b_E = params.n_b_E;
    n_b_I = params.n_b_I;
    n     = params.n;

    if iscolumn(X) % Single time point case, optimized for ODE solvers
        current_idx = 0;

        len_a_E = n_E * n_a_E;
        if len_a_E > 0
            a_E = reshape(X(current_idx + (1:len_a_E)), n_E, n_a_E);
        else
            a_E = [];
        end
        current_idx = current_idx + len_a_E;

        len_a_I = n_I * n_a_I;
        if len_a_I > 0
            a_I = reshape(X(current_idx + (1:len_a_I)), n_I, n_a_I);
        else
            a_I = [];
        end
        current_idx = current_idx + len_a_I;

        len_b_E = n_E * n_b_E;
        if len_b_E > 0
            b_E = reshape(X(current_idx + (1:len_b_E)), n_E, n_b_E);
        else
            b_E = [];
        end
        current_idx = current_idx + len_b_E;

        len_b_I = n_I * n_b_I;
        if len_b_I > 0
            b_I = reshape(X(current_idx + (1:len_b_I)), n_I, n_b_I);
        else
            b_I = [];
        end
        current_idx = current_idx + len_b_I;

        u_d = X(current_idx + (1:n));
    else % Multiple time points case
        nt = size(X, 1);
        current_idx = 0;

        len_a_E = n_E * n_a_E;
        if len_a_E > 0
            a_E = reshape(X(:, current_idx + (1:len_a_E))', n_E, n_a_E, nt);
        else
            a_E = [];
        end
        current_idx = current_idx + len_a_E;

        len_a_I = n_I * n_a_I;
        if len_a_I > 0
            a_I = reshape(X(:, current_idx + (1:len_a_I))', n_I, n_a_I, nt);
        else
            a_I = [];
        end
        current_idx = current_idx + len_a_I;

        len_b_E = n_E * n_b_E;
        if len_b_E > 0
            b_E = reshape(X(:, current_idx + (1:len_b_E))', n_E, n_b_E, nt);
        else
            b_E = [];
        end
        current_idx = current_idx + len_b_E;

        len_b_I = n_I * n_b_I;
        if len_b_I > 0
            b_I = reshape(X(:, current_idx + (1:len_b_I))', n_I, n_b_I, nt);
        else
            b_I = [];
        end
        current_idx = current_idx + len_b_I;

        u_d = X(:, current_idx + (1:n))';
    end
end