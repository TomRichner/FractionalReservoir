function params = package_params(n_E, n_I, E_indices, I_indices, n_a_E, n_a_I, n_b_E, n_b_I, tau_a_E, tau_a_I, tau_b_E, tau_b_I, tau_d, n, M, c_SFA, F_STD, tau_STD, EI_vec)
    % these params are sent into the ode solver via an anonymous function, nice to package them
    params.n_E = n_E;
    params.n_I = n_I;
    params.E_indices = E_indices;
    params.I_indices = I_indices;
    params.n_a_E = n_a_E;
    params.n_a_I = n_a_I;
    params.n_b_E = n_b_E;
    params.n_b_I = n_b_I;
    if params.n_a_E > 0
        params.tau_a_E = tau_a_E(:)';
    else
        params.tau_a_E = [];
    end
    if params.n_a_I > 0
        params.tau_a_I = tau_a_I(:)';
    else
        params.tau_a_I = [];
    end
    if params.n_b_E > 0
        params.tau_b_E = tau_b_E(:)';
    else
        params.tau_b_E = [];
    end
    if params.n_b_I > 0
        params.tau_b_I = tau_b_I(:)';
    else
        params.tau_b_I = [];
    end
    params.tau_d = tau_d;
    params.n = n;
    params.M = M;
    params.c_SFA = c_SFA(:);
    params.F_STD = F_STD(:);
    params.tau_STD = tau_STD;
    params.EI_vec = EI_vec(:);
end