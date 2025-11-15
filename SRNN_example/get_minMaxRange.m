function min_max_range = get_minMaxRange(params)
    % Returns bounds per state variable for the SRNN state vector used in
    % benettin_algorithm. By default the SRNN has no hard bounds, so all
    % entries are initialized to NaN, but the vector must match the full
    % state layout: S = [a_E(:); a_I(:); b_E(:); b_I(:); x(:)].
    %
    % You can edit the sections below to set bounds for specific groups.
    
    % Extract sizes
    n       = params.n;
    n_E     = params.n_E;
    n_I     = params.n_I;
    n_a_E   = params.n_a_E;
    n_a_I   = params.n_a_I;
    n_b_E   = params.n_b_E;
    n_b_I   = params.n_b_I;
    
    len_a_E = n_E * n_a_E;
    len_a_I = n_I * n_a_I;
    len_b_E = n_E * n_b_E;
    len_b_I = n_I * n_b_I;
    len_x   = n;
    
    N_sys_eqs = len_a_E + len_a_I + len_b_E + len_b_I + len_x;
    
    % Initialize with NaN (no bounds by default)
    min_max_range = nan(N_sys_eqs, 2);
    
    % Example of setting bounds per group (commented out):
    % min_max_range(1:len_a_E, :) = [... lower upper ...];
end