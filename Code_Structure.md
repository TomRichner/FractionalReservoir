# SRNN Reservoir Code Structure Documentation

## SRNN_reservoir.m

### Overview
`SRNN_reservoir.m` implements a rate network with spike-frequency adaptation (SFA) for use with MATLAB's ODE solvers (e.g., `ode45`). The function computes the time derivatives of all state variables according to the specified dynamical equations.

### Function Signature
```matlab
function [dS_dt] = SRNN_reservoir(t, S, t_ex, u_ex, params)
```

### Mathematical Model

The function implements the following system of differential equations:

```
dx_i/dt = -x_i/τ_d + Σ_j w_ij * r_j + u_i

r_i = φ(x_i - Σ_k a_{i,k})

da_{i,k}/dt = (-a_{i,k} + r_i) / τ_k
```

Where:
- `x_i`: Dendritic state of neuron i
- `r_i`: Firing rate of neuron i
- `a_{i,k}`: k-th adaptation variable for neuron i
- `φ`: Nonlinear activation function (e.g., ReLU, tanh)
- `w_ij`: Connection weight from neuron j to neuron i
- `u_i`: External input to neuron i
- `τ_d`: Dendritic time constant (scalar)
- `τ_k`: k-th adaptation time constant

---

## Input Arguments

### `t` (scalar)
Current time point requested by the ODE solver.

### `S` (N_sys_eqs × 1 column vector)
State vector containing all dynamic variables. **Total size**: `N_sys_eqs = n_E*n_a_E + n_I*n_a_I + n`

**Organization**: `S = [a_E(:); a_I(:); x(:)]`

### `t_ex` (nt × 1 column vector)
Time vector for external input, where `nt` is the number of time points.

### `u_ex` (n × nt matrix)
External input stimulus matrix, where each column corresponds to a time point in `t_ex`.

### `params` (struct)
Parameter structure containing all network parameters (see detailed breakdown below).

---

## Output

### `dS_dt` (N_sys_eqs × 1 column vector)
Time derivatives of all state variables, organized identically to `S`.

**Organization**: `dS_dt = [da_E_dt(:); da_I_dt(:); dx_dt]`

---

## State Vector Organization

### Full State Vector S
```
S = [a_E(:);      % Excitatory adaptation states (n_E * n_a_E × 1)
     a_I(:);      % Inhibitory adaptation states (n_I * n_a_I × 1)
     x]           % Dendritic states for all neurons (n × 1)
```

### Unpacked State Variables

#### `a_E` (n_E × n_a_E matrix, or empty if n_a_E = 0)
Adaptation variables for excitatory neurons.
- **Rows**: Individual excitatory neurons (indexed 1 to n_E)
- **Columns**: Different adaptation time constants (indexed 1 to n_a_E)
- **Example**: `a_E(i, k)` is the k-th adaptation variable for the i-th excitatory neuron

#### `a_I` (n_I × n_a_I matrix, or empty if n_a_I = 0)
Adaptation variables for inhibitory neurons.
- **Rows**: Individual inhibitory neurons (indexed 1 to n_I)
- **Columns**: Different adaptation time constants (indexed 1 to n_a_I)
- **Example**: `a_I(i, k)` is the k-th adaptation variable for the i-th inhibitory neuron

#### `x` (n × 1 column vector)
Dendritic states for all neurons (both excitatory and inhibitory).
- **Size**: n = n_E + n_I
- **Indexing**: First n_E elements are excitatory, next n_I elements are inhibitory

---

## Parameters Structure (`params`)

### Network Size Parameters

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `n` | scalar | 1×1 | Total number of neurons (n = n_E + n_I) |
| `n_E` | scalar | 1×1 | Number of excitatory neurons |
| `n_I` | scalar | 1×1 | Number of inhibitory neurons |
| `n_a_E` | scalar | 1×1 | Number of adaptation time constants for E neurons (0 to disable) |
| `n_a_I` | scalar | 1×1 | Number of adaptation time constants for I neurons (0 to disable) |

### Neuron Indices

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `E_indices` | vector | n_E×1 | Indices of excitatory neurons in the full network |
| `I_indices` | vector | n_I×1 | Indices of inhibitory neurons in the full network |

### Connectivity

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `W` | matrix | n×n | Connection weight matrix. `W(i,j)` is the weight from neuron j to neuron i. Should obey Dale's law (excitatory neurons have non-negative outgoing weights, inhibitory neurons have non-positive outgoing weights). |

### Time Constants

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `tau_d` | scalar | 1×1 | Dendritic time constant (seconds) |
| `tau_a_E` | vector | 1×n_a_E | Adaptation time constants for E neurons (seconds). Can be empty if n_a_E = 0. |
| `tau_a_I` | vector | 1×n_a_I | Adaptation time constants for I neurons (seconds). Can be empty if n_a_I = 0. |

### Activation Function

| Field | Type | Size | Description |
|-------|------|------|-------------|
| `activation_function` | function handle | - | Nonlinear activation function φ. Should accept a vector and return a vector of the same size. **Default**: `@(x) tanh(x)`. Common alternatives: `@(x) max(0, x)` (ReLU), `@(x) min(max(x, 0), 1)` (saturated ReLU). |

---

## Internal Variables

### Intermediate Computation Variables

#### `u` (n × 1 column vector)
External input at current time `t`, obtained by interpolating `u_ex` at time `t`.

#### `x_eff` (n × 1 column vector)
Effective dendritic potential after subtracting adaptation:
```matlab
x_eff(E_indices) = x(E_indices) - sum(a_E, 2)  % For E neurons
x_eff(I_indices) = x(I_indices) - sum(a_I, 2)  % For I neurons
```

#### `r` (n × 1 column vector)
Firing rate of all neurons: `r = activation_function(x_eff)`

### Derivative Variables

#### `dx_dt` (n × 1 column vector)
Time derivative of dendritic states:
```matlab
dx_dt = -x / tau_d + W * r + u
```

#### `da_E_dt` (n_E × n_a_E matrix, or empty)
Time derivatives of excitatory adaptation variables:
```matlab
da_E_dt = (r(E_indices) - a_E) ./ tau_a_E
```
Uses MATLAB broadcasting: `r(E_indices)` is n_E×1, `a_E` is n_E×n_a_E, `tau_a_E` is 1×n_a_E.

#### `da_I_dt` (n_I × n_a_I matrix, or empty)
Time derivatives of inhibitory adaptation variables:
```matlab
da_I_dt = (r(I_indices) - a_I) ./ tau_a_I
```
Uses MATLAB broadcasting: `r(I_indices)` is n_I×1, `a_I` is n_I×n_a_I, `tau_a_I` is 1×n_a_I.

---

## Performance Optimizations

### Persistent Variables
The function uses persistent variables for efficient input interpolation:
- `u_interpolant`: A `griddedInterpolant` object for fast interpolation of `u_ex`
- `t_ex_last`: Stores the last time vector to detect changes and rebuild the interpolant only when necessary

This avoids repeatedly creating the interpolant object on every function call during ODE integration.

---

## Usage Example

```matlab
% Setup parameters
params.n = 100;
params.n_E = 80;
params.n_I = 20;
params.E_indices = 1:80;
params.I_indices = 81:100;
params.n_a_E = 2;  % Two adaptation time constants for E neurons
params.n_a_I = 0;  % No adaptation for I neurons
params.W = randn(100, 100) * 0.1;  % Example connectivity
params.tau_d = 0.01;  % 10 ms
params.tau_a_E = [0.1, 1.0];  % Fast and slow adaptation
params.tau_a_I = [];  % Empty since n_a_I = 0
params.activation_function = @(x) tanh(x);

% Initial conditions
N_sys_eqs = params.n_E * params.n_a_E + params.n_I * params.n_a_I + params.n;
S0 = randn(N_sys_eqs, 1) * 0.01;  % Small random initial state

% External input
fs = 1000;  % Sampling frequency (Hz)
T = 1.0;    % Duration (s)
t_ex = linspace(0, T, T*fs+1)';
u_ex = randn(params.n, length(t_ex)) * 0.5;  % Random input

% Integrate
rhs = @(t, S) SRNN_reservoir(t, S, t_ex, u_ex, params);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[t_out, S_out] = ode45(rhs, t_ex, S0, opts);
```

---

## Notes

1. **Dale's Law**: The connection matrix `W` should respect Dale's law (excitatory neurons only make excitatory connections, inhibitory neurons only make inhibitory connections).

2. **Disabling Adaptation**: Set `n_a_E = 0` or `n_a_I = 0` to disable adaptation for excitatory or inhibitory neurons, respectively.

3. **State Vector Size**: The total size of the state vector is:
   ```
   N_sys_eqs = n_E * n_a_E + n_I * n_a_I + n
   ```

4. **Broadcasting**: The adaptation dynamics use MATLAB's implicit broadcasting to efficiently compute element-wise operations across multiple time constants.

5. **Interpolation**: External input is linearly interpolated between time points in `t_ex`. The interpolation method uses `'none'` for extrapolation, which returns NaN for out-of-bounds queries to catch errors if the ODE solver attempts to step outside the defined time range.

---

## SRNN_reservoir_caller.m

### Overview
`SRNN_reservoir_caller.m` is a complete example script demonstrating how to use `SRNN_reservoir.m` with MATLAB's `ode45` solver. It provides a simplified workflow compared to the full `SRNN_basic_example.m` (which includes Lyapunov exponent calculations), making it ideal for getting started with the SRNN reservoir model.

### Script Structure

#### 1. Initialization
- Clears workspace and sets random seed for reproducibility
- Uses `clear all` to reset persistent variables in `SRNN_reservoir.m`

#### 2. Network Setup
Uses helper functions from `reference_files/`:
- **`generate_M_no_iso(n, w, sparsity, EI)`**: Creates a sparse, strongly-connected connectivity matrix
  - `n = 10`: Total number of neurons
  - `EI = 0.7`: 70% excitatory, 30% inhibitory
  - `sparsity`: Controlled by mean in/out degree (5 connections)
  - `w`: Structure defining weight scaling for EE, EI, IE, II connections
- **`get_EI_indices(EI_vec)`**: Extracts indices of excitatory and inhibitory neurons

#### 3. Time Parameters
- Sampling frequency: `fs = 1000` Hz (1 ms resolution)
- Time interval: `T = [-2, 5]` seconds (negative time allows for initial transients)
- Time vector: `t = linspace(T(1), T(2), nt)'` (column vector)

#### 4. External Input Design
The script creates a 2-component input:
1. **Stimulus**: Sine wave applied to neuron 1
   - Baseline: 0.5, Amplitude: 0.5
   - Frequency: 1 Hz, Duration: 2 seconds
   - Starts at t = 1 second
2. **DC Offset**: Constant background input (0.1)
   - Ramps up smoothly over 1.5 seconds to avoid initial transients
   - Prevents network from settling at zero activity

#### 5. Adaptation Configuration
- **Excitatory neurons**: 3 adaptation time constants
  - `tau_a_E = logspace(log10(0.3), log10(15), 3)` 
  - Spans from 0.3s (fast) to 15s (slow)
- **Inhibitory neurons**: No adaptation (`n_a_I = 0`)
- **Dendritic time constant**: `tau_d = 0.025` seconds (25 ms)

#### 6. Integration with ode45
- Wraps `SRNN_reservoir` to include extra parameters: `t`, `u_ex`, `params`
- ODE options: `RelTol = 1e-6`, `AbsTol = 1e-8`
- Uses `ode45` (non-stiff solver, good for most networks)

#### 7. Post-Processing
Manually unpacks the state vector since `SRNN_reservoir` uses a simplified state organization:
- State organization: `S = [a_E(:); a_I(:); x(:)]`
- Extracts:
  - `a_E_ts`: n_E × n_a_E × nt (adaptation variables for E neurons)
  - `a_I_ts`: n_I × n_a_I × nt (adaptation variables for I neurons, empty if n_a_I = 0)
  - `x_ts`: n × nt (dendritic states)
- Computes firing rates using `compute_dependent_variables` (passes empty arrays for `b_E_ts`, `b_I_ts` since synaptic depression is not included)

#### 8. Visualization
Creates a 3-panel figure:
1. **External Input**: Shows `u_ex` for neurons 1-2
2. **Firing Rates**: Plots `r(t)` for all neurons
3. **Adaptation Variables**: Shows all adaptation variables for the first E neuron with time constants in legend

### Key Parameter Choices and Rationale

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `n = 10` | Small network | Fast simulation, easy visualization |
| `EI = 0.7` | 70% excitatory | Biologically realistic ratio |
| `mean_in_out_degree = 5` | Moderate connectivity | Ensures strong connectivity without being fully connected |
| `scale = 0.5/0.79782` | Weight scaling | Provides stable dynamics (not too weak, not too strong) |
| `tau_d = 0.025` s | Fast dendritic time constant | Typical for rate models, allows quick responses |
| `n_a_E = 3` | Multiple timescales | Captures rich adaptation dynamics |
| `DC = 0.1` | Small positive offset | Keeps network in sensitive regime |

### Modifying for Different Configurations

#### Change Network Size
```matlab
n = 100;  % Larger network
mean_in_out_degree = 10;  % Scale connectivity accordingly
```

#### Disable Adaptation
```matlab
n_a_E = 0;  % No adaptation
n_a_I = 0;
```

#### Use ReLU Activation Instead of Tanh
```matlab
params.activation_function = @(x) max(0, x);
```

#### Change Stimulus Pattern
```matlab
% Square wave instead of sine
u_ex(1, stim_start_idx:stim_end_idx) = stim_b0 + amp * sign(sin(2*pi*f_sin*t_stim));

% Multiple neurons stimulated
u_ex(1:3, stim_start_idx:stim_end_idx) = stim_b0 + amp * sin(2*pi*f_sin*t_stim);
```

#### Use Stiff Solver (for larger networks)
```matlab
[t_out, S_out] = ode15s(SRNN_wrapper, t, S0, ode_options);
```

### Expected Outputs

1. **Console Output**:
   - Integration progress message
   - Completion message with simulation time
   - Simulation speed (ratio of compute time to real time)

2. **Figure**:
   - Three-panel visualization showing input, firing rates, and adaptation
   - Clear legends and labels for interpretation

3. **Workspace Variables**:
   - `S_out`: Full state trajectory (nt × N_sys_eqs)
   - `t_out`: Output time vector from ODE solver
   - `r_ts`: Firing rates (n × nt)
   - `a_E_ts`, `a_I_ts`, `x_ts`: Unpacked state variables
   - `params`: Parameter structure (useful for further analysis)

### Dependencies

The script requires the following helper functions (located in `reference_files/`):
- `generate_M_no_iso.m`: Network connectivity generation
- `get_EI_indices.m`: Excitatory/inhibitory indexing
- `compute_dependent_variables.m`: Firing rate computation

Note: The script adds `../reference_files/` to the path automatically.

