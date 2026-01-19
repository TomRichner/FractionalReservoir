import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm, qr, pinv

# Set random seed for reproducibility
np.random.seed(42)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

def rk4_step_linear(A, x, dt):
    """Runge-Kutta 4 integration step for linear ODE dx/dt = Ax."""
    k1 = A @ x
    k2 = A @ (x + 0.5 * dt * k1)
    k3 = A @ (x + 0.5 * dt * k2)
    k4 = A @ (x + dt * k3)
    return x + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

def gram_schmidt_qr(V):
    """
    QR decomposition with a check to ensure positive diagonal R elements.
    This ensures the 'log' accumulation is consistent.
    """
    Q, R = qr(V)
    d = np.diag(np.sign(np.diag(R)))
    Q = np.dot(Q, d)
    R = np.dot(d, R)
    return Q, R

def generate_random_matrix_edge_of_stability(N, target_radius=1.0):
    """Generates a random matrix with max real eigenvalue approx 'target_radius - 1'."""
    # Start with random Gaussian matrix
    J = np.random.randn(N, N) / np.sqrt(N)
    # A = -I + gJ. If g ~ 1, max real part is near 0.
    # We adjust 'g' (gain) to match the target stability
    A = -np.eye(N) + target_radius * J
    return A

# ==============================================================================
# SECTION 1: VERIFICATION OF THEORETICAL DERIVATION
# Compares "Piecewise Formula" vs "Continuous RK4 Integration"
# ==============================================================================
print("-" * 60)
print("SECTION 1: Verification (Theory vs Numerical Simulation)")
print("-" * 60)

# Parameters
N_DIM = 20
DT = 0.01             # 100 Hz integration
SEGMENT_LEN = 0.5     # 500 ms
N_SEGMENTS = 100      # Total time = 50s
STEPS_PER_SEG = int(SEGMENT_LEN / DT)

# 1. Generate Sequence of LTI Systems
# We vary stability slowly to mimic a brain state transition
A_sequence = []
stability_profile = []

for i in range(N_SEGMENTS):
    # Oscillate between stable (0.8) and unstable (1.5)
    g = 1.15 + 0.5 * np.sin(i * 0.1) 
    A = generate_random_matrix_edge_of_stability(N_DIM, target_radius=g)
    A_sequence.append(A)
    stability_profile.append(g)

# --- Method A: Theoretical Calculation (Derived Formula) ---
# LE = (1/T) * Sum( log( diag( QR( Product( exp(A_k*tau) ) ) ) ) )
print("Computing Theoretical LEs (using Matrix Exponentials)...")
Q_theory = np.eye(N_DIM)
log_sum_theory = np.zeros(N_DIM)

for A in A_sequence:
    # 1. Exact Propagator for the segment
    M = expm(A * SEGMENT_LEN)
    
    # 2. Evolve Tangent Space
    Z = M @ Q_theory
    
    # 3. QR Decomposition
    Q_theory, R = gram_schmidt_qr(Z)
    
    # 4. Accumulate
    log_sum_theory += np.log(np.abs(np.diag(R)))

le_theory = log_sum_theory / (N_SEGMENTS * SEGMENT_LEN)
le_theory = np.sort(le_theory)[::-1] # Sort descending

# --- Method B: Numerical Integration (Simulation) ---
# Integrating dy/dt = A(t)y and dQ/dt = A(t)Q step-by-step
print("Computing Numerical LEs (using RK4 Integration)...")
Q_sim = np.eye(N_DIM)
log_sum_sim = np.zeros(N_DIM)

for A in A_sequence:
    for _ in range(STEPS_PER_SEG):
        # Integrate Tangent Space Matrix Q
        # dQ/dt = A * Q
        Q_sim = rk4_step_linear(A, Q_sim, DT)
        
        # Benettin's Method: Re-orthogonalize frequently (every step for stability)
        # Note: Mathematically, doing it at the end of the segment is equivalent,
        # but numerically, frequent QR is safer.
        Q_sim, R = gram_schmidt_qr(Q_sim)
        log_sum_sim += np.log(np.abs(np.diag(R)))

le_sim = log_sum_sim / (N_SEGMENTS * SEGMENT_LEN)
le_sim = np.sort(le_sim)[::-1]

# Comparison
print(f"\nAll LEs Comparison (N={N_DIM}):")
print(f"{'Idx':<4} | {'Theory':<12} | {'Simulation':<12} | {'Error':<10}")
for i in range(N_DIM):
    err = abs(le_theory[i] - le_sim[i])
    print(f"{i:<4} | {le_theory[i]:<12.5f} | {le_sim[i]:<12.5f} | {err:<10.2e}")

max_err = np.max(np.abs(le_theory - le_sim))
if max_err < 1e-3:
    print(f"\n>> SUCCESS: Theory matches Simulation (Max Err: {max_err:.2e})")
else:
    print(f"\n>> WARNING: Discrepancy detected (Max Err: {max_err:.2e})")


# ==============================================================================
# SECTION 2: SOMPOLINSKY NETWORK EXPLORATION
# Simulate Nonlinear -> Fit LTI -> Compute LEs
# ==============================================================================
print("\n" + "-" * 60)
print("SECTION 2: Sompolinsky Network (Nonlinear -> Fitted LTI)")
print("-" * 60)

# Parameters
N_NET = 30
TAU = 0.1          # Neuron time constant (100ms)
T_TOTAL = 60.0     # 60 seconds
STEPS = int(T_TOTAL / DT)
WINDOW_PTS = int(SEGMENT_LEN / DT) # Points per 500ms window

# Fixed random connectivity
J = np.random.randn(N_NET, N_NET) / np.sqrt(N_NET)

# Slowly varying 'Chaos Parameter' g(t)
t_space = np.linspace(0, T_TOTAL, STEPS)
g_t = 1.6 + 0.8 * np.sin(2 * np.pi * t_space / 20.0) # Period 20s

# Storage
X_data = np.zeros((N_NET, STEPS))
x = np.random.randn(N_NET) * 0.1

# Track True LEs (Jacobian) for reference
Q_true = np.eye(N_NET)
log_r_true = np.zeros(N_NET)
le_hist_true = []

print("Simulating Nonlinear Network and computing True LEs...")

for i in range(STEPS):
    g = g_t[i]
    phi = np.tanh(x)
    
    # 1. Update State (RK4 on Nonlinear System)
    # dx/dt = (-x + g * J * tanh(x)) / tau
    def f_nonlin(y): 
        return (-y + g * (J @ np.tanh(y))) / TAU
    
    k1 = f_nonlin(x)
    k2 = f_nonlin(x + 0.5*DT*k1)
    k3 = f_nonlin(x + 0.5*DT*k2)
    k4 = f_nonlin(x + DT*k3)
    x = x + (DT/6)*(k1 + 2*k2 + 2*k3 + k4)
    X_data[:, i] = x
    
    # 2. Update True LEs (Jacobian Evolution)
    # Jac = (-I + g * J * diag(phi')) / tau
    d_phi = 1.0 - phi**2
    Jac = (-np.eye(N_NET) + g * (J * d_phi)) / TAU
    
    # Simple Euler update for Q (sufficient for reference)
    Q_true = Q_true + (Jac @ Q_true) * DT
    Q_true, R = gram_schmidt_qr(Q_true)
    log_r_true += np.log(np.abs(np.diag(R)))
    
    if (i+1) % WINDOW_PTS == 0:
        le_hist_true.append(log_r_true / ((i+1)*DT))

# ------------------------------------------------------------------------------
# FITTING AND ANALYSIS
# ------------------------------------------------------------------------------
print("Fitting LTI models to 500ms segments...")

le_hist_fit = []
n_windows = int(STEPS / WINDOW_PTS)
Q_fit = np.eye(N_NET)
log_r_fit = np.zeros(N_NET)

for w in range(n_windows):
    # Extract Segment
    start = w * WINDOW_PTS
    end = start + WINDOW_PTS
    segment = X_data[:, start:end]
    
    # Fit VAR(1) Model: x(t+1) = M * x(t)
    # X_next = M * X_curr
    X_curr = segment[:, :-1]
    X_next = segment[:, 1:]
    
    # Least Squares Solution
    # M_fit represents the propagator for ONE time step (dt)
    M_fit = X_next @ pinv(X_curr)
    
    # Compute LEs for this Fitted Model
    # We apply M_fit 'WINDOW_PTS' times to advance time by 500ms
    for _ in range(WINDOW_PTS):
        Q_fit = M_fit @ Q_fit
        Q_fit, R = gram_schmidt_qr(Q_fit)
        log_r_fit += np.log(np.abs(np.diag(R)))
        
    le_hist_fit.append(log_r_fit / ((w+1)*SEGMENT_LEN))

# Convert to arrays for plotting
le_hist_true = np.array(le_hist_true)
le_hist_fit = np.array(le_hist_fit)

# Final Spectrum Comparison
final_true = np.sort(le_hist_true[-1])[::-1]
final_fit = np.sort(le_hist_fit[-1])[::-1]

print(f"\nFinal Spectrum (All {N_NET} Lyapunov Exponents):")
print(f"True Nonlinear:\n{final_true}")
print(f"Fitted LTI:\n{final_fit}")

# ==============================================================================
# PLOTTING
# ==============================================================================
plt.figure(figsize=(12, 10))

# Plot 1: Neural Activity
plt.subplot(3, 1, 1)
t_arr = np.linspace(0, T_TOTAL, STEPS)
for i in range(20):
    plt.plot(t_arr, X_data[i, :], alpha=0.5, linewidth=1.0)
plt.title("Sompolinsky Network Activity (First 20 Neurons)")
plt.ylabel("Firing Rate (tanh)")
# plt.legend(loc="upper right")

# Plot 2: Chaos Parameter
plt.subplot(3, 1, 2)
time_axis = np.linspace(0, T_TOTAL, len(le_hist_true))
plt.plot(time_axis, g_t[::WINDOW_PTS], 'r', label="Chaos Param g(t)")
plt.axhline(1.0, color='k', linestyle='--', alpha=0.5, label="Edge of Stability")
plt.ylabel("Gain g")
plt.legend(loc="upper right")
plt.title("Time-Varying Stability Parameter")

# Plot 3: LE Convergence
plt.subplot(3, 1, 3)
plt.plot(time_axis, le_hist_true[:, 0], 'b-', label="True LE 1")
plt.plot(time_axis, le_hist_fit[:, 0], 'c--', label="Fitted LE 1")
plt.plot(time_axis, le_hist_true[:, 1], 'g-', label="True LE 2")
plt.plot(time_axis, le_hist_fit[:, 1], 'y--', label="Fitted LE 2")
plt.xlabel("Time (s)")
plt.ylabel("Lyapunov Exponent")
plt.title("Convergence of Fitted vs True Lyapunov Exponents")
plt.legend()

plt.tight_layout()
plt.show()