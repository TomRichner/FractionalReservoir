# VAR_SRNN Implementation Notes

## Overview

This module applies Vector Autoregression (VAR) analysis to SRNN simulation output, using the same pipeline as the cscs_dynamics SEEG analysis. The goal is to compare stability metrics between simulated neural networks and real brain recordings.

## Current Implementation

### Simulation Conditions
Four SRNN conditions are treated as pseudo-subjects:
| Condition | Description | Expected Stability |
|-----------|-------------|-------------------|
| `SRNN_none` | No adaptation, no STD | Most unstable (eigenvalues near 1) |
| `SRNN_SFA` | Spike-frequency adaptation only | Moderately stable |
| `SRNN_STD` | Short-term depression only | Moderately stable |
| `SRNN_both` | Both SFA and STD | Most stable |

### Data Flow
1. **Simulate** → `VAR_SRNN_comparison.m` generates SRNN dynamics
2. **Export** → `export_SRNN_to_CSCS_format.m` saves as `.mat` files
3. **Preprocess** → `SRNN_VAR_load_simulated.m` applies filters and downsampling
4. **Analyze** → Standard cscs_dynamics pipeline (ICA → VAR → eigenvalue stats)

### Key Parameters
- **Neuron subset**: Option to use all 100 neurons or random subset (both E and I)
- **`max_minutes`**: Controls data duration for VAR fitting (matches `CSCS_run_VAR.m`)
- **ICA**: Run by default for consistency; `skip_ica` flag to disable

### Output Location
`FractionalResevoir/data/VAR_SRNN/`

---

## Future Enhancements

### 1. Add Realistic Noise (Priority: Medium)
Simulated data is "too clean" compared to SEEG. Future additions:
- **1/f (pink) noise**: Simulate distant neural activity
- **White noise floor**: Match recording system noise
- **Implementation**: Add noise before export in `export_SRNN_to_CSCS_format.m`

### 2. Stimulation Artifact Modeling (Priority: Low)
SEEG has stimulation artifacts that ICA removes. To test artifact removal:
- Add step/impulse artifacts at stimulation onset
- Verify ICA correctly identifies and removes them
- Compare VAR results with/without artifact removal

### 3. Extended Condition Comparisons (Priority: High)
After basic implementation works:
- Vary `level_of_chaos` parameter
- Explore different adaptation timescales
- Compare multiple network instantiations (different seeds)

### 4. Cross-Project Visualization (Priority: Medium)
Create unified plots comparing:
- Simulated vs SEEG eigenvalue distributions
- Effect of stimulation in both domains
- Lyapunov exponents (SRNN) vs VAR eigenvalues

---

## Testing Strategy

Use MATLAB testing framework via MCP server tools:
- `run_matlab_test_file`: For comprehensive test results
- `run_matlab_file`: For command window output
- `evaluate_matlab_code`: For inspecting variables (avoid large arrays)

---

## References
- `FractionalResevoir/Code_Structure.md`: SRNN implementation details
- `cscs_dynamics/Workflow.md`: 10-step SEEG analysis pipeline
