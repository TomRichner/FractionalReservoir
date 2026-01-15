# VAR_SRNN Implementation Notes

## Overview

This module applies Vector Autoregression (VAR) analysis to SRNN simulation output, using the same pipeline as the cscs_dynamics SEEG analysis. The goal is to compare stability metrics between simulated neural networks and real brain recordings.

---

## Status

### ✅ Complete

| Script | Purpose | Status |
|--------|---------|--------|
| `VAR_SRNN_comparison.m` | Simulate 4 conditions with Lyapunov | ✅ Working |
| `export_SRNN_to_CSCS_format.m` | Export to CSCS format | ✅ Working |
| `SRNN_VAR_subject_structure.m` | Generate subject struct | ✅ Working |
| `SRNN_load_and_preprocess.m` | Filter/downsample data | ✅ Working |
| `test_VAR_SRNN_export.m` | 7 unit tests | ✅ All pass |

### 🔄 To Do

- [ ] Run `CSCS_run_ICA.m` on preprocessed data  
- [ ] Run `CSCS_run_VAR.m` on ICA-cleaned data
- [ ] Compare eigenvalue distributions across conditions
- [ ] Validate expected stability pattern: none → SFA → STD → both

### 📋 Future Enhancements

1. **Higher Fs production runs** (5000 Hz simulation → 500 Hz after decimation)
2. **Add 1/f noise** to simulate distant neural activity
3. **Stimulation artifact modeling** for ICA validation
4. **Cross-project visualization** (SRNN vs SEEG eigenvalues)

---

## Simulation Conditions

| Condition | n_a_E | n_b_E | LLE (dev) | Expected Stability |
|-----------|-------|-------|-----------|-------------------|
| `SRNN_none` | 0 | 0 | -0.40 | Most unstable |
| `SRNN_SFA` | 3 | 0 | -0.04 | Moderately stable |
| `SRNN_STD` | 0 | 1 | — | Moderately stable |
| `SRNN_both` | 3 | 1 | — | Most stable |

---

## Data Flow

```
VAR_SRNN_comparison.m
    ↓ (simulate)
export_SRNN_to_CSCS_format.m
    ↓ (save raw .mat)
SRNN_load_and_preprocess.m
    ↓ (filter, detrend)
data/VAR_SRNN_processed/*.mat
    ↓
CSCS_run_ICA.m → CSCS_pick_ICA_components.m → CSCS_run_VAR.m
```

---

## Output Directories

- **Raw simulation**: `data/VAR_SRNN/`
- **Preprocessed**: `data/VAR_SRNN_processed/`
- **Figures**: `data/VAR_SRNN/figs/`

---

## Key Parameters

| Parameter | Dev Value | Prod Value |
|-----------|-----------|------------|
| `fs` | 200 Hz | 5000 Hz |
| `T_baseline` | 30 sec | 15 min |
| `T_stim` | 30 sec | 15 min |
| `Fs_final` | 200 Hz | 500 Hz |
| `deci_mode` | 'no_deci' | 'butter' |

---

## Testing

```matlab
% Run unit tests
cd scripts/VAR_SRNN
results = runtests('test_VAR_SRNN_export');
```

All 7 tests pass:
- Export format structure
- Data dimensions [time × channels]
- Channel count matches Files
- Sampling rate preservation
- Subject structure format
- Neuron subset option
- Baseline/stim split

---

## References

- `FractionalResevoir/Code_Structure.md`: SRNN model details
- `cscs_dynamics/Workflow.md`: 10-step SEEG analysis pipeline
- `cscs_dynamics/CSCSreview_loadfunction_revived.m`: Output format reference
