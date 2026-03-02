# Refactor: master_save_figs Override

**Date:** 2026-03-02

## Summary

Added `master_save_figs` support to `run_all_analyses.m` and all 4 analysis scripts, following the pattern from `ConnectivityAdaptation/run_all_figures.m`.

## How It Works

`run_all_analyses.m` sets a workspace variable:
```matlab
master_save_figs = 'save_all_figs';  % or 'save_no_figs' or 'follow_scripts_save_figs'
```

Each sub-script checks for it:
```matlab
if exist('master_save_figs', 'var')
    if strcmp(master_save_figs, 'save_all_figs'), save_figs = true;
    elseif strcmp(master_save_figs, 'save_no_figs'), save_figs = false;
    end
end
if ~exist('save_figs', 'var'), save_figs = false; end
```

When `save_figs = true`, scripts call `save_some_figs_to_folder_2` to save `.fig` and `.svg` to `<output_dir>/figures/`. PNGs are already auto-saved by `psa.plot_sensitivity` / `psa.plot`.

## Files Modified

| File                                   | Change                                                      |
| -------------------------------------- | ----------------------------------------------------------- |
| `run_all_analyses.m`                   | Added `master_save_figs = 'save_all_figs'`                  |
| `run_sensitivity_analysis.m`           | Added check block + save block after each parameter sweep   |
| `run_tau_sensitivity_analysis.m`       | Added check block + save blocks after tau_a and tau_b plots |
| `run_param_space_analysis2.m`          | Added check block + save block after `psa.plot`             |
| `Fig_2_fraction_excitatory_analysis.m` | Already had the check — no changes                          |
