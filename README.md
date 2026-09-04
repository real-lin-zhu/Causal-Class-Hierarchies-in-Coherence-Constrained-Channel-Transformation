# Causal-Class Hierarchies in Coherence-Constrained Channel Transformation

Minimal MATLAB/SeDuMi package for the numerical results in the paper by Lin Zhu, Benchi Zhao, Xuanqiang Zhao, Ranyiliu Chen, Shenggen Zheng, and Xin Wang.

This folder covers the amplitude-damping sweep and the 1000-pair random-channel study. The workspace did not contain the full-precision channel pair and solver output for the separate MISC--DISC example in Table 3, so the rounded matrices printed in the paper have not been substituted for that missing data.

## Requirements

- MATLAB (tested with R2025a)
- SeDuMi on the MATLAB path

CVX is not required. The `quantinf/` directory contains only the GPL-2 QuantInf functions used to generate the seeded random qubit channels; its license is included as `quantinf/COPYING`.

## Files

- `reduced_misc_primal_sdp.m`: reduced SeDuMi primal sweep for two amplitude-damping channels and the identity target at `epsilon = (0:19)/20`.
- `reduced_misc_dual_sdp.m`: corresponding reduced SeDuMi dual sweep.
- `solve_reduced_misc_sedumi.m`: reduced MISC primal/dual solver for an arbitrary pair of qubit resource channels.
- `run_random_misc_seq_ico_gap_distribution.m`: seeded 1000-pair fixed-order versus general-process experiment. The internal label `ICO` denotes the paper's general feasible class `Gen`; it does not assert a quantum-switch realization.
- `verify_saved_results.m`: checks the numerical statements directly from the supplied `.mat` files.
- `results/data/`: compact solver outputs used for the reported checks.
- `results/figures/Figure_combined.png`: the combined numerical figure used in the manuscript.

## Verify the supplied results

From this directory in MATLAB, run:

```matlab
verify_saved_results
```

The supplied data give:

- amplitude-damping primal/dual maximum absolute mismatch: `8.37455987889868e-09`;
- random pairs completed: `1000`, with no failed instances;
- fixed-order-to-general gaps at least `1e-5`: `989`;
- instances below `1e-5`: `11`;
- random-study primal/dual maximum absolute mismatch: `1.89943245090807e-08`.

## Recompute

Amplitude-damping sweep:

```matlab
reduced_misc_primal_sdp
reduced_misc_dual_sdp
```

Full seeded 1000-pair primal/dual run:

```matlab
setenv('RANDOM_MISC_NUM_SETS', '1000');
setenv('RANDOM_MISC_SOLVE_ALL_DUALS', 'true');
setenv('RANDOM_MISC_NUM_DUAL_SAMPLES', '5');
setenv('RANDOM_MISC_PROGRESS_EVERY', '20');
setenv('RANDOM_MISC_OUTPUT_TAG', 'sedumi_1000_primal_dual');
run_random_misc_seq_ico_gap_distribution
```

The random experiment uses the fixed seed `20260623`.
