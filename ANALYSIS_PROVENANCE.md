# Analysis provenance

The report is a historical group design document. Public-release analysis caveats are kept in this companion document so the report retains its original formatting. Organisation and compilation fixes do not establish the validity of its engineering conclusions.

## Released models

- The report's final mission altitude is **520 km**. The supplied MATLAB studies use **510 km**. They are retained as development studies; silently changing altitude would alter the analysis rather than restore the original evidence.
- The PD Monte Carlo script defaults to **10 runs**, not 1,000. It draws independent Gaussian per-run scale factors for gravity-gradient, drag, SRP and residual dipole torques, and computes an offset-sensor acceleration with a simple vibration term.
- The integral-augmented variant defaults to **1,000 runs** and uses the same style of Gaussian scaling.
- The stochastic variant defaults to **100 runs**, uses fixed-step RK4, Ornstein-Uhlenbeck correlated environmental disturbances and lognormal static uncertainties, and plots translational drag/SRP acceleration at the centre of mass.
- These distributions do not match the report's original description of triangular disturbance distributions. That historical claim has not been retroactively attributed to the released scripts.
- The saved microgravity figure's title and legend match the stochastic variant's COM metric. This is a candidate source, not proof of the exact generating revision or settings. It does not demonstrate offset-payload microgravity compliance.
- The actuator comparison script uses a 15 A m^2 per-axis MTQ setting; the nominal controller studies use 10 A m^2. Other trade studies contain additional assumptions and hardware variants.

## Repeating the checks

In a fresh MATLAB session, from the repository root:

```matlab
addpath('matlab');
manifest = run_analysis('smoke');
```

This runs all 12 original analyses, reduces only `N_runs` to 2 in temporary execution copies, seeds the RNG and saves figures, workspaces and a JSON manifest under `outputs/matlab/smoke`. A smoke run establishes that code executes; two samples cannot establish Monte Carlo convergence or requirement compliance.

```matlab
manifest = run_analysis('full');
```

Full mode uses the released scripts' unmodified run counts. The scripts clear the base workspace and close figures. Both modes record the original and executed source hashes, the exact executed `.m` source, the initial RNG state, installed MATLAB products, elapsed time and errors. Explicit per-script RNG calls remain effective and are preserved in the recorded source. Each result directory also contains the saved final workspace. Generated files are ignored by Git.

MATLAB R2025b passed all 12 smoke executions on 30 October 2026. Full Monte Carlo studies have not been rerun for this release. The minimum MATLAB release and minimal toolbox set are not yet established; the manifest records the environment actually tested.

Re-establishing final-design compliance requires a separately reviewed model configuration at 520 km, agreed actuator/inertia inputs, a defined payload acceleration metric and an adequately sized convergence study. The archived figures have been preserved rather than presented as newly verified results.
