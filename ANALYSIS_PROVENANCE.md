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

This discovers the 12 original analyses plus the added LQR/Kalman example, reduces only `N_runs` to 2 in temporary execution copies, seeds the RNG and saves figures, workspaces and a JSON manifest under `outputs/matlab/smoke`. A smoke run establishes that code executes; two samples cannot establish Monte Carlo convergence or requirement compliance.

```matlab
manifest = run_analysis('full');
```

Full mode uses the released scripts' unmodified run counts. The scripts clear the base workspace and close figures. Both modes record the original and executed source hashes, the exact executed `.m` source, the initial RNG state, installed MATLAB products, elapsed time and errors. Explicit per-script RNG calls remain effective and are preserved in the recorded source. Each result directory also contains the saved final workspace. Generated files are ignored by Git.

MATLAB R2025b passed all 12 smoke executions on 30 October 2026. Full Monte Carlo studies have not been rerun for this release. The minimum MATLAB release and minimal toolbox set are not yet established; the manifest records the environment actually tested.

Re-establishing final-design compliance requires a separately reviewed model configuration at 520 km, agreed actuator/inertia inputs, a defined payload acceleration metric and an adequately sized convergence study. The archived figures have been preserved rather than presented as newly verified results.

## Added LQR/Kalman example

`sun_tracking_lqr_kalman_simulation.m` is a new development example, not a recovered source for the report's archived figures. It retains the single-run PD plant's parameters, uses a consistent scalar-first quaternion/inertial-to-body DCM convention, and starts with a small deliberate attitude/rate offset. Historical scripts were not changed to match it.

Its feedback design is a six-state small-angle periodic magnetic LQR; estimation is a linear Kalman filter using synthetic full three-axis attitude and gyro measurements. A sun vector alone would not supply these full-attitude measurements. The nonlinear truth plant, input saturation and the linear estimator are separate. Optional nominal disturbance preview is enabled by default and assumes a known deterministic reference-orbit torque model. The filter's modeled covariance is not proof of statistical consistency under nonlinear model mismatch.

The tuned example includes reference-orbit environmental-torque attitude sensitivity in its scheduled linear model. Its roll attitude weight is 0.125 times the transverse attitude weights, and its input scale is 0.25 A m^2. The same nonlinear truth plant, initial attitude, synthetic sensor noise and actuator limits are retained. The linear Kalman filter is scheduled offline about the reference rather than relinearized about estimated states.

On 31 October 2026, seven held-out paired noise seeds (4:10) gave second-orbit mean RMS Sun errors of 0.904 degrees for the tuned LQR/KF and 1.024 degrees for the original PD gains, with 2.4% more RMS dipole effort. Each seed improved both RMS and peak Sun error. `scripts/benchmark_pd_lqr.m` reproduces the comparison. The original PD source is preserved; the comparison runs its gains and allocation on the common sampled nonlinear plant with raw noisy feedback.

A stable unsaturated linear periodic model does not establish nonlinear stability, saturated stability, pointing compliance or microgravity compliance. This file does not calculate the offset-payload microgravity metric. `scripts/validate_lqr_kalman.m` passed disturbed operation, filtering improvement, unforced convergence, equilibrium and timestep sensitivity checks. Gains and noise assumptions still require mission-specific validation.
