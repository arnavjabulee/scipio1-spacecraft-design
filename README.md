# SCIPIO-1 spacecraft design and analysis

Space Forge and University of Bristol AVDASI 2025 Space Group 4 (SS4 SerenSpace) final engineering design report and supporting MATLAB analyses for AOCS. The report covers mission analysis, configuration, structures, propulsion, attitude and orbit control (AOCS), communications, power and thermal design. Original author credits remain in the report; student emails are redacted, and the group photograph has been restored at the project owner's request.

## Contents

- `latex/main.tex`: report entry point; compile from the `latex` directory.
- `latex/figures/`: report figures and drawings, named by their content.
- `latex/references.bib`: merged bibliography with unique citation keys.
- `matlab/aocs/`: analysis organised around the report's AOCS topics.
- `matlab/communications/`: UHF downlink opportunity plot.

## MATLAB analyses

Open a script and run it in MATLAB. Scripts contain their own inputs and helper functions; numerical models were preserved during organisation. Several scripts clear the workspace and close figures. The 12 original analyses and the added LQR/Kalman example were checked and executed on MATLAB R2025b. These examples use base MATLAB; compatibility with older releases has not been established.

| Report topic | Script under `matlab/` | Purpose |
| --- | --- | --- |
| AOCS actuator selection and sizing | `aocs/actuator_sizing/actuator_torque_authority.m` | Reaction-wheel and magnetorquer torque authority over an orbit |
| AOCS actuator selection and sizing | `aocs/actuator_sizing/disturbance_torque_comparison.m` | Disturbance torques and available directional MTQ torque |
| Microgravity trade | `aocs/microgravity/magnetorquer_microgravity_trade.m` | Dipole/bandwidth trade and per-axis magnetic controllability |
| Microgravity trade | `aocs/microgravity/reaction_wheel_microgravity_estimate.m` | Reaction-wheel disturbance estimate at the payload |
| Control architecture | `aocs/nominal_control/sun_tracking_pd_simulation.m` | Single closed-loop PD sun-tracking simulation |
| LQR and state estimation | `aocs/nominal_control/sun_tracking_lqr_kalman_simulation.m` | Periodic magnetic LQR with a linear Kalman filter and nonlinear truth dynamics |
| Nominal operations verification | `aocs/nominal_control/sun_tracking_pd_monte_carlo.m` | PD Monte Carlo study; 10 runs by default |
| Alternative controller | `aocs/nominal_control/sun_tracking_pi_monte_carlo.m` | Integral-augmented controller variant; 1,000 runs by default |
| Alternative disturbance model | `aocs/nominal_control/sun_tracking_stochastic_monte_carlo.m` | RK4 with correlated stochastic disturbances; 100 runs by default |
| LEOP detumbling | `aocs/detumbling/mtq_vs_reaction_wheel_detumbling.m` | MTQ-only versus wheel-only detumbling from 2 degrees/s |
| Operational modes / power budget | `aocs/duty_cycles/mission_mode_duty_cycles.m` | Duty estimates across mission modes |
| Operational modes / power budget | `aocs/duty_cycles/nominal_actuator_duty_cycles.m` | MTQ-first operation with wheel backup and power estimates |
| Communications ground stations | `communications/uhf_downlink_opportunities.m` | Scatter plot of embedded inter-pass interval data |


## Running the analyses

In a fresh MATLAB session, from the repository root:

```matlab
addpath('matlab');
run_analysis('smoke'); % all scripts, two Monte Carlo runs per study
% run_analysis('full'); % original run counts; substantially slower
```

The runner discovers all 13 analyses and saves figures, workspaces, executed source, hashes and environment information in `outputs/matlab/`. The 12 original analyses passed the smoke check on MATLAB R2025b; the new LQR/Kalman example has a separate nonlinear validation suite. See `ANALYSIS_PROVENANCE.md` for the 510/520 km distinction and the limits of the archived figures.

## LQR with linear Kalman estimation

```matlab
run('matlab/aocs/nominal_control/sun_tracking_lqr_kalman_simulation.m');
% Optional validation, run from a fresh session:
addpath('scripts');
report = validate_lqr_kalman;
```

This additional example retains the single-run PD model's 510 km orbit, inertia, gravity-gradient/drag/SRP/residual-dipole disturbances and 10 A m^2 per-axis magnetic dipole limit. It uses six small-angle attitude/rate states, an approximately two-second sample time, and noisy synthetic star-tracker attitude and gyro-rate measurements. A full sun-aligned attitude reference fixes the otherwise free roll angle. All sensors are continuously available; bias and outages are not modeled.

The magnetic input matrix changes along the orbit. A periodic Riccati recursion designs the LQR gains across that orbit, respecting the actuator's instantaneous rank-two torque geometry. A linear Kalman prediction/update uses these state-space matrices and the actual saturated dipole command. The nonlinear quaternion plant retains the changing field, disturbances and gyroscopic dynamics. Riccati and filter steps are implemented directly, without Control System Toolbox functions or an EKF.

The default also enables model-based periodic disturbance feedforward. Set `use_disturbance_feedforward = false` to examine pure LQR regulation without that preview. This feedforward assumes the nominal deterministic orbit/disturbance model is known; it does not use truth-state feedback. Weights, sensor covariances, sample time and initial errors are editable near the top of the file.

Plots show pointing error, true/estimated rates, disturbance/control torques, dipole limits and attitude-estimation error. The `metrics` structure records numerical diagnostics. This is a local-model demonstration, not flight gain tuning or final mission verification: magnetic underactuation leaves a periodic pointing error of several degrees in the default disturbed case. The linear stability diagnostic applies to the unsaturated reference model, and its guarantee does not extend to large rotations or saturation.

`validate_lqr_kalman` checks disturbed operation, filtering improvement over raw measurements, unforced convergence, reference equilibrium and step-size sensitivity. Its generated execution copies, metrics and nominal plots are saved under ignored `outputs/lqr_kalman/validation/`.

## Building the report

Use TeX Live or MiKTeX with the packages declared in `latex/main.tex`, plus BibTeX. From the repository root, in PowerShell:

```powershell
./scripts/build_report.ps1
# MiKTeX only: permit missing package installation when necessary
./scripts/build_report.ps1 -InstallMissingPackages
```
