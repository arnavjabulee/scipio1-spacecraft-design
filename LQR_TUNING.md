# LQR tuning and comparison with PD

The tuned periodic LQR with linear Kalman estimation improves Sun pointing over the original PD gains in the matched nominal simulation. This comparison concerns body +X alignment with the Sun. It does not establish payload microgravity performance.

## Changes

The original LQR example used a double-integrator attitude model, equal attitude weights and an input scale of 5 A m^2. Tuning explored roll weights, rate weights and input penalties on training noise seeds 1:3. Softer roll penalties improved Sun pointing but could exceed the existing 10-degree full-attitude range check. The selected design retains the original reference attitude and that check.

The scheduled state model now includes the environmental-torque attitude Jacobian about the known reference orbit, calculated by central differences with a 1e-5 radian perturbation. An augmented matrix exponential discretizes the frozen state, magnetic-input and affine-disturbance matrices. The reference schedule is computed before simulation; estimation remains a linear Kalman filter with a fixed measurement matrix.

Selected cost parameters are:

- Attitude scale: 1 degree, with roll weight 0.125 and transverse weights 1.
- Angular-rate scale: 0.05 degrees/s on all axes.
- Dipole scale: 0.25 A m^2, giving input weight matrix R = 16 I.
- Nominal periodic disturbance feedforward remains enabled.

The nonlinear plant, physical initial attitude/rates, 510 km orbit, two-orbit duration, approximately 2 s sampling, sensor noise and 10 A m^2 per-axis actuator limits are unchanged. The filter noise assumptions are unchanged.

## Held-out results

MATLAB R2025b Update 3, 31 October 2026. Seven paired trials use noise seeds 4:10, which were excluded from tuning. Metrics below are arithmetic means of each trial's second-orbit metrics; RMS values are time weighted, and peak means the mean of trial maxima.

| Metric | Original PD gains | Tuned LQR/KF | Change |
| --- | ---: | ---: | ---: |
| RMS Sun error | 1.0235 degrees | 0.9037 degrees | 11.7% lower |
| Peak Sun error | 3.0981 degrees | 2.4933 degrees | 19.5% lower |
| RMS dipole vector magnitude | 0.12866 A m^2 | 0.13176 A m^2 | 2.4% higher |
| Saturated samples | 0% | 0% | None |

Every paired seed improved RMS and peak Sun error, used at most 110% of PD's RMS dipole effort and stayed below 10 degrees full-attitude error for the LQR. The largest LQR full-attitude error was 9.958 degrees across these trials, so there is little margin to the local-model check.

The PD comparison reads Kp = 5e-4 and Kd = 1e-2 from the original source, and uses its cross-product Sun error and magnetic dipole allocation. Both controllers use the common nonlinear quaternion plant, integration, sampling, physical initial state, noise draws and actuator clipping. PD receives raw noisy full-attitude/gyro measurements; LQR receives the linear filter estimate. This evaluates the complete controller/estimator designs. It does not isolate the contribution of LQR from filtering or nominal disturbance preview. The original PD script itself is unchanged; its original continuous, truth-feedback execution is a different experiment.

PD does not constrain roll about the Sun axis. LQR retains a softly weighted full-attitude reference. The stronger Sun-pointing result therefore trades roll precision against transverse pointing. Model uncertainty, sensor biases/outages, other initial conditions and long-duration operation have not been covered by this seven-seed check. The known deterministic disturbance preview remains an assumption.

## Reproduce

From the repository root in a fresh MATLAB session:

```matlab
addpath('scripts');
validation = validate_lqr_kalman;
comparison = benchmark_pd_lqr; % held-out seeds 4:10
```

The benchmark asserts per-seed improvement in both pointing metrics, the effort bound and the existing local-attitude bound. It records executed sources, paired initial states, per-seed metrics and traces under ignored `outputs/lqr_kalman/benchmark/`. Validation records its results and nominal plots under ignored `outputs/lqr_kalman/validation/`.

Validation passed nominal finite-state, actuator, quaternion and covariance checks; filter error was lower than raw measurement error on every attitude/rate axis. Unforced last-orbit RMS Sun error was 0.00054 degrees, reference equilibrium was preserved, and halving the timestep changed the noiseless disturbed last-orbit RMS by 0.00133 degrees.
