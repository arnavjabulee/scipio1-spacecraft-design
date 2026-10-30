# SCIPIO-1 spacecraft design and analysis

Space Forge and University of Bristol AVDASI 2025 Space Group 4 (SS4 SerenSpace) final engineering design report and supporting MATLAB analyses for AOCS. The report covers mission analysis, configuration, structures, propulsion, attitude and orbit control (AOCS), communications, power and thermal design. Original author credits remain in the report; student emails are redacted, and the group photograph has been restored at the project owner's request.

## Contents

- `latex/main.tex`: report entry point; compile from the `latex` directory.
- `latex/figures/`: report figures and drawings, named by their content.
- `latex/references.bib`: merged bibliography with unique citation keys.
- `matlab/aocs/`: analysis organised around the report's AOCS topics.
- `matlab/communications/`: UHF downlink opportunity plot.

## MATLAB analyses

Open a script and run it in MATLAB. Scripts contain their own inputs and helper functions; numerical models were preserved during organisation. Several scripts clear the workspace and close figures. All 12 scripts were checked and executed on MATLAB R2025b. MATLAB dependency analysis identifies MATLAB as the only required product on that release; compatibility with older releases has not been established.

| Report topic | Script under `matlab/` | Purpose |
| --- | --- | --- |
| AOCS actuator selection and sizing | `aocs/actuator_sizing/actuator_torque_authority.m` | Reaction-wheel and magnetorquer torque authority over an orbit |
| AOCS actuator selection and sizing | `aocs/actuator_sizing/disturbance_torque_comparison.m` | Disturbance torques and available directional MTQ torque |
| Microgravity trade | `aocs/microgravity/magnetorquer_microgravity_trade.m` | Dipole/bandwidth trade and per-axis magnetic controllability |
| Microgravity trade | `aocs/microgravity/reaction_wheel_microgravity_estimate.m` | Reaction-wheel disturbance estimate at the payload |
| Control architecture | `aocs/nominal_control/sun_tracking_pd_simulation.m` | Single closed-loop PD sun-tracking simulation |
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

The runner saves figures, workspaces, executed source, hashes and environment information in `outputs/matlab/`. All 12 analyses passed the smoke check on MATLAB R2025b. See `ANALYSIS_PROVENANCE.md` for the 510/520 km distinction and the limits of the archived figures. The original report is a historical design document; the companion analysis notes identify the metric mismatch in the saved Monte Carlo plot.

## Building the report

Use TeX Live or MiKTeX with the packages declared in `latex/main.tex`, plus BibTeX. From the repository root, in PowerShell:

```powershell
./scripts/build_report.ps1
# MiKTeX only: permit missing package installation when necessary
./scripts/build_report.ps1 -InstallMissingPackages
```

The script resolves the bibliography and cross-references, and writes `output/pdf/scipio1_design_report.pdf`. Build files remain in `latex/build`. Shell escape is no longer required. The report was built locally. The report retains its original formatting, including Tables 9 and 10 side by side. Page flow from Thermal through Summary uses compact spacing and flexible table placement to avoid stretched gaps. The group photo is included from `latex/figures/group_photo.png` in its original full-width position. On other systems, run `pdflatex -output-directory=build main.tex`, `bibtex build/main`, then two more pdflatex passes from `latex` (create `build` first).

## Licence and attribution

Contributor permission was confirmed for publication. Project-owned code is [MIT licensed](LICENSE-CODE); original report content and project documentation use [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). See [LICENSE](LICENSE) for scope, attribution and third-party exclusions. The original author credits and upstream IEEEtran notices remain.
