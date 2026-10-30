%% Monte Carlo: MTQ sun tracking (from-scratch, RK4 + colored disturbances)
clear; clc; close all;

%% ========= CONSTANTS & MODEL PARAMETERS =========
Re   = 6378.137e3;                % Earth radius [m]
muE  = 3.986004418e14;            % Earth mu [m^3/s^2]
h0   = 510e3;                     % nominal altitude [m]
a    = Re + h0;                   % semi-major axis [m]
incl = deg2rad(97.4);             % SSO-ish inclination [rad]

n    = sqrt(muE/a^3);             % mean motion [rad/s]
Torb = 2*pi/n;                    % orbital period [s]
Tsim = 2*Torb;                    % simulate 2 orbits

% Fixed-step integrator
dt   = 0.5;                       % [s] (band-limited env -> this is plenty)
t    = 0:dt:Tsim; Nt = numel(t);

omega_E = [0; 0; 7.2921159e-5];   % Earth rotation [rad/s]

%% Spacecraft geometry/mass
m_sat = 40;                       % [kg]
J = [  1.991   -0.04488  -0.03172 ; ...
      -0.04488  2.721     0.03037 ; ...
      -0.03172  0.03037   1.329  ];  % kg m^2

% Surfaces (simple refs)
A_array = 0.85*0.40;                        % each wing
A_arrays_total = 2*A_array;
A_drag_ref = 0.85*0.40 + A_arrays_total;
A_srp_ref  = 0.40*0.50 + A_arrays_total;

% CP offsets (body frame, nominal guesses)
r_cp_drag_body = [ 0.20;  0.00;  0.10 ];
r_cp_srp_body  = [-0.10;  0.00;  0.20 ];

% Coefficients (nominal means)
Cd_mu = 2.2;
Cr_mu = 1.3;

% SRP
P_srp = 4.5e-6;                   % [N/m^2] (absorbing baseline; reflectivity in Cr)

% Exponential atmosphere baseline
rho_ref = 1e-12;    % [kg/m^3] at h_ref
h_ref   = 400e3;    % [m]
H_scale = 60e3;     % [m]

% Controller
Kp = 6e-3;  Kd = 2e-2;

% MTQ actuator
m_max_axis = 10;                               % [A m^2] per axis
m_res_dir  = (1/sqrt(3))*[1;1;1];              % fixed direction; magnitude per run

% Orbit frame (perifocal->ECI)
R_oi = [1 0 0; 0 cos(incl) -sin(incl); 0 sin(incl) cos(incl)];

% Sun (fixed direction in this simplified model)
s_eci = [1;0;0];

% Accelerometer (for micro-g metric)
r_acc_body = [0.493; 0; 0];
k_vib = 0.01;                                    % accel noise scale

%% ========= STOCHASTIC MODELS =========
% OU parameters (tune envelope)
tau_rho = 600;  sig_rho = 0.30;    % density multiplier 1Â±30%, ~10 min corr
tau_w   = 300;  sig_w   = 150;     % along-track wind [m/s], ~5 min corr
tau_srp = 3600; sig_srp = 0.03;    % SRP slow factor Â±3%
tau_mj  = 1000; sig_mj  = 0.05;    % residual dipole jitter Â±5%

% Static per-run uncertainties (relative)
Cd_rel  = 0.20;     % ~20% geometric std
Cr_rel  = 0.05;     % ~5%
mres_rel= 0.50;     % wide spread around 1 A m^2
cp_jit  = 0.02;     % Â±2 cm CP perturbation

%% ========= MONTE CARLO =========
rng(1);
N_runs = 100;       % set to 1000 if you want; dt=0.5 keeps it tractable

% Output logs
final_sun_err = zeros(1,N_runs);
max_sun_err   = zeros(1,N_runs);
t_runs         = cell(1,N_runs);
sun_err_runs   = cell(1,N_runs);
omega_runs     = cell(1,N_runs);
omega_mag_runs = cell(1,N_runs);
micro_g_runs   = cell(1,N_runs);

for r = 1:N_runs
    fprintf('Run %d / %d\n', r, N_runs);
    % ----- Per-run STATIC draws -----
    Cd = logn_from_mu_rel(Cd_mu, Cd_rel);
    Cr = logn_from_mu_rel(Cr_mu, Cr_rel);

    mres_mag   = logn_from_mu_rel(1.0, mres_rel);
    m_res_body = mres_mag * m_res_dir;

    r_cp_drag_b = r_cp_drag_body + cp_jit*(2*rand(3,1)-1);
    r_cp_srp_b  = r_cp_srp_body  + cp_jit*(2*rand(3,1)-1);

    % ----- Colored time series (all on the sim grid) -----
    rho_mult = 1 + ou_series(Nt, dt, tau_rho, sig_rho, 0);
    wind_al  =     ou_series(Nt, dt, tau_w,   sig_w,   0);
    srp_mult = 1 + ou_series(Nt, dt, tau_srp, sig_srp, 0);
    mres_jit = 1 + ou_series(Nt, dt, tau_mj,  sig_mj,  0);

    % ----- State init -----
    q = [1;0;0;0];    % body->ECI quaternion
    w = [0;0;0];

    % Logs
    sun_err = zeros(1,Nt);
    omega_mag = zeros(1,Nt);
    micro_g   = zeros(1,Nt);

    % ----- Fixed vectors for orbit propagation -----
    % Simple circular orbit param form, evaluated as needed inside dyn

    % ----- Main propagation (RK4) -----
    for k = 1:Nt
        tk = t(k);

        % Metrics (computed after we have state)
        C_bi = quat_to_dcm_bi(q);
        s_b  = C_bi*s_eci; s_b = s_b/norm(s_b);
        sun_err(k) = acos(max(-1,min(1, dot(s_b,[1;0;0]))))*180/pi;
        omega_mag(k) = norm(w);

        % Micro-g proxy computed after force calculation (below)

        % Precompute env scalars at this step
        rho_m = rho_mult(k);
        w_al  = wind_al(k);
        srp_m = srp_mult(k);
        mj    = mres_jit(k);

        % One RK4 step for rotational dynamics only (translation is analytic)
        f = @(q_,w_) rhs_rot(tk, q_, w_, ...
                             a, n, R_oi, ...
                             J, muE, Re, ...
                             rho_ref, h_ref, H_scale, ...
                             omega_E, ...
                             Cd, Cr, P_srp, A_drag_ref, A_srp_ref, ...
                             r_cp_drag_b, r_cp_srp_b, ...
                             s_eci, m_max_axis, ...
                             m_res_body, mj, ...
                             rho_m, w_al, srp_m);

        [q, w, a_trans_body] = rk4_step_rot(q, w, dt, f);

        % Add accel-derived vib noise and form micro-g metric
        a_vib = k_vib*norm(a_trans_body)*randn(3,1);
        micro_g(k) = norm(a_trans_body + a_vib)/9.81;

        % Normalize quaternion (cheap)
        q = q / norm(q);
    end

    % Store
    t_runs{r}         = t;
    sun_err_runs{r}   = sun_err;
    omega_runs{r}     = zeros(3,Nt);   % if you want per-axis, adapt rk4 to return
    omega_mag_runs{r} = omega_mag;
    micro_g_runs{r}   = micro_g;

    final_sun_err(r) = sun_err(end);
    max_sun_err(r)   = max(sun_err);
end

%% ========= PLOTS =========
% 1) Sun alignment error
figure; hold on;
for r = 1:N_runs
    plot(t_runs{r}/60, sun_err_runs{r}, 'Color', [0 0 1 0.02]);
end
[xm, xX, xu] = env_stats(t_runs, sun_err_runs, t, 2000);
plot(xu/60, xm, 'k', 'LineWidth', 1.5);
plot(xu/60, xX, 'r--', 'LineWidth', 1);
xlabel('Time [min]'); ylabel('Sun alignment error [deg]');
title('Sun alignment error (transparent runs, mean, max)'); grid on;
legend('Runs','Mean','Max');

% 2) |omega|
figure; hold on;
for r = 1:N_runs
    plot(t_runs{r}/60, omega_mag_runs{r}, 'Color', [0 0 1 0.02]);
end
[wm, wX, wu] = env_stats(t_runs, omega_mag_runs, t, 2000);
plot(wu/60, wm, 'k', 'LineWidth', 1.5);
plot(wu/60, wX, 'r--', 'LineWidth', 1);
xlabel('Time [min]'); ylabel('|\omega| [rad/s]');
title('Angular rate magnitude'); grid on; legend('Runs','Mean','Max');

% 3) Micro-g
figure; hold on;
for r = 1:N_runs
    plot(t_runs{r}/60, micro_g_runs{r}, 'Color', [0 0 1 0.02]);
end
[mm, mX, mu] = env_stats(t_runs, micro_g_runs, t, 2000);
plot(mu/60, mm, 'k', 'LineWidth', 1.5);
plot(mu/60, mX, 'r--', 'LineWidth', 1);
xlabel('Time [min]'); ylabel('Microgravity [g]');
title('Non-gravitational acceleration level'); grid on; legend('Runs','Mean','Max');

% 4) Histogram of max sun error per run
figure; histogram(max_sun_err, 'BinWidth', 2);
xlabel('Max sun alignment error [deg]'); ylabel('Count');
title('Max sun error during 2 orbits'); grid on;


%% ==================== LOCAL FUNCTIONS ====================
function [q_next, w_next, a_trans_body] = rk4_step_rot(q, w, h, f)
    % Single RK4 step; f returns [qdot, wdot, a_trans_body] at given (q,w)
    [k1q, k1w, ~] = f(q, w);
    [k2q, k2w, ~] = f(q + 0.5*h*k1q, w + 0.5*h*k1w);
    [k3q, k3w, ~] = f(q + 0.5*h*k2q, w + 0.5*h*k2w);
    [k4q, k4w, a4] = f(q +     h*k3q, w +     h*k3w);

    q_next = q + (h/6)*(k1q + 2*k2q + 2*k3q + k4q);
    w_next = w + (h/6)*(k1w + 2*k2w + 2*k3w + k4w);

    % Return last a_trans_body (good enough for metric)
    a_trans_body = a4;
end

function [qdot, wdot, a_trans_body] = rhs_rot(tk, q, w, ...
    a, n, R_oi, J, muE, Re, ...
    rho_ref, h_ref, H_scale, ...
    omega_E, ...
    Cd, Cr, P_srp, A_drag_ref, A_srp_ref, ...
    r_cp_drag_b, r_cp_srp_b, ...
    s_eci, m_max_axis, ...
    m_res_body, mj, ...
    rho_mult_k, wind_al_k, srp_mult_k)

    % Orientation stuff
    C_bi = quat_to_dcm_bi(q);

    % Orbit position/velocity (circular, perifocal -> ECI)
    theta = n*tk;
    r_orb = [a*cos(theta); a*sin(theta); 0];
    v_orb = [-a*n*sin(theta); a*n*cos(theta); 0];
    r_eci = R_oi * r_orb;
    v_eci = R_oi * v_orb;

    % Magnetic field (dipole approx) in body
    B_eci  = earth_dipole_B(r_eci);
    B_body = C_bi * B_eci;

    % Gravity-gradient torque (deterministic)
    r_norm = norm(r_eci);
    r_hat_b = C_bi * (-r_eci/r_norm);
    tau_gg  = 3*muE/r_norm^3 * cross(r_hat_b, J*r_hat_b);

    % Atmosphere density (baseline * colored multiplier)
    h_inst = r_norm - Re;
    rho0   = rho_ref * exp(-(h_inst - h_ref)/H_scale);
    rho    = rho0 * rho_mult_k;

    % Relative velocity incl. co-rotating atm + colored along-track wind
    v_atm_eci = cross(omega_E, r_eci);
    v_rel_eci = v_eci - v_atm_eci;

    vhat_orb  = v_eci / max(1e-9, norm(v_eci));
    v_rel_eci = v_rel_eci - wind_al_k * vhat_orb;

    v_rel_b = C_bi * v_rel_eci;
    vmag    = norm(v_rel_b);

    if vmag < 1e-6
        Fd_b = [0;0;0];
        tau_drag = [0;0;0];
    else
        v_hat = v_rel_b / vmag;
        Fd_b  = -0.5 * rho * vmag^2 * Cd * A_drag_ref * v_hat;
        tau_drag = cross(r_cp_drag_b, Fd_b);
    end

    % SRP torque with slow colored multiplier
    s_b = C_bi * s_eci; s_b = s_b / norm(s_b);
    F_srp_b = -P_srp * Cr * A_srp_ref * srp_mult_k * s_b;
    tau_srp = cross(r_cp_srp_b, F_srp_b);

    % Residual magnetic torque (with small jitter)
    tau_res = cross(mj * m_res_body, B_body);

    % Total disturbance
    tau_dist = tau_gg + tau_drag + tau_srp + tau_res;

    % Controller (sun tracking)
    e_sun  = cross(s_b, [1;0;0]);
    tau_cmd = -6e-3 * e_sun - 2e-2 * w;

    % MTQ actuation with saturation
    if norm(B_body) < 1e-12
        m_cmd = [0;0;0];
    else
        m_cmd = cross(B_body, tau_cmd) / (norm(B_body)^2);
        m_cmd = max(min(m_cmd, m_max_axis), -m_max_axis);
    end
    tau_mtq = cross(m_cmd, B_body);

    tau_total = tau_dist + tau_mtq;

    % Rigid-body dynamics
    wdot = J \ (tau_total - cross(w, J*w));

    % Quaternion kinematics
    wx=w(1); wy=w(2); wz=w(3);
    Omega = [ 0   -wx  -wy  -wz;
              wx   0    wz  -wy;
              wy  -wz   0    wx;
              wz   wy  -wx   0  ];
    qdot = 0.5 * Omega * q;

    % Return COM translational non-grav accel for micro-g metric
    a_trans_body = (Fd_b + F_srp_b) / 40;  % m_sat is 40 kg; inline for speed
end

function x = ou_series(N, dt, tau, sigma, x0)
    phi = exp(-dt/tau);
    q   = sigma*sqrt(1 - phi^2);
    x   = zeros(N,1); if nargin<5, x0=0; end; x(1)=x0;
    for k = 1:N-1, x(k+1) = phi*x(k) + q*randn; end
end

function m = logn_from_mu_rel(mu, rel)
    % Draw lognormal with mean 'mu' and relative std 'rel'
    s2 = log(1 + rel^2); s = sqrt(s2); mlog = log(mu) - 0.5*s2;
    m = exp(mlog + s*randn);
end

function C = quat_to_dcm_bi(q)
    q = q(:)/norm(q); q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    C = [ 1-2*(q2^2+q3^2), 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2);
          2*(q1*q2-q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3+q0*q1);
          2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), 1-2*(q1^2+q2^2)];
end

function B = earth_dipole_B(r_eci)
    M = 7.96e15; r = norm(r_eci); er = r_eci/r; m_hat=[0;0;1];
    coef = M/r^3; mdot = dot(m_hat, er);
    B = coef*(3*mdot*er - m_hat);
end

function [m, Mx, xgrid] = env_stats(t_runs, series_runs, xgrid, Ngrid)
    if nargin<4, Ngrid = numel(xgrid); end
    xgrid = linspace(xgrid(1), xgrid(end), Ngrid);
    Nr = numel(series_runs);
    G = zeros(Nr, Ngrid);
    for i=1:Nr
        G(i,:) = interp1(t_runs{i}, series_runs{i}, xgrid, 'linear','extrap');
    end
    m  = mean(G,1);
    Mx = max(G,[],1);
end
