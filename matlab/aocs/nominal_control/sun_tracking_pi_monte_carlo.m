%% Monte Carlo: MTQ-based sun tracking controller with PI control
%  1000 runs, 2 orbits each

clear; clc; close all;

%% ========= CONSTANTS & MODEL PARAMETERS =========
Re   = 6378.137e3;           % Earth radius [m]
muE  = 3.986004418e14;       % Earth mu [m^3/s^2]
h0   = 510e3;                % nominal altitude [m]
a    = Re + h0;              % semi-major axis [m]
incl = deg2rad(97.4);        % SSO-ish inclination [rad]

n    = sqrt(muE/a^3);        % mean motion [rad/s]
Torb = 2*pi/n;               % orbital period [s]
Tsim = 2*Torb;               % simulate 2 orbits

% Earth rotation
omega_E = [0; 0; 7.2921159e-5];  % [rad/s]

%% Satellite properties
m_sat = 40;                 % [kg]

% Accelerometer location in body frame (relative to COM) [m]
r_acc_body = [0.493; 0; 0];

% MTQ vibration coefficient (simple scalar model)
k_vib = 0.01;                 % dimensionless, ~1% of translational accel as noise

% Bus dimensions
Lx = 0.85; Ly = 0.40; Lz = 0.50;

% Solar arrays (two wings, each 0.85 x 0.40)
A_array = 0.85 * 0.40;
A_arrays_total = 2 * A_array;

% Reference areas
A_drag_ref = 0.85 * 0.40 + A_arrays_total;  % rough drag area
A_srp_ref  = 0.40 * 0.50 + A_arrays_total;  % rough SRP area

% CP offsets (body frame, guesses)
r_cp_drag_body = [ 0.20;  0.00;  0.10 ];
r_cp_srp_body  = [-0.10;  0.00;  0.20 ];

% Coefficients
Cd = 2.2;
Cr = 1.3;

% Solar radiation pressure
P_srp = 4.5e-6;       % [N/m^2]

% Exponential atmosphere model
rho_ref = 1e-12;      % [kg/m^3] at h_ref
h_ref   = 400e3;      % [m]
H_scale = 60e3;       % [m]

%% Inertia (kg·m^2)  (given)
J = [  1.991   -0.04488  -0.03172 ; ...
      -0.04488  2.721     0.03037 ; ...
      -0.03172  0.03037   1.329  ];

%% Orbit frame rotation (perifocal -> ECI)
R_oi = [1 0 0; ...
        0 cos(incl) -sin(incl); ...
        0 sin(incl)  cos(incl)];

%% Sun direction (ECI, fixed for this simplified model)
s_eci_const = [1; 0; 0];
s_eci_const = s_eci_const / norm(s_eci_const);

%% MTQ properties
m_max_axis = 15;                         % max dipole per axis [A·m^2]
m_res_body = (1/sqrt(3)) * [1;1;1];      % residual dipole, |m_res| ≈ 1 A·m^2

%% Controller gains (PD kept as-is, add small integral)
Kp = 5e-4;   % proportional on sun-direction error
Kd = 1e-2;   % derivative on angular rate
Ki = 1e-5;   % small integral gain on sun-direction error

%% Pack parameters into struct
p.muE    = muE;
p.Re     = Re;
p.a      = a;
p.n      = n;
p.R_oi   = R_oi;
p.J      = J;
p.omegaE = omega_E;

p.Cd = Cd; p.Cr = Cr;
p.P_srp = P_srp;
p.A_drag_ref = A_drag_ref;
p.A_srp_ref  = A_srp_ref;
p.r_cp_drag_body = r_cp_drag_body;
p.r_cp_srp_body  = r_cp_srp_body;
p.rho_ref  = rho_ref;
p.h_ref    = h_ref;
p.H_scale  = H_scale;

p.m_max_axis   = m_max_axis;
p.m_res_body   = m_res_body;
p.s_eci_const  = s_eci_const;

p.Kp = Kp;
p.Kd = Kd;
p.Ki = Ki;

% disturbance scale factors
p.gg_scale   = 1.0;   % gravity-gradient
p.drag_scale = 1.0;   % aerodynamic drag
p.srp_scale  = 1.0;   % solar radiation pressure
p.mres_scale = 1.0;   % residual dipole

%% ========= MONTE CARLO SETTINGS =========
N_runs = 1000;

% Storage for summary metrics
final_sun_err = zeros(1, N_runs);
max_sun_err   = zeros(1, N_runs);

% For plotting, store time + sun error + omega for each run
t_runs         = cell(1, N_runs);
sun_err_runs   = cell(1, N_runs);
omega_runs     = cell(1, N_runs);
omega_mag_runs = cell(1, N_runs);
micro_g_runs   = cell(1, N_runs);   % microgravity (|a_non-grav| / g)

rng(1);   % fix seed for reproducibility

%% ========= MONTE CARLO LOOP =========
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

for rIdx = 1:N_runs
    fprintf('Run %d / %d\n', rIdx, N_runs);

    % Initial state: perfect alignment, zero rate, zero integral
    q0 = [1; 0; 0; 0];     % identity quaternion
    w0 = [0; 0; 0];        % zero angular rate
    z0 = [0; 0; 0];        % integral of e_sun
    x0 = [q0; w0; z0];     % 10x1 state

    % per-run disturbance uncertainty
    p_run = p;
    p_run.gg_scale   = 1 + 0.01*randn;
    p_run.drag_scale = 1 + 0.15*randn;
    p_run.srp_scale  = 1 + 0.06*randn;
    p_run.mres_scale = 1 + 0.20*randn;

    % Integrate dynamics for this run
    [t, x] = ode45(@(t,x) sat_dyn(t, x, p_run), [0 Tsim], x0, opts);

    q_hist = x(:,1:4).';
    w_hist = x(:,5:7).';

    Nt = length(t);
    sun_err_deg = zeros(1,Nt);
    micro_g     = zeros(1,Nt);

    for k = 1:Nt
        tk = t(k);
        qk = q_hist(:,k) / norm(q_hist(:,k));   % normalize
        C_bi = quat_to_dcm_bi(qk);

        %% Sun alignment error
        s_b = C_bi * p_run.s_eci_const;
        s_b = s_b / norm(s_b);
        s_des_b = [1; 0; 0];                    % body +x should point to Sun

        cosang = dot(s_b, s_des_b);
        cosang = max(min(cosang,1),-1);
        sun_err_deg(k) = acos(cosang) * 180/pi;

        %% Orbit position & velocity at time tk
        theta = p_run.n * tk;
        r_orb = [p_run.a*cos(theta); p_run.a*sin(theta); 0];
        v_orb = [-p_run.a*p_run.n*sin(theta); p_run.a*p_run.n*cos(theta); 0];

        r_eci = p_run.R_oi * r_orb;
        v_eci = p_run.R_oi * v_orb;

        %% Translational non-grav forces (drag + SRP) at COM
        % Atmosphere density
        r_norm = norm(r_eci);
        h_inst = r_norm - p_run.Re;
        rho    = p_run.rho_ref * exp(-(h_inst - p_run.h_ref)/p_run.H_scale);

        % Relative velocity
        v_atm_eci = cross(p_run.omegaE, r_eci);
        v_rel_eci = v_eci - v_atm_eci;
        v_rel_body = C_bi * v_rel_eci;

        v_mag = norm(v_rel_body);
        Fd_body = [0;0;0];
        if v_mag >= 1e-6
            v_hat   = v_rel_body / v_mag;
            Fd_nom  = -0.5 * rho * v_mag^2 * p_run.Cd * p_run.A_drag_ref * v_hat;
            Fd_body = p_run.drag_scale * Fd_nom;
        end

        % SRP force in body
        s_body = C_bi * p_run.s_eci_const;
        s_body = s_body / norm(s_body);
        F_srp_nom  = -p_run.P_srp * p_run.Cr * p_run.A_srp_ref * s_body;
        F_srp_body = p_run.srp_scale * F_srp_nom;

        % Translational non-grav acceleration at COM in body frame
        a_trans_body = (Fd_body + F_srp_body) / m_sat;   % [m/s^2]

        %% Rotational kinematics at time tk
        w_body = w_hist(:,k);      % angular velocity in body frame

        if k == 1
            if Nt > 1
                dt    = t(2) - t(1);
                alpha = (w_hist(:,2) - w_hist(:,1)) / dt;
            else
                alpha = [0;0;0];
            end
        else
            dt    = t(k) - t(k-1);
            alpha = (w_hist(:,k) - w_hist(:,k-1)) / dt;
        end

        %% Rotational acceleration contributions at accelerometer
        % α × r term
        a_rot_body = cross(alpha, r_acc_body);

        % Centrifugal ω × (ω × r) term
        a_cent_body = cross(w_body, cross(w_body, r_acc_body));

        % Coriolis term is zero here because r_acc_body is fixed in the body
        a_cor_body = [0;0;0];

        %% Simple MTQ vibration model (random accel noise proportional to a_trans)
        vib_scale   = k_vib * norm(a_trans_body);   % [m/s^2]
        a_vib_body  = vib_scale * randn(3,1);

        %% Total specific force measured by accelerometer (in body frame)
        a_meas_body = a_trans_body + a_rot_body + a_cent_body + a_cor_body + a_vib_body;

        %% Microgravity level [g]
        micro_g(k) = norm(a_meas_body) / 9.81;
    end

    % Save metrics
    t_runs{rIdx}         = t;
    sun_err_runs{rIdx}   = sun_err_deg;
    omega_runs{rIdx}     = w_hist;
    omega_mag_runs{rIdx} = sqrt(sum(w_hist.^2, 1));
    micro_g_runs{rIdx}   = micro_g;

    final_sun_err(rIdx) = sun_err_deg(end);
    max_sun_err(rIdx)   = max(sun_err_deg);
end

%% ========= PLOTS =========

% 1) Sun alignment error for all runs
figure; hold on;

for rIdx = 1:N_runs
    plot(t_runs{rIdx}/60, sun_err_runs{rIdx}, ...
        'Color', [0 0 1 0.01]);
end

t_uniform = linspace(0, Tsim, 2000);
sun_grid  = zeros(N_runs, length(t_uniform));

for rIdx = 1:N_runs
    sun_grid(rIdx,:) = interp1(t_runs{rIdx}, sun_err_runs{rIdx}, ...
                               t_uniform, 'linear', 'extrap');
end

sun_mean = mean(sun_grid,1);
sun_max  = max(sun_grid,[],1);

plot(t_uniform/60, sun_mean, 'k', 'LineWidth', 1.5);
plot(t_uniform/60, sun_max, 'r--', 'LineWidth', 1);

xlabel('Time [min]');
ylabel('Sun alignment error [deg]');
title('Sun alignment error (transparent runs, mean, max)');
legend('Runs (transparent)','Mean','Max');
grid on;

% 2) Angular rates per axis
figure;

% omega_x
subplot(3,1,1); hold on;
for rIdx = 1:N_runs
    plot(t_runs{rIdx}/60, omega_runs{rIdx}(1,:), ...
        'Color', [0 0 1 0.01]);
end

omega_grid_x = zeros(N_runs, length(t_uniform));
for rIdx = 1:N_runs
    omega_grid_x(rIdx,:) = interp1(t_runs{rIdx}, omega_runs{rIdx}(1,:), ...
                                   t_uniform, 'linear', 'extrap');
end

plot(t_uniform/60, mean(omega_grid_x,1), 'k', 'LineWidth', 1.5);
plot(t_uniform/60, max(omega_grid_x,[],1), 'r--', 'LineWidth', 1);

ylabel('\omega_x [rad/s]');
title('Angular rates (transparent runs, mean, max)');
grid on;

% omega_y
subplot(3,1,2); hold on;
for rIdx = 1:N_runs
    plot(t_runs{rIdx}/60, omega_runs{rIdx}(2,:), ...
        'Color', [0 0 1 0.01]);
end

omega_grid_y = zeros(N_runs, length(t_uniform));
for rIdx = 1:N_runs
    omega_grid_y(rIdx,:) = interp1(t_runs{rIdx}, omega_runs{rIdx}(2,:), ...
                                   t_uniform, 'linear', 'extrap');
end

plot(t_uniform/60, mean(omega_grid_y,1), 'k', 'LineWidth', 1.5);
plot(t_uniform/60, max(omega_grid_y,[],1), 'r--', 'LineWidth', 1);

ylabel('\omega_y [rad/s]');
grid on;

% omega_z
subplot(3,1,3); hold on;
for rIdx = 1:N_runs
    plot(t_runs{rIdx}/60, omega_runs{rIdx}(3,:), ...
        'Color', [0 0 1 0.01]);
end

omega_grid_z = zeros(N_runs, length(t_uniform));
for rIdx = 1:N_runs
    omega_grid_z(rIdx,:) = interp1(t_runs{rIdx}, omega_runs{rIdx}(3,:), ...
                                   t_uniform, 'linear', 'extrap');
end

plot(t_uniform/60, mean(omega_grid_z,1), 'k', 'LineWidth', 1.5);
plot(t_uniform/60, max(omega_grid_z,[],1), 'r--', 'LineWidth', 1);

ylabel('\omega_z [rad/s]');
xlabel('Time [min]');
grid on;

% 2b) Total angular rate magnitude |omega|
figure; hold on;

h_runs = plot(nan, nan, 'Color', [0 0 1 0.05]);  % dummy for legend

for rIdx = 1:N_runs
    plot(t_runs{rIdx}/60, omega_mag_runs{rIdx}, ...
         'Color', [0 0 1 0.01]);
end

omega_mag_grid = zeros(N_runs, length(t_uniform));
for rIdx = 1:N_runs
    omega_mag_grid(rIdx,:) = interp1(t_runs{rIdx}, omega_mag_runs{rIdx}, ...
                                     t_uniform, 'linear', 'extrap');
end

omega_mag_mean = mean(omega_mag_grid, 1);
omega_mag_max  = max(omega_mag_grid, [], 1);
omega_mag_p25  = prctile(omega_mag_grid, 25, 1);
omega_mag_p75  = prctile(omega_mag_grid, 75, 1);

h_mean = plot(t_uniform/60, omega_mag_mean, 'k', 'LineWidth', 1.5);
h_max  = plot(t_uniform/60, omega_mag_max,  'r--','LineWidth', 1);
h_p25  = plot(t_uniform/60, omega_mag_p25, 'Color',[0 0.6 0],   'LineWidth',1);
h_p75  = plot(t_uniform/60, omega_mag_p75, 'Color',[0.7 0 0.7], 'LineWidth',1);

xlabel('Time [min]');
ylabel('|\omega| [rad/s]');
title('Total angular rate magnitude (transparent runs, mean, max, 25%, 75%)');

leg = legend([h_runs h_mean h_max h_p25 h_p75], ...
             {'Runs','Mean','Max','25%','75%'});
leg.FontSize = 7;
leg.Box = 'on';
leg.NumColumns = 1;

grid on;

% 3) Microgravity (non-grav accel magnitude, in g)
figure; hold on;

h_runs_mg = plot(nan, nan, 'Color', [0 0 1 0.05]);

for rIdx = 1:N_runs
    plot(t_runs{rIdx}/60, micro_g_runs{rIdx}, ...
         'Color', [0 0 1 0.01]);
end

micro_g_grid = zeros(N_runs, length(t_uniform));
for rIdx = 1:N_runs
    micro_g_grid(rIdx,:) = interp1(t_runs{rIdx}, micro_g_runs{rIdx}, ...
                                   t_uniform, 'linear', 'extrap');
end

mg_mean = mean(micro_g_grid,1);
mg_max  = max(micro_g_grid,[],1);

h_mean_mg = plot(t_uniform/60, mg_mean, 'k', 'LineWidth', 1.5);
h_max_mg  = plot(t_uniform/60, mg_max,  'r--','LineWidth', 1);

xlabel('Time [min]');
ylabel('Microgravity [g]');
title('Microgravity level (drag + SRP + rotation + vib)');

leg2 = legend([h_runs_mg h_mean_mg h_max_mg], ...
              {'Runs','Mean','Max'});
leg2.FontSize = 7;
leg2.Box = 'on';
leg2.NumColumns = 1;

grid on;

% 4) Histogram of max sun error during each run
figure;
histogram(max_sun_err, 'BinWidth', 2);
xlabel('Max sun alignment error [deg]');
ylabel('Count');
title('Max sun error during 2 orbits (Monte Carlo)');
grid on;

%% ================= LOCAL FUNCTIONS =================

function dx = sat_dyn(t, x, p)
    % x = [q(1:4); w(5:7); z(8:10)]
    q = x(1:4);
    w = x(5:7);
    z = x(8:10);

    q = q / norm(q);

    C_bi = quat_to_dcm_bi(q);

    theta = p.n * t;
    r_orb = [p.a*cos(theta); p.a*sin(theta); 0];
    v_orb = [-p.a*p.n*sin(theta); p.a*p.n*cos(theta); 0];

    r_eci = p.R_oi * r_orb;
    v_eci = p.R_oi * v_orb;

    B_eci  = earth_dipole_B(r_eci);
    B_body = C_bi * B_eci;

    [tau_gg, tau_drag, tau_srp] = dist_torques(r_eci, v_eci, C_bi, p);
    m_res_eff = p.mres_scale * p.m_res_body;
    tau_res   = cross(m_res_eff, B_body);
    tau_dist  = tau_gg + tau_drag + tau_srp + tau_res;

    s_b = C_bi * p.s_eci_const;
    s_b = s_b / norm(s_b);
    s_des = [1;0;0];

    e_sun  = cross(s_b, s_des);

    % PI controller on sun-direction error
    tau_cmd = -p.Kp * e_sun - p.Kd * w - p.Ki * z;

    if norm(B_body) < 1e-10
        m_cmd = [0;0;0];
    else
        m_cmd = cross(B_body, tau_cmd) / (norm(B_body)^2);
        m_cmd = max(min(m_cmd, p.m_max_axis), -p.m_max_axis);
    end

    tau_mtq = cross(m_cmd, B_body);
    tau_total = tau_dist + tau_mtq;

    wdot = p.J \ (tau_total - cross(w, p.J*w));

    wx = w(1); wy = w(2); wz = w(3);
    Omega = [ 0   -wx  -wy  -wz;
              wx   0    wz  -wy;
              wy  -wz   0    wx;
              wz   wy  -wx   0  ];
    qdot = 0.5 * Omega * q;

    % Integral state dynamics
    zdot = e_sun;

    dx = [qdot; wdot; zdot];
end

function [tau_gg, tau_drag, tau_srp] = dist_torques(r_eci, v_eci, C_bi, p)
    r_norm = norm(r_eci);
    r_hat_i = -r_eci / r_norm;
    r_hat_b = C_bi * r_hat_i;

    tau_gg_nom = 3*p.muE/r_norm^3 * cross(r_hat_b, p.J*r_hat_b);
    tau_gg     = p.gg_scale * tau_gg_nom;

    h_inst = r_norm - p.Re;
    rho    = p.rho_ref * exp(-(h_inst - p.h_ref)/p.H_scale);

    v_atm_eci = cross(p.omegaE, r_eci);
    v_rel_eci = v_eci - v_atm_eci;
    v_rel_body = C_bi * v_rel_eci;

    v_mag = norm(v_rel_body);
    if v_mag < 1e-6
        tau_drag = [0;0;0];
    else
        v_hat = v_rel_body / v_mag;
        Fd_nom = -0.5 * rho * v_mag^2 * p.Cd * p.A_drag_ref * v_hat;
        tau_drag_nom = cross(p.r_cp_drag_body, Fd_nom);
        tau_drag     = p.drag_scale * tau_drag_nom;
    end

    s_body = C_bi * p.s_eci_const;
    s_body = s_body / norm(s_body);
    F_srp_nom = -p.P_srp * p.Cr * p.A_srp_ref * s_body;
    tau_srp_nom = cross(p.r_cp_srp_body, F_srp_nom);
    tau_srp     = p.srp_scale * tau_srp_nom;
end

function B = earth_dipole_B(r_eci)
    M   = 7.96e15;  % [T·m^3]
    r   = norm(r_eci);
    er  = r_eci / r;
    m_hat = [0;0;1];

    coef = M / r^3;
    mdot = dot(m_hat, er);
    B = coef * (3*mdot*er - m_hat);
end

function C = quat_to_dcm_bi(q)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    C = [ 1-2*(q2^2+q3^2),   2*(q1*q2+q0*q3),   2*(q1*q3-q0*q2);
          2*(q1*q2-q0*q3),   1-2*(q1^2+q3^2),   2*(q2*q3+q0*q1);
          2*(q1*q3+q0*q2),   2*(q2*q3-q0*q1),   1-2*(q1^2+q2^2)];
end
