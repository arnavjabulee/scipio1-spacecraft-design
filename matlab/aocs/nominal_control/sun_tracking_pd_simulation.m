%% SSO satellite, MTQ-based sun tracking controller simulation
clear; clc; close all;

%% CONSTANTS
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

%% SATELLITE PROPERTIES
m_sat = 40;                 % mass [kg]

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

%% INERTIA (kg·m^2)  (given)
J = [  1.991   -0.04488  -0.03172 ; ...
      -0.04488  2.721     0.03037 ; ...
      -0.03172  0.03037   1.329  ];

%% ORBIT FRAME ROTATION (perifocal -> ECI)
R_oi = [1 0 0; ...
        0 cos(incl) -sin(incl); ...
        0 sin(incl)  cos(incl)];

%% SUN DIRECTION (ECI, fixed)
s_eci_const = [1; 0; 0];
s_eci_const = s_eci_const / norm(s_eci_const);

%% MTQ PROPERTIES
m_max_axis = 10;                  % max dipole per axis [A·m^2]
m_res_body = [0.05; 0.05; 0.05];  % residual dipole [A·m^2]

%% CONTROLLER GAINS (tune as needed)
Kp = 5e-4;   % proportional on sun-direction error
Kd = 1e-2;   % derivative on angular rate

%% INITIAL ATTITUDE: sun-aligned as before
% First compute orbit position & velocity at t = 0
theta0 = 0;
r_orb0 = [a*cos(theta0); a*sin(theta0); 0];
v_orb0 = [-a*n*sin(theta0); a*n*cos(theta0); 0];
r_eci0 = R_oi * r_orb0;
v_eci0 = R_oi * v_orb0;

% Orbit angular momentum
h_eci = cross(r_eci0, v_eci0); h_eci = h_eci / norm(h_eci);

% Build body frame in ECI:
%   b_x aligned with Sun
%   b_z close to orbit normal and orthogonal to b_x
b_x = s_eci_const;
tmp = h_eci - dot(h_eci, b_x)*b_x;
if norm(tmp) < 1e-8
    tmp = [0;1;0] - dot([0;1;0], b_x)*b_x;
end
b_z = tmp / norm(tmp);
b_y = cross(b_z, b_x); b_y = b_y / norm(b_y);

C_bi0 = [b_x.'; b_y.'; b_z.'];        % inertial -> body at t=0
q0    = dcm_to_quat(C_bi0);           % quaternion (scalar-first)

% Add a small initial error if you want to see convergence:
% (uncomment if desired)
% dq = axang_to_quat([0 1 0 deg2rad(5)]);  % 5 deg about body y
% q0 = quat_multiply(dq, q0);
% q0 = q0 / norm(q0);

omega0 = [0; 0; 0];                   % initial angular rate [rad/s]

x0 = [q0; omega0];                    % state: [q0;q1;q2;q3; wx;wy;wz]

%% PACK PARAMETERS
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

p.m_max_axis = m_max_axis;
p.m_res_body = m_res_body;
p.s_eci_const = s_eci_const;

p.Kp = Kp;
p.Kd = Kd;

%% INTEGRATE ATTITUDE DYNAMICS
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t, x] = ode45(@(t,x) sat_dyn(t, x, p), [0 Tsim], x0, opts);

q_hist = x(:,1:4).';
w_hist = x(:,5:7).';

% Renormalize quaternions just in case
for k = 1:length(t)
    q_hist(:,k) = q_hist(:,k) / norm(q_hist(:,k));
end

%% POST-PROCESS: compute errors, torques, commanded dipoles
Nt = length(t);
sun_err_deg = zeros(1,Nt);
tau_dist_body = zeros(3,Nt);
tau_mtq_body  = zeros(3,Nt);
m_cmd_body    = zeros(3,Nt);

for k = 1:Nt
    tk = t(k);
    qk = q_hist(:,k);
    wk = w_hist(:,k);

    C_bi = quat_to_dcm_bi(qk);  % inertial -> body

    % Orbit position and velocity at time tk
    theta = p.n * tk;
    r_orb = [p.a*cos(theta); p.a*sin(theta); 0];
    v_orb = [-p.a*p.n*sin(theta); p.a*p.n*cos(theta); 0];

    r_eci = p.R_oi * r_orb;
    v_eci = p.R_oi * v_orb;

    % Magnetic field
    B_eci  = earth_dipole_B(r_eci);
    B_body = C_bi * B_eci;

    % Disturbances (same as in dynamics)
    [tau_gg, tau_drag, tau_srp] = dist_torques(r_eci, v_eci, C_bi, p);
    tau_res = cross(p.m_res_body, B_body);

    tau_dist = tau_gg + tau_drag + tau_srp + tau_res;
    tau_dist_body(:,k) = tau_dist;

    % Controller (same law as in dynamics)
    s_b = C_bi * p.s_eci_const;
    s_b = s_b / norm(s_b);
    s_des = [1;0;0];

    e_sun = cross(s_b, s_des);
    tau_cmd = -p.Kp * e_sun - p.Kd * wk;

    if norm(B_body) < 1e-10
        m_cmd = [0;0;0];
    else
        m_cmd = cross(B_body, tau_cmd) / (norm(B_body)^2);
        m_cmd = max(min(m_cmd, p.m_max_axis), -p.m_max_axis);  % saturate per axis
    end

    tau_mtq = cross(m_cmd, B_body);

    tau_mtq_body(:,k) = tau_mtq;
    m_cmd_body(:,k)   = m_cmd;

    % Sun alignment error
    s_des_b = [1;0;0];
    cosang  = dot(s_b, s_des_b);
    cosang  = max(min(cosang,1),-1);
    sun_err_deg(k) = acos(cosang) * 180/pi;
end

%% PLOTS

% Sun alignment error
figure;
plot(t/60, sun_err_deg);
xlabel('Time [min]');
ylabel('Sun alignment error [deg]');
title('Sun alignment error (body +x vs Sun)');
grid on;

% Angular rates
figure;
plot(t/60, w_hist(1,:), ...
     t/60, w_hist(2,:), ...
     t/60, w_hist(3,:));
xlabel('Time [min]');
ylabel('\omega [rad/s]');
title('Body angular rates');
legend('\omega_x','\omega_y','\omega_z');
grid on;

% Disturbance vs MTQ torque per axis
figure;
subplot(3,1,1);
plot(t/60, tau_dist_body(1,:), 'r', t/60, tau_mtq_body(1,:), 'b');
ylabel('\tau_x [N·m]');
title('Torque per axis: disturbance vs MTQ');
legend('Disturbance','MTQ');
grid on;

subplot(3,1,2);
plot(t/60, tau_dist_body(2,:), 'r', t/60, tau_mtq_body(2,:), 'b');
ylabel('\tau_y [N·m]');
legend('Disturbance','MTQ');
grid on;

subplot(3,1,3);
plot(t/60, tau_dist_body(3,:), 'r', t/60, tau_mtq_body(3,:), 'b');
ylabel('\tau_z [N·m]');
xlabel('Time [min]');
legend('Disturbance','MTQ');
grid on;

% Commanded dipoles
figure;
plot(t/60, m_cmd_body(1,:), ...
     t/60, m_cmd_body(2,:), ...
     t/60, m_cmd_body(3,:));
xlabel('Time [min]');
ylabel('m_{cmd} [A·m^2]');
title('Commanded MTQ dipoles');
legend('m_x','m_y','m_z');
grid on;

%% Total torque magnitudes
tau_dist_mag = sqrt(sum(tau_dist_body.^2, 1));     % |disturbance|
tau_mtq_mag  = sqrt(sum(tau_mtq_body.^2, 1));      % |MTQ|
tau_total_mag = sqrt(sum((tau_dist_body + tau_mtq_body).^2, 1));   % |dist + MTQ|

figure; hold on;
h1 = plot(t/60, tau_dist_mag,  'r',      'LineWidth', 1.5);
h2 = plot(t/60, tau_mtq_mag,   'b',      'LineWidth', 1.5);

xlabel('Time [min]');
ylabel('Torque magnitude [N·m]');
title('Torque magnitudes: disturbance, MTQ');
legend([h1 h2], {'|τ_{dist}|','|τ_{MTQ}|'});
grid on;
hold off;



%% ================= LOCAL FUNCTIONS =================

function dx = sat_dyn(t, x, p)
    % State
    q = x(1:4);
    w = x(5:7);

    q = q / norm(q);  % keep unit

    % Attitude matrix
    C_bi = quat_to_dcm_bi(q);  % inertial -> body

    % Orbit position & velocity (circular)
    theta = p.n * t;
    r_orb = [p.a*cos(theta); p.a*sin(theta); 0];
    v_orb = [-p.a*p.n*sin(theta); p.a*p.n*cos(theta); 0];

    r_eci = p.R_oi * r_orb;
    v_eci = p.R_oi * v_orb;

    % Magnetic field
    B_eci  = earth_dipole_B(r_eci);
    B_body = C_bi * B_eci;

    % Disturbance torques
    [tau_gg, tau_drag, tau_srp] = dist_torques(r_eci, v_eci, C_bi, p);
    tau_res = cross(p.m_res_body, B_body);
    tau_dist = tau_gg + tau_drag + tau_srp + tau_res;

    % Controller: sun tracking
    s_b = C_bi * p.s_eci_const;
    s_b = s_b / norm(s_b);
    s_des = [1;0;0];               % want +x axis to face the Sun

    e_sun  = cross(s_b, s_des);    % direction to rotate
    tau_cmd = -p.Kp * e_sun - p.Kd * w;

    % Map commanded torque to MTQ dipole
    if norm(B_body) < 1e-10
        m_cmd = [0;0;0];
    else
        m_cmd = cross(B_body, tau_cmd) / (norm(B_body)^2);
        m_cmd = max(min(m_cmd, p.m_max_axis), -p.m_max_axis);  % saturate per axis
    end

    tau_mtq = cross(m_cmd, B_body);

    tau_total = tau_dist + tau_mtq;

    % Rotational dynamics
    wdot = p.J \ (tau_total - cross(w, p.J*w));

    % Quaternion kinematics (scalar-first)
    wx = w(1); wy = w(2); wz = w(3);
    Omega = [ 0   -wx  -wy  -wz;
              wx   0    wz  -wy;
              wy  -wz   0    wx;
              wz   wy  -wx   0  ];
    qdot = 0.5 * Omega * q;

    dx = [qdot; wdot];
end

function [tau_gg, tau_drag, tau_srp] = dist_torques(r_eci, v_eci, C_bi, p)
    r_norm = norm(r_eci);
    r_hat_i = -r_eci / r_norm;     % nadir
    r_hat_b = C_bi * r_hat_i;

    % Gravity gradient
    tau_gg = 3*p.muE/r_norm^3 * cross(r_hat_b, p.J*r_hat_b);

    % Atmosphere density
    h_inst = r_norm - p.Re;
    rho    = p.rho_ref * exp(-(h_inst - p.h_ref)/p.H_scale);

    % Relative velocity
    v_atm_eci = cross(p.omegaE, r_eci);
    v_rel_eci = v_eci - v_atm_eci;
    v_rel_body = C_bi * v_rel_eci;

    v_mag = norm(v_rel_body);
    if v_mag < 1e-6
        tau_drag = [0;0;0];
    else
        v_hat = v_rel_body / v_mag;
        Fd = -0.5 * rho * v_mag^2 * p.Cd * p.A_drag_ref * v_hat;
        tau_drag = cross(p.r_cp_drag_body, Fd);
    end

    % Solar radiation pressure
    s_body = C_bi * p.s_eci_const;
    s_body = s_body / norm(s_body);
    F_srp = -p.P_srp * p.Cr * p.A_srp_ref * s_body;
    tau_srp = cross(p.r_cp_srp_body, F_srp);
end

function B = earth_dipole_B(r_eci)
    % Simple centered dipole aligned with +Z, magnitude tuned for LEO
    M   = 7.96e15;  % [T·m^3]
    r   = norm(r_eci);
    er  = r_eci / r;
    m_hat = [0;0;1];

    coef = M / r^3;
    mdot = dot(m_hat, er);
    B = coef * (3*mdot*er - m_hat);
end

function C = quat_to_dcm_bi(q)
    % Convert quaternion (scalar-first) to DCM inertial->body
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    C = [ 1-2*(q2^2+q3^2),   2*(q1*q2+q0*q3),   2*(q1*q3-q0*q2);
          2*(q1*q2-q0*q3),   1-2*(q1^2+q3^2),   2*(q2*q3+q0*q1);
          2*(q1*q3+q0*q2),   2*(q2*q3-q0*q1),   1-2*(q1^2+q2^2)];
end

function q = dcm_to_quat(C)
    % Convert DCM inertial->body to quaternion (scalar-first)
    tr = trace(C);
    if tr > 0
        S = sqrt(tr+1.0)*2;
        q0 = 0.25 * S;
        q1 = (C(3,2) - C(2,3)) / S;
        q2 = (C(1,3) - C(3,1)) / S;
        q3 = (C(2,1) - C(1,2)) / S;
    else
        if (C(1,1) > C(2,2)) && (C(1,1) > C(3,3))
            S = sqrt(1.0 + C(1,1) - C(2,2) - C(3,3))*2;
            q0 = (C(3,2) - C(2,3)) / S;
            q1 = 0.25 * S;
            q2 = (C(1,2) + C(2,1)) / S;
            q3 = (C(1,3) + C(3,1)) / S;
        elseif C(2,2) > C(3,3)
            S = sqrt(1.0 + C(2,2) - C(1,1) - C(3,3))*2;
            q0 = (C(1,3) - C(3,1)) / S;
            q1 = (C(1,2) + C(2,1)) / S;
            q2 = 0.25 * S;
            q3 = (C(2,3) + C(3,2)) / S;
        else
            S = sqrt(1.0 + C(3,3) - C(1,1) - C(2,2))*2;
            q0 = (C(2,1) - C(1,2)) / S;
            q1 = (C(1,3) + C(3,1)) / S;
            q2 = (C(2,3) + C(3,2)) / S;
            q3 = 0.25 * S;
        end
    end
    q = [q0;q1;q2;q3] / norm([q0;q1;q2;q3]);
end
