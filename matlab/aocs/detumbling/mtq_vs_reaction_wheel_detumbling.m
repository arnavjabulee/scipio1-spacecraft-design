%% Detumble from 2 deg/s: MTQs-only vs RWs-only (only angular velocity plots)
clear; clc;

%% ORBIT & CONSTANTS
Re   = 6371e3; mu = 3.986e14; B_eq = 30e-6;
altitude = 510e3; r = Re + altitude;
B_mag = B_eq * (Re/r)^3;         % [T]
n_orbit = sqrt(mu/r^3); T_orbit = 2*pi/n_orbit;

%% INERTIA (kg·m^2)
J = [  1.991   -0.04488  -0.03172 ; ...
      -0.04488  2.721     0.03037 ; ...
      -0.03172  0.03037   1.329  ];
Javg = trace(J)/3;                % use average inertia for gain sizing

%% ACTUATORS
m_mtq_max = [10; 10; 10];        % [A·m^2] MTQ limits per axis
T_rw_peak = 20e-3;               % [N·m] RW torque limit per axis

%% CONTROLLER GAINS
% Target: 2 deg/s -> 0.01 deg/s in <= 20 min
omega_ratio = deg2rad(2)/deg2rad(0.01);   % 200
t_target = 20*60;                          % 1200 s
tau_target = t_target / log(omega_ratio);  % ~226 s

% B-dot gain from linearized tau = J/(k B^2); add 10% margin
k_Bdot = 1.1 * (Javg / (tau_target * B_mag^2));  % ~1.6e7 A·m^2 per (T/s)

% For reference, print expected initial saturation check
fprintf('B = %.2e T, tau_target ≈ %.1f s, k_Bdot = %.2e\n', B_mag, tau_target, k_Bdot);

K_rw = 0.1;                        % RW rate-damping gain [N·m·s/rad]

%% INITIAL RATES
omega0_mag = deg2rad(2);
dir0 = [1;1;0.5]; dir0 = dir0/norm(dir0);
omega0 = omega0_mag * dir0;
omega_tol = deg2rad(0.01);

%% SIM SETTINGS
dt = 0.1; t_max = 10*T_orbit; Nmax = ceil(t_max/dt);

%% ------------------------ CASE 1 — MTQs ONLY (B-dot) ------------------------
omega = omega0;
B_body = [0;0;B_mag];

t_mtq = zeros(1,Nmax); w_mtq = zeros(3,Nmax);
m_cmd_log = zeros(3,Nmax); m_clip_log = zeros(3,Nmax);
sat_flags  = false(3,Nmax);

kend_mtq = Nmax;
for k = 1:Nmax
    t = (k-1)*dt;
    t_mtq(k) = t; w_mtq(:,k) = omega;

    if norm(omega) <= omega_tol, kend_mtq = k; break; end

    % B-dot: dB/dt = -omega x B
    Bdot_body = -cross(omega, B_body);

    % Commanded dipole and per-axis clipping
    m_cmd = -k_Bdot * Bdot_body;  m_cmd_log(:,k) = m_cmd;
    m_clip = min(max(m_cmd, -m_mtq_max), m_mtq_max);  % clip component-wise
    sat_flags(:,k) = (abs(m_clip) == m_mtq_max);       % saturation flags
    m_clip_log(:,k) = m_clip;

    % Torque and rigid-body dynamics
    T = cross(m_clip, B_body);
    omega_dot = J \ (T - cross(omega, J*omega));
    omega = omega + omega_dot*dt;

    % Update B in body frame (keep |B| fixed)
    B_body = B_body + Bdot_body*dt;
    nb = norm(B_body); if nb>0, B_body = B_body * (B_mag/nb); else, B_body = [0;0;B_mag]; end
end

t_mtq = t_mtq(1:kend_mtq); w_mtq = w_mtq(:,1:kend_mtq);
m_cmd_log = m_cmd_log(:,1:kend_mtq); m_clip_log = m_clip_log(:,1:kend_mtq);
sat_flags  = sat_flags(:,1:kend_mtq);

%% -------------------- CASE 2 — REACTION WHEELS ONLY --------------------
omega = omega0;
t_rw = zeros(1,Nmax); w_rw = zeros(3,Nmax);
kend_rw = Nmax;
for k = 1:Nmax
    t = (k-1)*dt;
    t_rw(k) = t; w_rw(:,k) = omega;
    if norm(omega) <= omega_tol, kend_rw = k; break; end

    T_cmd = -K_rw * omega;
    T_cmd = min(max(T_cmd, -T_rw_peak), T_rw_peak);
    omega_dot = J \ (T_cmd - cross(omega,J*omega));
    omega = omega + omega_dot*dt;
end
t_rw = t_rw(1:kend_rw); w_rw = w_rw(:,1:kend_rw);

%% ------------------------------ PLOTS ------------------------------
figure; plot(t_mtq/60, vecnorm(w_mtq)*180/pi, 'LineWidth', 1.6);
xlabel('Time [min]'); ylabel('|ω| [deg/s]'); grid on;
title('Detumble with MTQs only (B-dot)');

figure; plot(t_rw/60, vecnorm(w_rw)*180/pi, 'LineWidth', 1.6);
xlabel('Time [min]'); ylabel('|ω| [deg/s]'); grid on;
title('Detumble with Reaction Wheels only (rate damping)');

sat_any = any(sat_flags,2);
fprintf('MTQ axis saturation detected?  x:%d  y:%d  z:%d\n', sat_any(1), sat_any(2), sat_any(3));
