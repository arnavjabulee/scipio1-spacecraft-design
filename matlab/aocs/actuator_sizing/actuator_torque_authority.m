clear; clc;

%% ORBIT & CONSTANTS (same flavour as your existing script)
mu  = 3.986e14;        % [m^3/s^2]
Re  = 6371e3;          % [m]
alt = 510e3;           % [m]
r   = Re + alt;

incl     = deg2rad(97.5);   % SSO inclination
tilt_mag = deg2rad(11);     % geomagnetic dipole tilt

% Approximate |B| at this altitude (constant magnitude, varying direction)
B_mag = 3e-5;               % [T] ~30 µT

% For now: body frame aligned with ECI
R_body_to_ECI = eye(3);

%% ACTUATORS
% ---------- Magnetorquers: 3 x CR0100 ----------
% CubeSpace CR0100 max dipole ≈ 10 A·m² at full drive (adjust if needed)
m_CR0100 = 10;              % [A·m²] per rod
n_MTQ    = 3;               % 3 orthogonal rods (X, Y, Z)

% ---------- Reaction wheels: 4 x CW0500 ----------
% CW0500 max torque ≈ 16 mN·m per wheel
tau_RW_single = 16e-3;      % [N·m]

% Assume 3 RW giving torque about X, Y, Z (4th is redundant)
n_RW_axes   = 3;
tau_RW_axis = tau_RW_single;   % [N·m] max about each body axis

%% ORBIT SAMPLING & MAGNETIC FIELD
n_points    = 500;
theta_orbit = linspace(0, 2*pi, n_points);   % true anomaly

% Geomagnetic dipole axis in ECI (unit-ish vector)
m_eci = [sin(tilt_mag); 0; cos(tilt_mag)];

% Preallocate
tau_mag_x = zeros(1, n_points);
tau_mag_y = zeros(1, n_points);
tau_mag_z = zeros(1, n_points);

tau_tot_x = zeros(1, n_points);
tau_tot_y = zeros(1, n_points);
tau_tot_z = zeros(1, n_points);
tau_tot_min_axis = zeros(1, n_points);

for i = 1:n_points
    u = theta_orbit(i);

    % Orbit position unit vector in ECI (like your earlier code)
    r_orb = [cos(u);
             sin(u)*cos(incl);
             sin(u)*sin(incl)];

    % Dipole magnetic field direction in ECI (up to a scalar)
    Bvec_eci  = 3*(dot(m_eci, r_orb))*r_orb - m_eci;
    Bvec_body = R_body_to_ECI' * Bvec_eci;   % into body frame

    % Normalise
    Bhat_body = Bvec_body / norm(Bvec_body);

    % --- Magnetorquer torque about each body axis ---
    % Each CR0100 rod aligned with +X, +Y, +Z.
    cos_th_x = abs(Bhat_body(1));
    cos_th_y = abs(Bhat_body(2));
    cos_th_z = abs(Bhat_body(3));

    sin_th_x = sqrt(1 - cos_th_x^2);
    sin_th_y = sqrt(1 - cos_th_y^2);
    sin_th_z = sqrt(1 - cos_th_z^2);

    tau_mag_x(i) = m_CR0100 * B_mag * sin_th_x;
    tau_mag_y(i) = m_CR0100 * B_mag * sin_th_y;
    tau_mag_z(i) = m_CR0100 * B_mag * sin_th_z;

    % --- Total torque per axis (RW + MTQ) ---
    tau_tot_x(i) = tau_RW_axis + tau_mag_x(i);
    tau_tot_y(i) = tau_RW_axis + tau_mag_y(i);
    tau_tot_z(i) = tau_RW_axis + tau_mag_z(i);

    % Worst axis at this orbit position
    tau_tot_min_axis(i) = min([tau_tot_x(i), tau_tot_y(i), tau_tot_z(i)]);
end

theta_deg = rad2deg(theta_orbit);

%% PLOTS
figure('Color','w');

subplot(2,1,1);
plot(theta_deg, tau_mag_x*1e3, 'LineWidth', 1.5); hold on;
plot(theta_deg, tau_mag_y*1e3, 'LineWidth', 1.5);
plot(theta_deg, tau_mag_z*1e3, 'LineWidth', 1.5);
xlabel('True anomaly [deg]');
ylabel('Magnetorquer torque [mN·m]');
title('Magnetorquer Max Torque vs Orbit Position (3 x CR0100)');
legend('X-axis','Y-axis','Z-axis','Location','best');
grid on;

subplot(2,1,2);
plot(theta_deg, tau_tot_x*1e3, 'LineWidth', 1.5); hold on;
plot(theta_deg, tau_tot_y*1e3, 'LineWidth', 1.5);
plot(theta_deg, tau_tot_z*1e3, 'LineWidth', 1.5);
plot(theta_deg, tau_tot_min_axis*1e3, '--k', 'LineWidth', 1.5);
xlabel('True anomaly [deg]');
ylabel('Total available torque [mN·m]');
title('Total Available Torque (RW + MTQ) vs Orbit Position');
legend('X-axis','Y-axis','Z-axis','Min over axes','Location','best');
grid on;

%% SUMMARY
fprintf('--- Torque Summary (per axis) ---\n');
fprintf('Reaction wheel torque per axis (CW0500): %.2f mN·m\n', ...
        tau_RW_axis*1e3);
fprintf('Magnetorquer peak torque about an axis (CR0100, best case): %.2f mN·m\n', ...
        max([tau_mag_x, tau_mag_y, tau_mag_z])*1e3);
fprintf('Minimum total torque across all axes over the orbit: %.2f mN·m\n', ...
        min(tau_tot_min_axis)*1e3);
