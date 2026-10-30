%% SSO satellite, inertially fixed & sun-aligned
%% Disturbance torques + MTQ authority + sun alignment error
%% Simulation: 2 orbits

clear; clc; close all;

%% CONSTANTS
Re   = 6378.137e3;           % Earth radius [m]
muE  = 3.986004418e14;       % Earth mu [m^3/s^2]
h0   = 510e3;                % nominal altitude [m]
a    = Re + h0;              % semi-major axis [m]
incl = deg2rad(97.4);        % SSO-ish inclination [rad]

% Orbit period
n    = sqrt(muE/a^3);        % mean motion [rad/s]
Torb = 2*pi/n;               % orbital period [s]

% Simulation length: 2 complete orbits
Tsim = 2*Torb;

% Time grid
Npts  = 3000;
tspan = linspace(0, Tsim, Npts);

% Earth rotation
omega_E = [0; 0; 7.2921159e-5];   % [rad/s]

%% SATELLITE PROPERTIES
m_sat = 40;                 % mass [kg]

% Bus dimensions
Lx = 0.85;
Ly = 0.40;
Lz = 0.50;

% Solar arrays (two wings, each 0.85 x 0.40)
A_array = 0.85 * 0.40;
A_arrays_total = 2 * A_array;

% Reference areas (rough)
A_drag_ref = 0.85 * 0.40 + A_arrays_total;
A_srp_ref  = 0.40 * 0.50 + A_arrays_total;

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

%% ORBIT IN ECI (circular)
theta = n * tspan;

r_orb = [ a*cos(theta) ; ...
          a*sin(theta) ; ...
          zeros(size(theta)) ];

v_orb = [ -a*n*sin(theta) ; ...
           a*n*cos(theta) ; ...
           zeros(size(theta)) ];

R_oi = [1 0 0; ...
        0 cos(incl) -sin(incl); ...
        0 sin(incl)  cos(incl)];

r_eci = R_oi * r_orb;
v_eci = R_oi * v_orb;

%% SUN DIRECTION (ECI)
s_eci_const = [1; 0; 0];
s_eci_const = s_eci_const / norm(s_eci_const);

%% ATTITUDE (INERTIALLY FIXED, SUN-ALIGNED)
h_eci = cross(r_eci(:,1), v_eci(:,1));
h_eci = h_eci / norm(h_eci);

b_x = s_eci_const;
tmp = h_eci - dot(h_eci, b_x)*b_x;
if norm(tmp) < 1e-8
    tmp = [0;1;0] - dot([0;1;0], b_x)*b_x;
end
b_z = tmp / norm(tmp);
b_y = cross(b_z, b_x);
b_y = b_y / norm(b_y);

% inertial -> body DCM
C_bi = [b_x.'; b_y.'; b_z.'];

%% MTQs
m_max_axis = 15;                  % max dipole per axis [A·m^2]
m_res_body = [0.05; 0.05; 0.05];  % residual dipole [A·m^2]

%% STORAGE
B_eci  = zeros(3, Npts);
B_body = zeros(3, Npts);

alpha_max_body = zeros(3, Npts);

tau_gg_body   = zeros(3, Npts);
tau_drag_body = zeros(3, Npts);
tau_srp_body  = zeros(3, Npts);
tau_res_body  = zeros(3, Npts);

%% MAIN LOOP
for k = 1:Npts
    r_k = r_eci(:,k);
    v_k = v_eci(:,k);
    r_norm = norm(r_k);

    %% Magnetic field
    B_eci(:,k)  = earth_dipole_B(r_k);
    B_body(:,k) = C_bi * B_eci(:,k);

    %% MTQ maximum authority (per-axis alpha magnitude, not directional yet)
    for i = 1:3
        m_body = zeros(3,1);
        m_body(i) = m_max_axis;
        tau_body  = cross(m_body, B_body(:,k));
        alpha     = J \ tau_body;
        alpha_max_body(i,k) = norm(alpha);
    end

    %% Gravity-gradient torque
    r_hat_i = -r_k / r_norm;        % nadir direction (ECI)
    r_hat_b = C_bi * r_hat_i;       % in body frame
    tau_gg_body(:,k) = 3*muE/r_norm^3 * cross(r_hat_b, J*r_hat_b);

    %% Aerodynamic drag torque
    h_inst  = r_norm - Re;
    rho     = rho_ref * exp(-(h_inst - h_ref)/H_scale);
    rho;

    v_atm_eci = cross(omega_E, r_k);
    v_rel_eci = v_k - v_atm_eci;
    v_rel_body = C_bi * v_rel_eci;

    v_mag = norm(v_rel_body);
    if v_mag < 1e-6
        tau_drag_body(:,k) = [0;0;0];
    else
        v_hat = v_rel_body / v_mag;
        Fd = -0.5 * rho * v_mag^2 * Cd * A_drag_ref * v_hat;
        tau_drag_body(:,k) = cross(r_cp_drag_body, Fd);
    end

    %% Solar radiation pressure torque
    s_body = C_bi * s_eci_const;
    s_body = s_body / norm(s_body);
    F_srp = -P_srp * Cr * A_srp_ref * s_body;
    tau_srp_body(:,k) = cross(r_cp_srp_body, F_srp);

    %% Residual magnetic dipole torque
    tau_res_body(:,k) = cross(m_res_body, B_body(:,k));
end

%% TOTAL DISTURBANCE TORQUE
tau_dist_body = tau_gg_body + tau_drag_body + tau_srp_body + tau_res_body;
tau_dist_mag  = vecnorm(tau_dist_body);

%% NEW: MAX MTQ TORQUE ALONG DISTURBANCE DIRECTION
tau_mtq_max_along_dist = zeros(1, Npts);

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

for k = 1:Npts
    Bd = B_body(:,k);
    Bd_norm = norm(Bd);
    if Bd_norm < 1e-12 || tau_dist_mag(k) < 1e-16
        tau_mtq_max_along_dist(k) = 0;
        continue;
    end

    % Unit vector along total disturbance torque
    u = tau_dist_body(:,k) / tau_dist_mag(k);

    % Each axis contributes tau_i = m_i * (e_i × B)
    c1 = cross(ex, Bd);
    c2 = cross(ey, Bd);
    c3 = cross(ez, Bd);

    % Optimal m_i for maximizing projection along u:
    % maximize u·(Σ m_i c_i) with |m_i| <= m_max_axis
    % => m_i = m_max_axis * sign(u·c_i)
    m1 = m_max_axis * sign(dot(u, c1));
    m2 = m_max_axis * sign(dot(u, c2));
    m3 = m_max_axis * sign(dot(u, c3));

    tau_max_vec = m1*c1 + m2*c2 + m3*c3;

    % Component along disturbance direction
    tau_mtq_max_along_dist(k) = max(0, dot(u, tau_max_vec));
end

%% NEW: MTQ max torque per body axis
tau_mtq_max_axis = zeros(3, Npts);   % rows: x,y,z

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

for k = 1:Npts
    Bd = B_body(:,k);

    % If field is too small, no torque possible
    if norm(Bd) < 1e-12
        tau_mtq_max_axis(:,k) = 0;
        continue;
    end

    % Basis torque directions for 1 A·m² along each coil axis
    c1 = cross(ex, Bd);
    c2 = cross(ey, Bd);
    c3 = cross(ez, Bd);

    % Optimal ±m to maximize torque along each body axis
    % axis x
    u = ex;                  % desired torque direction
    m1 = m_max_axis * sign(dot(u, c1));
    m2 = m_max_axis * sign(dot(u, c2));
    m3 = m_max_axis * sign(dot(u, c3));
    tau_vec = m1*c1 + m2*c2 + m3*c3;
    tau_mtq_max_axis(1,k) = dot(u, tau_vec);

    % axis y
    u = ey;
    m1 = m_max_axis * sign(dot(u, c1));
    m2 = m_max_axis * sign(dot(u, c2));
    m3 = m_max_axis * sign(dot(u, c3));
    tau_vec = m1*c1 + m2*c2 + m3*c3;
    tau_mtq_max_axis(2,k) = dot(u, tau_vec);

    % axis z
    u = ez;
    m1 = m_max_axis * sign(dot(u, c1));
    m2 = m_max_axis * sign(dot(u, c2));
    m3 = m_max_axis * sign(dot(u, c3));
    tau_vec = m1*c1 + m2*c2 + m3*c3;
    tau_mtq_max_axis(3,k) = dot(u, tau_vec);
end

% Disturbance per axis
tau_dist_axis = tau_dist_body;   % just rename for clarity


%% Magnitudes for other plots
Bmag_body     = vecnorm(B_body);
tau_gg_mag    = vecnorm(tau_gg_body);
tau_drag_mag  = vecnorm(tau_drag_body);
tau_srp_mag   = vecnorm(tau_srp_body);
tau_res_mag   = vecnorm(tau_res_body);
tau_total_mag = tau_dist_mag;      % same as |τ_gg + τ_drag + τ_srp + τ_res|

%% SUN ALIGNMENT ERROR (still kinematic)
b1_eci = C_bi.' * [1;0;0];
b1_eci = b1_eci / norm(b1_eci);
s_eci  = s_eci_const;

err_sun = zeros(1, Npts);
for k = 1:Npts
    cosang = dot(b1_eci, s_eci);
    cosang = max(min(cosang,1),-1);
    err_sun(k) = acos(cosang);
end
err_sun_deg = err_sun * 180/pi;

%% PLOTS

% Magnetic field magnitude
figure;
plot(tspan/60, Bmag_body*1e6);
xlabel('Time [min]');
ylabel('|B_{body}| [\muT]');
title('Magnetic field magnitude in body frame');
grid on;

% MTQ max angular acceleration per axis
figure;
plot(tspan/60, alpha_max_body(1,:), ...
     tspan/60, alpha_max_body(2,:), ...
     tspan/60, alpha_max_body(3,:));
xlabel('Time [min]');
ylabel('|\alpha_{max}| [rad/s^2]');
title('MTQ maximum angular acceleration (15 A·m^2 per axis)');
legend('Body X','Body Y','Body Z');
grid on;

% Old breakdown of disturbance torques
figure;
plot(tspan/60, tau_gg_mag, ...
     tspan/60, tau_drag_mag, ...
     tspan/60, tau_srp_mag, ...
     tspan/60, tau_res_mag, ...
     tspan/60, tau_total_mag);
xlabel('Time [min]');
ylabel('|Torque| [N·m]');
title('Disturbance torque magnitudes');
legend('Gravity gradient','Aerodynamic drag','SRP','Residual dipole','Total');
grid on;

%% NEW: Combined disturbance vs max MTQ torque *in that direction*
figure;
plot(tspan/60, tau_dist_mag, ...
     tspan/60, tau_mtq_max_along_dist, 'LineWidth', 1.2);
xlabel('Time [min]');
ylabel('Torque [N·m]');
title('Disturbance torque vs max MTQ torque along disturbance direction');
legend('|τ_{dist}|','Max |τ_{MTQ}| along τ_{dist}');
grid on;

% Sun alignment error
figure;
plot(tspan/60, err_sun_deg);
xlabel('Time [min]');
ylabel('Sun alignment error [deg]');
title('Sun alignment error');
grid on;

%% ====== LOCAL FUNCTION ======
function B = earth_dipole_B(r_eci)
    % Simple centered dipole field aligned with ECI z-axis
    %   r_eci : 3x1 position [m]
    %   B     : 3x1 magnetic field [T]
    %
    % M here is in T·m^3, so DO NOT multiply by mu0/(4*pi).

    M   = 7.96e15;          % [T·m^3], tuned so |B| ~ 3e-5 T at Earth's surface

    r  = norm(r_eci);
    er = r_eci / r;

    m_hat = [0; 0; 1];      % dipole aligned with +Z

    coef  = M / r^3;        % <--- removed mu0/(4*pi)
    mdot  = dot(m_hat, er);

    B = coef * (3*mdot*er - m_hat);  % [T]
end
figure;
subplot(3,1,1);
plot(tspan/60, tau_dist_axis(1,:), 'r', 'LineWidth', 1); hold on;
plot(tspan/60, tau_mtq_max_axis(1,:), 'b', 'LineWidth', 1.2);
ylabel('\tau_x [N·m]');
legend('Disturbance','MTQ max');
title('Torque comparison per axis');
grid on;

subplot(3,1,2);
plot(tspan/60, tau_dist_axis(2,:), 'r', 'LineWidth', 1); hold on;
plot(tspan/60, tau_mtq_max_axis(2,:), 'b', 'LineWidth', 1.2);
ylabel('\tau_y [N·m]');
legend('Disturbance','MTQ max');
grid on;

subplot(3,1,3);
plot(tspan/60, tau_dist_axis(3,:), 'r', 'LineWidth', 1); hold on;
plot(tspan/60, tau_mtq_max_axis(3,:), 'b', 'LineWidth', 1.2);
ylabel('\tau_z [N·m]');
xlabel('Time [min]');
legend('Disturbance','MTQ max');
grid on;
