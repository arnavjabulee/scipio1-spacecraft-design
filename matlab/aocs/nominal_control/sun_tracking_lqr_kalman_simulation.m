%% MTQ sun tracking: periodic state-space LQR + linear Kalman filter
% Same 510 km circular orbit, spacecraft inertia, environmental torque models
% and 10 A m^2 per-axis MTQ limit as sun_tracking_pd_simulation.m.
% Nonlinear quaternion truth model; linear, discrete-time estimator/controller.
% State x = [theta_x;theta_y;theta_z;omega_x;omega_y;omega_z], SI units.
% theta is a small rotation error about a fixed, sun-aligned reference.
% A fixed roll reference is chosen as well as body +X Sun alignment.
%
% Sensors: synthetic three-axis star-tracker attitude-error measurements and
% three-axis gyro rates, available every sample, independent zero-mean noise.
% This is NOT a sun-vector-only observer, an EKF, or an outage/bias model.
% The linear Kalman filter uses reference-orbit schedules A(k), B(k) and fixed H.
% A periodic Riccati recursion accounts for the rank-two magnetic actuator;
% an instantaneous frozen-field three-axis LQR would not be stabilizable.
% All Riccati/filter steps use base MATLAB; no Control System Toolbox needed.
% References for the discrete LQR cost and linear-estimator conventions:
% https://www.mathworks.com/help/control/ref/lti.dlqr.html
% https://www.mathworks.com/help/control/ref/ss.kalman.html
clear; clc; close all;
rng(1,'twister');

%% Scenario and parameters (matching the single-run PD baseline)
p.Re = 6378.137e3;
p.muE = 3.986004418e14;
h0 = 510e3;
p.a = p.Re + h0;
incl = deg2rad(97.4);
p.n = sqrt(p.muE/p.a^3);
Torb = 2*pi/p.n;
p.R_oi = [1 0 0; 0 cos(incl) -sin(incl); 0 sin(incl) cos(incl)];
p.J = [1.991 -0.04488 -0.03172; ...
      -0.04488 2.721 0.03037; -0.03172 0.03037 1.329];
p.mass = 40;
p.omegaE = [0;0;7.2921159e-5];
p.Cd = 2.2; p.Cr = 1.3; p.P_srp = 4.5e-6;
p.A_drag_ref = 0.85*0.40 + 2*0.85*0.40;
p.A_srp_ref = 0.40*0.50 + 2*0.85*0.40;
p.r_cp_drag_body = [0.20;0;0.10];
p.r_cp_srp_body = [-0.10;0;0.20];
p.rho_ref = 1e-12; p.h_ref = 400e3; p.H_scale = 60e3;
p.m_max_axis = 10;
p.m_res_body = [0.05;0.05;0.05];
p.s_eci_const = [1;0;0];

% Sample grid repeats exactly after one orbit. Commands are held per sample.
sample_time_requested = 2; % seconds; editable
samples_per_orbit = round(Torb/sample_time_requested);
dt = Torb/samples_per_orbit;
t = (0:2*samples_per_orbit)*dt;
Nt = numel(t);

% Select the same geometric reference as the PD script, with a consistent
% scalar-first quaternion/DCM convention: C_bi maps inertial vectors to body.
[r0,v0] = orbit_state(0,p);
orbit_normal = cross(r0,v0); orbit_normal = orbit_normal/norm(orbit_normal);
bx = p.s_eci_const;
bz = orbit_normal - dot(orbit_normal,bx)*bx; bz = bz/norm(bz);
by = cross(bz,bx);
C_ref = [bx.';by.';bz.'];
q_ref = dcm_to_quat(C_ref);
assert(norm(quat_to_dcm(q_ref)-C_ref,'fro') < 1e-12, ...
    'Quaternion/DCM convention mismatch.');

%% Periodic small-angle state-space model and magnetic input
% Linearize environmental torques about the known reference orbit/attitude.
% A_c(k)=[0 I; J^-1*d(tau_dist)/d(theta) 0],
% B_c(k)=[0; -J^-1*skew(B_ref(k))], d_c(k)=[0; J^-1*tau_ref(k)].
% Gyroscopic terms vanish at zero reference rate. Truth keeps all nonlinear
% terms. This is an offline scheduled LINEAR model, not state-based EKF
% relinearization. Exact frozen-coefficient ZOH uses one augmented expm.
A = zeros(6,6,samples_per_orbit);
B = zeros(6,3,samples_per_orbit);
B_reference = zeros(3,samples_per_orbit);
nominal_disturbance_step = zeros(6,samples_per_orbit);
G_torque = [0.5*dt^2*(p.J\eye(3)); dt*(p.J\eye(3))];
H = eye(6);
for k = 1:samples_per_orbit
    reference_time = (k-1)*dt;
    [nominal_torque,B_reference(:,k)] = environment(reference_time,C_ref,p);
    torque_jacobian = zeros(3);
    perturbation_size = 1e-5; % radians; central difference about reference
    for axis_index = 1:3
        perturbation = zeros(3,1);
        perturbation(axis_index) = perturbation_size;
        torque_plus = environment(reference_time, ...
            expm(-skew(perturbation))*C_ref,p);
        torque_minus = environment(reference_time, ...
            expm(skew(perturbation))*C_ref,p);
        torque_jacobian(:,axis_index) = ...
            (torque_plus-torque_minus)/(2*perturbation_size);
    end
    A_continuous = [zeros(3) eye(3);p.J\torque_jacobian zeros(3)];
    B_continuous = [zeros(3);-(p.J\skew(B_reference(:,k)))];
    d_continuous = [zeros(3,1);p.J\nominal_torque];
    transition = expm([A_continuous B_continuous d_continuous;zeros(4,10)]*dt);
    A(:,:,k) = transition(1:6,1:6);
    B(:,:,k) = transition(1:6,7:9);
    nominal_disturbance_step(:,k) = transition(1:6,10);
end

% Tuned for Sun-axis pointing on the common PD benchmark. A softer roll
% penalty reflects that roll about the Sun axis is not the primary objective.
% Increasing the input penalty avoids the original example's excessive
% dipole effort. The physical plant, initial attitude and limits are unchanged.
attitude_scale = deg2rad(1);
rate_scale = deg2rad(0.05);
roll_weight = 0.125;
dipole_scale = 0.25;
Q_lqr = diag([[roll_weight 1 1]/attitude_scale^2,ones(1,3)/rate_scale^2]);
R_lqr = eye(3)/dipole_scale^2;
[K_lqr,P_next_schedule,riccati_info] = periodic_lqr(A,B,Q_lqr,R_lqr);
assert(riccati_info.spectral_radius < 1, ...
    'The unsaturated linear periodic closed loop is not stable.');
% Optional model-based periodic disturbance preview. The affine Bellman
% recursion adds feedforward to the LQR feedback without using truth states.
% This assumes the deterministic reference-orbit torque model is known;
% disable it to examine pure regulation with unmodeled disturbances.
use_disturbance_feedforward = true;
[m_feedforward,nominal_disturbance_step] = periodic_feedforward( ...
    A,B,K_lqr,P_next_schedule,R_lqr,nominal_disturbance_step, ...
    use_disturbance_feedforward);

%% Linear Kalman-filter noise model
attitude_noise_std = deg2rad(0.03); % per axis, per sample
gyro_noise_std = deg2rad(0.003);    % rad/s, per axis, per sample
measurement_std = [repmat(attitude_noise_std,3,1); ...
                   repmat(gyro_noise_std,3,1)];
R_kf = diag(measurement_std.^2);
torque_noise_std = 5e-6; % Nm, piecewise-constant unknown torque per sample
Q_kf = G_torque*(torque_noise_std^2*eye(3))*G_torque.';
Q_kf = Q_kf + diag([repmat(1e-12,1,3),repmat(1e-14,1,3)]);
P_prior = diag([repmat(deg2rad(1)^2,1,3), ...
                repmat(deg2rad(0.02)^2,1,3)]);
x_prior = zeros(6,1);

%% Initial truth: a small, deliberately nonzero pointing/rate error
theta_initial = deg2rad([1;-1.5;0.5]);
omega_initial = deg2rad([0.005;-0.005;0.0025]);
q_initial = dcm_to_quat(expm(-skew(theta_initial))*C_ref);
truth = [q_initial;omega_initial];

%% Logs
q_hist = zeros(4,Nt); w_hist = zeros(3,Nt);
x_true = zeros(6,Nt); x_est = zeros(6,Nt); measurements = zeros(6,Nt);
P_diag = zeros(6,Nt); innovations = zeros(6,Nt);
sun_err_deg = zeros(1,Nt); attitude_err_deg = zeros(1,Nt);
tau_dist_body = zeros(3,Nt); tau_mtq_body = zeros(3,Nt);
m_cmd_body = zeros(3,Nt); B_body_hist = zeros(3,Nt);
saturation = false(3,Nt);
min_covariance_eigenvalue = inf;
max_covariance_asymmetry = 0;

%% Measurement update -> estimated-state control -> truth/filter prediction
for k = 1:Nt
    phase = mod(k-1,samples_per_orbit)+1;
    q = truth(1:4)/norm(truth(1:4));
    omega = truth(5:7);
    C_bi = quat_to_dcm(q);
    theta_true = attitude_error(C_bi,C_ref);
    x_true(:,k) = [theta_true;omega];
    measurements(:,k) = H*x_true(:,k) + measurement_std.*randn(6,1);

    % Linear KF with a current (posterior) estimate and Joseph covariance
    % update; do not use nonlinear truth or noise-free rates for feedback.
    S = H*P_prior*H.' + R_kf;
    L_kf = (P_prior*H.')/S;
    innovations(:,k) = measurements(:,k)-H*x_prior;
    x_posterior = x_prior + L_kf*innovations(:,k);
    I_LH = eye(6)-L_kf*H;
    P_posterior = I_LH*P_prior*I_LH.' + L_kf*R_kf*L_kf.';
    max_covariance_asymmetry = max(max_covariance_asymmetry, ...
        norm(P_posterior-P_posterior.','fro'));
    P_posterior = 0.5*(P_posterior+P_posterior.');
    min_covariance_eigenvalue = min(min_covariance_eigenvalue, ...
        min(eig(P_posterior)));

    m_requested = -K_lqr(:,:,phase)*x_posterior+m_feedforward(:,phase);
    m_applied = max(min(m_requested,p.m_max_axis),-p.m_max_axis);
    [tau_dist,B_body] = environment(t(k),C_bi,p);
    q_hist(:,k) = q; w_hist(:,k) = omega;
    x_est(:,k) = x_posterior; P_diag(:,k) = diag(P_posterior);
    m_cmd_body(:,k) = m_applied;
    saturation(:,k) = abs(m_requested)>p.m_max_axis;
    B_body_hist(:,k) = B_body;
    tau_dist_body(:,k) = tau_dist;
    tau_mtq_body(:,k) = cross(m_applied,B_body);
    sun_err_deg(k) = rad2deg(acos(max(-1,min(1, ...
        dot(C_bi*p.s_eci_const,[1;0;0])))));
    attitude_err_deg(k) = rad2deg(norm(theta_true));

    if k < Nt
        % Actual saturated dipole, not the unattainable LQR request, is the
        % known input to the filter. The plant uses the actual attitude and
        % field; the linear filter uses the reference-field input schedule.
        truth = rk4_step(t(k),truth,dt,m_applied,p);
        truth(1:4) = truth(1:4)/norm(truth(1:4));
        x_prior = A(:,:,phase)*x_posterior+B(:,:,phase)*m_applied ...
            +nominal_disturbance_step(:,phase);
        P_prior = A(:,:,phase)*P_posterior*A(:,:,phase).'+Q_kf;
        P_prior = 0.5*(P_prior+P_prior.');
    end
end

%% Numerical checks and summary (stability check excludes saturation/noise)
assert(all(isfinite(x_true(:))) && all(isfinite(x_est(:))), ...
    'Non-finite simulation states.');
assert(max(abs(vecnorm(q_hist)-1)) < 1e-10,'Quaternion norm drift.');
assert(max(abs(m_cmd_body(:))) <= p.m_max_axis+1e-12,'MTQ limit exceeded.');
assert(min_covariance_eigenvalue >= -1e-12,'Invalid filter covariance.');
orthogonality = abs(sum(tau_mtq_body.*B_body_hist,1));
assert(max(orthogonality) < 1e-15,'Magnetic torque must be normal to B.');
if max(attitude_err_deg) > 10
    warning('LQRKF:SmallAngleRange', ...
        'Attitude error exceeded 10 degrees; the linear design is approximate.');
end
last_orbit = t>=Torb;
metrics.initial_sun_error_deg = sun_err_deg(1);
metrics.final_sun_error_deg = sun_err_deg(end);
metrics.last_orbit_rms_sun_error_deg = sqrt(mean(sun_err_deg(last_orbit).^2));
metrics.max_attitude_error_deg = max(attitude_err_deg);
metrics.attitude_estimation_rmse_deg = rad2deg(sqrt(mean( ...
    (x_est(1:3,:)-x_true(1:3,:)).^2,2)));
metrics.rate_estimation_rmse_rad_s = sqrt(mean( ...
    (x_est(4:6,:)-x_true(4:6,:)).^2,2));
metrics.raw_attitude_measurement_rmse_deg = rad2deg(sqrt(mean( ...
    (measurements(1:3,:)-x_true(1:3,:)).^2,2)));
metrics.raw_rate_measurement_rmse_rad_s = sqrt(mean( ...
    (measurements(4:6,:)-x_true(4:6,:)).^2,2));
metrics.saturated_sample_fraction = mean(any(saturation(:,1:end-1),1));
metrics.min_covariance_eigenvalue = min_covariance_eigenvalue;
metrics.max_covariance_asymmetry = max_covariance_asymmetry;
metrics.linear_periodic_spectral_radius = riccati_info.spectral_radius;
fprintf('Periodic LQR + linear KF: dt = %.4f s, two orbits\n',dt);
fprintf('Riccati sweeps: %d; periodic residual: %.3g; spectral radius: %.6g\n', ...
    riccati_info.iterations,riccati_info.relative_residual, ...
    riccati_info.spectral_radius);
fprintf('Sun error: initial %.4f deg, final %.4f deg, last-orbit RMS %.4f deg\n', ...
    metrics.initial_sun_error_deg,metrics.final_sun_error_deg, ...
    metrics.last_orbit_rms_sun_error_deg);
fprintf('Samples with a saturated dipole axis: %.2f %%\n', ...
    100*metrics.saturated_sample_fraction);

%% Plots (same control outputs as the PD example, plus estimator diagnostics)
figure;
plot(t/60,sun_err_deg,t/60,attitude_err_deg);
xlabel('Time [min]'); ylabel('Error [deg]'); grid on;
legend('Body +X vs Sun','Full reference attitude');
title('Periodic LQR + linear KF: pointing error');

figure;
for axis_index = 1:3
    subplot(3,1,axis_index);
    plot(t/60,w_hist(axis_index,:),t/60,x_est(axis_index+3,:),'--');
    ylabel(sprintf('omega_%d [rad/s]',axis_index)); grid on;
    legend('Truth','Kalman estimate');
end
xlabel('Time [min]');

figure;
for axis_index = 1:3
    subplot(3,1,axis_index);
    plot(t/60,tau_dist_body(axis_index,:),t/60,tau_mtq_body(axis_index,:));
    ylabel(sprintf('tau_%d [N m]',axis_index)); grid on;
    legend('Disturbance','MTQ');
end
xlabel('Time [min]');

figure;
plot(t/60,m_cmd_body.'); hold on;
yline(p.m_max_axis,':','Color',[0.5 0.5 0.5]);
yline(-p.m_max_axis,':','Color',[0.5 0.5 0.5]);
xlabel('Time [min]'); ylabel('Dipole [A m^2]'); grid on;
legend('m_x','m_y','m_z'); title('Applied MTQ dipoles');

figure;
for axis_index = 1:3
    subplot(3,1,axis_index);
    sigma = sqrt(P_diag(axis_index,:));
    plot(t/60,rad2deg(x_est(axis_index,:)-x_true(axis_index,:)));
    hold on;
    plot(t/60,rad2deg(3*sigma),'--','Color',[0.85 0.33 0.10]);
    plot(t/60,rad2deg(-3*sigma),'--','Color',[0.85 0.33 0.10]);
    ylabel(sprintf('theta_%d error [deg]',axis_index)); grid on;
    legend('Estimation error','Modeled +3 sigma','Modeled -3 sigma');
end
xlabel('Time [min]');
sgtitle('Linear KF attitude error and modeled +/-3 sigma');

%% Local functions
function [gains,P_schedule,info] = periodic_lqr(A,B,Q,R)
    % Periodic Riccati recursion with A(k), B(k), P(N+1)=P(1).
    count = size(B,3); gains = zeros(3,6,count); boundary = Q;
    P_schedule = zeros(6,6,count);
    tolerance = 1e-10; max_iterations = 100;
    for sweep = 1:max_iterations
        P_next = boundary;
        for phase = count:-1:1
            input = B(:,:,phase);
            P_schedule(:,:,phase) = P_next;
            gain = (R+input.'*P_next*input)\(input.'*P_next*A(:,:,phase));
            closed_loop = A(:,:,phase)-input*gain;
            % Equivalent Riccati form, avoiding subtractive cancellation.
            P_now = Q+gain.'*R*gain+closed_loop.'*P_next*closed_loop;
            P_now = 0.5*(P_now+P_now.');
            gains(:,:,phase) = gain;
            P_next = P_now;
        end
        residual = norm(P_next-boundary,'fro')/max(1,norm(P_next,'fro'));
        boundary = P_next;
        if residual<tolerance, break; end
    end
    assert(residual<tolerance,'Periodic LQR Riccati iteration did not converge.');
    monodromy = eye(6);
    for phase = 1:count
        monodromy = (A(:,:,phase)-B(:,:,phase)*gains(:,:,phase))*monodromy;
    end
    info.iterations = sweep;
    info.relative_residual = residual;
    info.spectral_radius = max(abs(eig(monodromy)));
end

function [feedforward,disturbance] = periodic_feedforward( ...
        A,B,K,P_schedule,R,disturbance,enabled)
    count = size(B,3); feedforward = zeros(3,count);
    if ~enabled
        disturbance(:) = 0;
        return;
    end
    boundary = zeros(6,1);
    for sweep = 1:100
        linear_next = boundary;
        for phase = count:-1:1
            input = B(:,:,phase); P_next = P_schedule(:,:,phase);
            offset = P_next*disturbance(:,phase)+linear_next;
            feedforward(:,phase) = -(R+input.'*P_next*input)\(input.'*offset);
            linear_next = (A(:,:,phase)-input*K(:,:,phase)).'*offset;
        end
        residual = norm(linear_next-boundary)/max(1,norm(linear_next));
        boundary = linear_next;
        if residual<1e-10, break; end
    end
    assert(residual<1e-10,'Periodic disturbance feedforward did not converge.');
end

function next = rk4_step(t,state,dt,dipole,p)
    k1 = dynamics(t,state,dipole,p);
    k2 = dynamics(t+dt/2,state+dt*k1/2,dipole,p);
    k3 = dynamics(t+dt/2,state+dt*k2/2,dipole,p);
    k4 = dynamics(t+dt,state+dt*k3,dipole,p);
    next = state+dt*(k1+2*k2+2*k3+k4)/6;
end

function derivative = dynamics(t,state,dipole,p)
    q = state(1:4)/norm(state(1:4)); omega = state(5:7);
    C_bi = quat_to_dcm(q);
    [disturbance,B_body] = environment(t,C_bi,p);
    acceleration = p.J\(disturbance+cross(dipole,B_body) ...
        -cross(omega,p.J*omega));
    qdot = 0.5*[-dot(q(2:4),omega); q(1)*omega+cross(q(2:4),omega)];
    derivative = [qdot;acceleration];
end

function [disturbance,B_body] = environment(t,C_bi,p)
    [r,v] = orbit_state(t,p);
    B_body = C_bi*earth_dipole_B(r);
    nadir = C_bi*(-r/norm(r));
    gravity_gradient = 3*p.muE/norm(r)^3*cross(nadir,p.J*nadir);
    density = p.rho_ref*exp(-(norm(r)-p.Re-p.h_ref)/p.H_scale);
    relative_velocity = C_bi*(v-cross(p.omegaE,r));
    drag_force = -0.5*density*p.Cd*p.A_drag_ref ...
        *norm(relative_velocity)*relative_velocity;
    srp_force = -p.P_srp*p.Cr*p.A_srp_ref*(C_bi*p.s_eci_const);
    disturbance = gravity_gradient+cross(p.r_cp_drag_body,drag_force) ...
        +cross(p.r_cp_srp_body,srp_force)+cross(p.m_res_body,B_body);
end

function [r,v] = orbit_state(t,p)
    angle = p.n*t;
    r = p.R_oi*[p.a*cos(angle);p.a*sin(angle);0];
    v = p.R_oi*[-p.a*p.n*sin(angle);p.a*p.n*cos(angle);0];
end

function B = earth_dipole_B(r)
    er = r/norm(r); north = [0;0;1];
    B = 7.96e15/norm(r)^3*(3*dot(north,er)*er-north);
end

function S = skew(v)
    S = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
end

function theta = attitude_error(C_bi,C_ref)
    error_matrix = C_bi*C_ref.';
    angle = acos(max(-1,min(1,(trace(error_matrix)-1)/2)));
    antisymmetric = error_matrix.'-error_matrix;
    vector = [antisymmetric(3,2);antisymmetric(1,3);antisymmetric(2,1)];
    if angle<1e-7
        theta = vector/2;
    else
        assert(angle<pi-1e-4,'Attitude error is outside the local chart.');
        theta = angle/(2*sin(angle))*vector;
    end
end

function C = quat_to_dcm(q)
    scalar = q(1); vector = q(2:4);
    C = (scalar^2-dot(vector,vector))*eye(3) ...
        +2*(vector*vector.')-2*scalar*skew(vector);
end

function q = dcm_to_quat(C)
    % Standard active-matrix conversion is applied to C transpose, since C is passive.
    R = C.';
    if trace(R)>0
        scale = 2*sqrt(1+trace(R));
        q = [scale/4;(R(3,2)-R(2,3))/scale; ...
            (R(1,3)-R(3,1))/scale;(R(2,1)-R(1,2))/scale];
    else
        [~,axis_index] = max(diag(R));
        j = mod(axis_index,3)+1; k = mod(j,3)+1;
        scale = 2*sqrt(1+R(axis_index,axis_index)-R(j,j)-R(k,k));
        vector = zeros(3,1); vector(axis_index) = scale/4;
        vector(j) = (R(j,axis_index)+R(axis_index,j))/scale;
        vector(k) = (R(k,axis_index)+R(axis_index,k))/scale;
        q = [(R(k,j)-R(j,k))/scale;vector];
    end
    q = q/norm(q);
    if q(1)<0, q=-q; end
end
