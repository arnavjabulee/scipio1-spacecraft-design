clear; clc;

%% BASIC ORBIT / MASS PROPS
mu = 3.986e14; Re = 6371e3; alt = 510e3; r = Re + alt;
rho = 2e-12; Cd = 2.2; A = 1.097; m = 40;
B = 3e-5; 
v_orbit = sqrt(mu/r);

a = 0.85; b = 0.5; c = 0.4;
Ixx = (1/12)*m*(b^2 + c^2);
Iyy = (1/12)*m*(a^2 + c^2);
Izz = (1/12)*m*(a^2 + b^2);
I = mean([Ixx Iyy Izz]); 

tau_dist   = [1e-4 5e-4 1e-3]; % use 1e-5 and 1e-6, higher values for aero drag. re-entry ~1e-3
bandwidths = [0.001 0.01 0.05];
dipoles    = [5 10 20 40];

L          = a/2;
a_gg       = 3*mu/r^3 * L;
a_drag     = 0.5 * rho * Cd * A * v_orbit^2 / m; 
v_internal = 0.001;
r_payload  = 0.3;

% %% MICRO-G ANALYSIS (unchanged)
% results = [];
% 
% for m_dip = dipoles
%     tau_mag_max = m_dip * B; 
%     for bw = bandwidths
%         omega_bw = 2*pi*bw;
%         for td = tau_dist
%             tau_res      = max(td - 0.5*tau_mag_max, 0.1*td);
%             alpha_res    = tau_res / I;
%             omega_res    = alpha_res / omega_bw;
%             a_centrifugal = omega_res^2 * r_payload;
%             a_coriolis     = 2 * omega_res * v_internal;
% 
%             a_total_rms = sqrt(a_gg^2 + a_drag^2 + ...
%                                a_centrifugal^2 + a_coriolis^2);
%             micro_g = a_total_rms / 9.81;
% 
%             results = [results; m_dip bw td omega_res*180/pi a_total_rms micro_g];
%         end
%     end
% end
% 
% colNames = {'Dipole_A_m2','Bandwidth_Hz','DistTorque_Nm', ...
%             'ResidualRate_deg_s','TotalAccel_m_s2','MicroG'};
% T = array2table(results, 'VariableNames', colNames);
% 
% disp(T);
% 
% figure('Color','w');
% 
% subplot(2,1,1);
% hold on;
% unique_bw = unique(T.Bandwidth_Hz);
% colors = lines(length(unique_bw));
% 
% for i = 1:length(unique_bw)
%     bw = unique_bw(i);
%     idx = T.Bandwidth_Hz == bw;
%     plot(T.Dipole_A_m2(idx), T.MicroG(idx)*1e6, '-o', ...
%          'Color', colors(i,:), 'LineWidth', 1.5, ...
%          'DisplayName', sprintf('BW = %.3f Hz', bw));
% end
% 
% yline(100, '--r', 'Payload Limit (100 μg)', 'LineWidth', 1.2);
% 
% xlabel('Dipole Moment (A·m²)');
% ylabel('Microgravity (μg)');
% title('Microgravity vs Dipole Moment');
% legend('Location','best');
% grid on;
% ylim([0 150]);
% 
% 
% subplot(2,1,2);
% hold on;
% unique_dip = unique(T.Dipole_A_m2);
% colors = lines(length(unique_dip));
% 
% for i = 1:length(unique_dip)
%     dip = unique_dip(i);
%     idx = T.Dipole_A_m2 == dip;
%     plot(T.Bandwidth_Hz(idx), T.MicroG(idx)*1e6, '-s', ...
%          'Color', colors(i,:), 'LineWidth', 1.5, ...
%          'DisplayName', sprintf('Dipole = %.0f A·m²', dip));
% end
% 
% yline(100, '--r', 'Payload Limit (100 μg)', 'LineWidth', 1.2);
% 
% xlabel('Control Bandwidth (Hz)');
% ylabel('Microgravity (μg)');
% title('Microgravity vs Control Bandwidth');
% legend('Location','best');
% grid on;
% ylim([0 150])

%% PER-AXIS MAGNETORQUER CONTROLLABILITY

incl     = deg2rad(97.5);
tilt_mag = deg2rad(11);
n_points = 500;
theta_orbit = linspace(0, 2*pi, n_points);

CBx = zeros(1, n_points);   % controllability index about X
CBy = zeros(1, n_points);   % controllability index about Y
CBz = zeros(1, n_points);   % controllability index about Z

m_eci = [sin(tilt_mag); 0; cos(tilt_mag)];
R_body_to_ECI = eye(3);  % body fixed to ECI (global-Z pointing)

for i = 1:n_points
    u = theta_orbit(i);
    r_orb = [cos(u); ...
             sin(u)*cos(incl); ...
             sin(u)*sin(incl)];
    
    Bvec_eci  = 3*(dot(m_eci,r_orb))*r_orb - m_eci;
    Bvec_body = R_body_to_ECI' * Bvec_eci;
    
    Bnorm = norm(Bvec_body);
    Bx = Bvec_body(1) / Bnorm;
    By = Bvec_body(2) / Bnorm;
    Bz = Bvec_body(3) / Bnorm;
    
    % theta_i = angle between axis i and B, so cos(theta_i) = |component|
    thetaBx = acos(abs(Bx));
    thetaBy = acos(abs(By));
    thetaBz = acos(abs(Bz));
    
    % Controllability index: sin(theta_i) (0 = bad, 1 = best)
    CBx(i) = sin(thetaBx);
    CBy(i) = sin(thetaBy);
    CBz(i) = sin(thetaBz);
end

theta_deg = rad2deg(theta_orbit);

figure('Color','w');

subplot(3,1,1);
plot(theta_deg, CBx, 'LineWidth', 1.5);
xlabel('Orbital True Anomaly (deg)');
ylabel('C_{Bx} = sin(\theta_{Bx})');
title('Magnetorquer Controllability about X-axis');
grid on;
hold on;
yline(0.8, '--g', 'Good');
yline(0.5, '--', 'Moderate', 'Color', [1, 0.5, 0]);
yline(0.2, '--r', 'Poor');
hold off;

subplot(3,1,2);
plot(theta_deg, CBy, 'LineWidth', 1.5);
xlabel('Orbital True Anomaly (deg)');
ylabel('C_{By} = sin(\theta_{By})');
title('Magnetorquer Controllability about Y-axis');
grid on;
hold on;
yline(0.8, '--g', 'Good');
yline(0.5, '--', 'Moderate', 'Color', [1, 0.5, 0]);
yline(0.2, '--r', 'Poor');
hold off;

subplot(3,1,3);
plot(theta_deg, CBz, 'LineWidth', 1.5);
xlabel('Orbital True Anomaly (deg)');
ylabel('C_{Bz} = sin(\theta_{Bz})');
title('Magnetorquer Controllability about Z-axis');
grid on;
hold on;
yline(0.8, '--g', 'Good');
yline(0.5, '--', 'Moderate', 'Color', [1, 0.5, 0]);
yline(0.2, '--r', 'Poor');
hold off;

fprintf('Average controllability index (X): %.3f, min: %.3f\n', mean(CBx), min(CBx));
fprintf('Average controllability index (Y): %.3f, min: %.3f\n', mean(CBy), min(CBy));
fprintf('Average controllability index (Z): %.3f, min: %.3f\n', mean(CBz), min(CBz));

% %% 3D VISUALISATION (unchanged except variable reuse)
% 
% figure('Color','w');
% hold on;
% axis equal;
% grid on;
% xlabel('X (ECI)');
% ylabel('Y (ECI)');
% zlabel('Z (ECI)');
% title('Magnetic Field & Satellite Orientation');
% view(35,25);
% 
% [Xe,Ye,Ze] = sphere(50);
% surf(Re*Xe, Re*Ye, Re*Ze, ...
%     'FaceAlpha', 0.25, 'EdgeColor', 'none', 'FaceColor', [0.4 0.6 1]);
% 
% incl     = deg2rad(97.5);
% tilt_mag = deg2rad(11);
% n_points = 100;
% theta_orbit = linspace(0, 2*pi, n_points);
% 
% m_eci = [sin(tilt_mag); 0; cos(tilt_mag)]; 
% r_orbit = zeros(3,n_points);
% B_dir   = zeros(3,n_points);
% 
% for i = 1:n_points
%     u = theta_orbit(i);
%     r_orb = [cos(u); sin(u)*cos(incl); sin(u)*sin(incl)];
%     r_orbit(:,i) = r_orb;
% 
%     Bvec     = 3*(dot(m_eci,r_orb))*r_orb - m_eci;
%     B_dir(:,i) = Bvec / norm(Bvec);
% end
% 
% plot3(r_orbit(1,:)*Re, r_orbit(2,:)*Re, r_orbit(3,:)*Re, ...
%     'Color','k','LineWidth',1.2,'DisplayName','Orbit Path');
% 
% skip = 8; 
% scale = 0.4*Re;
% for i = 1:skip:n_points
%     quiver3(r_orbit(1,i)*Re, r_orbit(2,i)*Re, r_orbit(3,i)*Re, ...
%             B_dir(1,i)*scale, B_dir(2,i)*scale, B_dir(3,i)*scale, ...
%             'Color',[0 0.6 0],'LineWidth',1.5,'MaxHeadSize',0.9, ...
%             'DisplayName','Magnetic Field Vectors');
% end
% 
% h_orbit = plot3(r_orbit(1,:)*Re, r_orbit(2,:)*Re, r_orbit(3,:)*Re, ...
%     'Color','k','LineWidth',1.2);
% 
% for i = 1:skip:n_points
%     quiver3(r_orbit(1,i)*Re, r_orbit(2,i)*Re, r_orbit(3,i)*Re, ...
%             B_dir(1,i)*scale, B_dir(2,i)*scale, B_dir(3,i)*scale, ...
%             'Color',[0 0.6 0],'LineWidth',1.5,'MaxHeadSize',0.9);
% end
% h_B = quiver3(r_orbit(1,1)*Re, r_orbit(2,1)*Re, r_orbit(3,1)*Re, ...
%             B_dir(1,1)*scale, B_dir(2,1)*scale, B_dir(3,1)*scale, ...
%             'Color',[0 0.6 0],'LineWidth',1.5,'MaxHeadSize',0.9);
% 
% h_Z = quiver3(0,0,0,0,0,Re*1.5, ...
%     'Color','r','LineWidth',2,'MaxHeadSize',0.5);
% text(0,0,Re*1.55,'Body Z-axis','Color','r','FontWeight','bold');
% 
% h_dip = quiver3(0,0,0,m_eci(1)*Re*1.5,m_eci(2)*Re*1.5,m_eci(3)*Re*1.5, ...
%     'Color','b','LineWidth',2,'MaxHeadSize',0.5);
% text(m_eci(1)*Re*1.55,m_eci(2)*Re*1.55,m_eci(3)*Re*1.55,'Geomagnetic Axis', ...
%     'Color','b','FontWeight','bold');
% 
% h_earth = surf(Re*Xe, Re*Ye, Re*Ze, ...
%     'FaceAlpha', 0.25, 'EdgeColor', 'none', 'FaceColor', [0.4 0.6 1]);
% 
% legend([h_earth, h_orbit, h_B, h_Z, h_dip], ...
%     {'Earth','Orbit Path','Magnetic Field Vectors','Body Z-axis','Geomagnetic Axis'}, ...
%     'Location','northoutside');
