clear; clc;

mu = 3.986e14; Re = 6371e3; alt = 510e3; r = Re + alt;
rho = 2e-12; Cd = 2.2; A = 1.097; m = 40;  
v_orbit = sqrt(mu/r);

a = 0.85; b = 0.5; c = 0.4;
Ixx = (1/12)*m*(b^2 + c^2);
Iyy = (1/12)*m*(a^2 + c^2);
Izz = (1/12)*m*(a^2 + b^2);
I = mean([Ixx Iyy Izz]);

wheel_mass = 0.21;
unbalance_static = 1.5e-6;
unbalance_dynamic = 20e-9;
r_payload = 0.3;

L = a/2;
a_gg = 3*mu/r^3 * L;          
a_drag = 0.5 * rho * Cd * A * v_orbit^2 / m;
v_internal = 0.001;           

tau_w = 0.005;       
tau_dist = 5e-4;     
bw = 0.05;           
omega_bw = 2*pi*bw;

tau_res = max(tau_dist - 0.8*tau_w, 0.1*tau_dist);
alpha_res = tau_res / I;
omega_res = alpha_res / omega_bw;

a_centrifugal = omega_res^2 * r_payload;
a_coriolis = 2 * omega_res * v_internal;

rpm_vals = 0:100:10000;
k_vals = [1e3, 3e3, 4e3, 3e4, 1e5];
zeta_vals = [0.01, 0.05, 0.1];
colors = lines(length(k_vals));

figure('Color','w');
tiledlayout(3,1);

for z = 1:length(zeta_vals)
    zeta = zeta_vals(z);
    MicroG = zeros(length(k_vals), length(rpm_vals));

    for i = 1:length(k_vals)
        k_damper = k_vals(i);
        omega_n = sqrt(k_damper / wheel_mass);

        for j = 1:length(rpm_vals)
            rpm = rpm_vals(j);
            omega_wheel = rpm * 2*pi/60;

            a_static = (unbalance_static * omega_wheel^2) / m;
            alpha_vib = omega_wheel;
            tau_dyn = unbalance_dynamic * alpha_vib;
            a_dynamic = tau_dyn / (m * r_payload);
            a_unbalance = sqrt(a_static^2 + a_dynamic^2);

            r = omega_wheel / omega_n;
            T_omega = sqrt( (1 + (2*zeta*r)^2) / ((1 - r^2)^2 + (2*zeta*r)^2) );
            a_unbalance_damped = a_unbalance * T_omega;

            a_total_rms = sqrt(a_gg^2 + a_drag^2 + ...
                               a_centrifugal^2 + a_coriolis^2 + a_unbalance_damped^2);
            MicroG(i,j) = a_total_rms / 9.81 * 1e6;
        end
    end

    nexttile; hold on;
    plotHandles = gobjects(length(k_vals),1);

    for i = 1:length(k_vals)
        plotHandles(i) = plot(rpm_vals, MicroG(i,:), 'LineWidth', 2, 'Color', colors(i,:), ...
                              'DisplayName', sprintf('k = %.1e N/m', k_vals(i)));
    end

    hLimit = yline(100, '--r', 'LineWidth', 1.2);
    hLimit.DisplayName = 'Payload Limit (100 µg)';

    xlabel('Reaction Wheel Speed (rpm)');
    ylabel('Microgravity (\mug)');
    title(sprintf('Damping Ratio ζ = %.2f', zeta));
    grid on;
    xlim([0 10000]);
    ylim([0 150]);

    if abs(zeta - 0.1) < 1e-3
        legend([plotHandles; hLimit], 'Location','northwest');
    end

    fprintf('\n--- ζ = %.2f ---\n', zeta);
    for i = 1:length(k_vals)
        fprintf('k = %.1e N/m → Micro-g range: %.3f – %.3f µg\n', ...
            k_vals(i), min(MicroG(i,:)), max(MicroG(i,:)));
    end
end

sgtitle('Microgravity with RW Speed Under Different Dynamic Conditions');
