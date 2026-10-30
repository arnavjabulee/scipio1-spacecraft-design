%% MTQ-FIRST ACTUATION WITH RW BACKUP (NOMINAL MODE, 510 km SSO)

clear; clc;

%% ORBIT / ENVIRONMENT
altitude  = 510e3;
Re        = 6371e3;
mu        = 3.986e14;
B_eq      = 30e-6;
r         = Re + altitude;
B_mag     = B_eq * (Re/r)^3;

n_orbit   = sqrt(mu/r^3);
T_orbit   = 2*pi/n_orbit;

%% INERTIA MATRIX
J = [  1.991   -0.04488  -0.03172 ; ...
      -0.04488  2.721     0.03037 ; ...
      -0.03172  0.03037   1.329  ];

%% ACTUATORS
m_max = [15;15;15];         % MTQ [A·m^2]
rw_peak = 20e-3;            % [N·m]
rw_nom  = 0.5*rw_peak;

P_MTQ = [0.5 1.25 0.1];     % [nom peak idle]
P_RW  = [3.5 32.0 0.1];

omega_thr = deg2rad(0.1);
rw_deadband = 0.1 * rw_nom; % only use wheels above this residual

%% SIMULATION GRID
N_orbits = 10;
dt = 0.1;
N  = round(N_orbits*T_orbit/dt);

omega = zeros(3,N);
Tdist = zeros(3,N);

%% DISTURBANCES
Ivals = eig(J);
DeltaI = max(Ivals)-min(Ivals);
Tgg = 3*mu/(r^3)*DeltaI * 0.3;

Tmisc = 5e-8;

for k = 1:N
    t = (k-1)*dt;
    Tdist(:,k) = ...
        [0.3*sin(2*pi*t/T_orbit);
         0.7*sin(2*pi*t/T_orbit+pi/4);
         0.5*sin(4*pi*t/T_orbit+pi/3)] * Tgg ...
      + Tmisc*[sin(2*pi*t/0.6/T_orbit);
               sin(2*pi*t/0.8/T_orbit);
               sin(2*pi*t/1.2/T_orbit)];
end

%% DUTY COUNTERS
mtq_nom = 0; mtq_peak = 0; mtq_idle = 0;
rw_nom  = 0; rw_peak  = 0; rw_idle  = 0;

%% MAIN LOOP
for k = 1:N-1
    t  = (k-1)*dt;
    w  = omega(:,k);
    Td = Tdist(:,k);

    B = B_mag * [cos(n_orbit*t); sin(n_orbit*t); 0];
    S = [0 -B(3) B(2); B(3) 0 -B(1); -B(2) B(1) 0];

    %% MTQ ALLOCATION
    m_cmd = zeros(3,1);
    Tmtq  = zeros(3,1);

    if norm(w) <= omega_thr
        m_unc = pinv(S)*(-Td);
        for i=1:3
            m_cmd(i) = max(min(m_unc(i), m_max(i)), -m_max(i));
        end
        Tmtq = S*m_cmd;
        Tres = Td + Tmtq;
    else
        Tres = Td;
    end

    %% RW ALLOCATION WITH DEADBAND
    Trw = zeros(3,1);
    if norm(Tres) > rw_deadband || norm(w) > omega_thr
        for i=1:3
            Treq = -Tres(i);
            if abs(Treq) <= rw_nom
                Trw(i)=Treq;
            elseif abs(Treq)<= rw_peak
                Trw(i)=Treq;
            else
                Trw(i)=sign(Treq)*rw_peak;
            end
        end
    end

    %% DUTY STATE: MTQs
    if norm(m_cmd) < 1e-9
        mtq_idle = mtq_idle + dt;
    else
        if norm(Tmtq) <= 0.5*(norm(m_max)*B_mag)
            mtq_nom = mtq_nom + dt;
        else
            mtq_peak = mtq_peak + dt;
        end
    end

    %% DUTY STATE: RWs
    if norm(Trw) < 1e-12
        rw_idle = rw_idle + dt;
    else
        if max(abs(Trw)) <= 0.5*20e-3
            rw_nom = rw_nom + dt;
        else
            rw_peak = rw_peak + dt;
        end
    end

    %% DYNAMICS
    Ttot = Td + Tmtq + Trw;
    omega(:,k+1) = omega(:,k) + dt*(J\(Ttot - cross(w,J*w)));
end

total = N*dt;

%% FRACTIONS
mtq_f = [mtq_nom mtq_peak mtq_idle]/total;
rw_f  = [rw_nom  rw_peak  rw_idle]/total;

%% POWER
P_mtq = mtq_f(1)*P_MTQ(1) + mtq_f(2)*P_MTQ(2) + mtq_f(3)*P_MTQ(3);
P_rw  = rw_f(1)*P_RW(1)  + rw_f(2)*P_RW(2)  + rw_f(3)*P_RW(3);

fprintf("\n--- NOMINAL MODE DUTY CYCLES ---\n")
fprintf("MTQ  [nom peak idle] = [%5.3f  %5.3f  %5.3f]\n", mtq_f)
fprintf("RW   [nom peak idle] = [%5.3f  %5.3f  %5.3f]\n", rw_f)
fprintf("Total MTQ power = %.3f W (3x)\n", 3*P_mtq)
fprintf("Total RW power  = %.3f W (4x)\n", 4*P_rw)
