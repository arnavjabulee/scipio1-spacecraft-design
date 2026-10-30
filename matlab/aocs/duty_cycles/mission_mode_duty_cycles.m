%% Duty cycles generator (edit variables below)

% ---- Inputs
h_km                 = 510;
mission_days         = 200;
leop_days            = 7;
safe_days            = 6;
cam_events           = 3;
cam_hours_each       = 4;
decom_days           = 4;

st_availability      = 0.95;   % fraction of time star tracker usable
leop_st_frac         = 1;   % fraction of LEOP with ST on
cam_st_frac          = 1;   % fraction of CAM with ST on

sun_sensor_eclipse_on= 0;      % 1=use SS in eclipse in ops
mtq_dipole_Am2       = 12;     % per-axis dipole capability
B_T                  = 35e-6;  % mean field [T]
T_disturb_Nm         = 5e-6;   % lumped disturbance torque
mtq_efficiency       = 0.4;    % bang-bang proxy
science_wheel_frac   = 0.10;   % one wheel assist fraction in Science

tot_h   = mission_days*24;
F.LEOP  = leop_days*24 / tot_h;
F.Safe  = safe_days*24 / tot_h;
F.Collision = cam_events*cam_hours_each / tot_h;
F.Decommission = decom_days*24 / tot_h;
F.Science = max(1 - (F.LEOP+F.Safe+F.Collision+F.Decommission), 0);

sunlit = sunlit_fraction_seasonal(h_km,181);
st_science = st_availability * sunlit;

D.CubeStar     = [leop_st_frac, st_science, 0.00, cam_st_frac, 0.00];
ss_science     = clamp(0.85 + 0.05*(sun_sensor_eclipse_on~=0), 0, 1);
D.CubeSense    = [0.80, ss_science, 1.00, 0.60, 0.70];
D.CubeMag      = [1.00, 1.00, 1.00, 1.00, 1.00];

T_mtq_max      = mtq_dipole_Am2 * B_T;
need_frac      = clamp(T_disturb_Nm / max(T_mtq_max,eps), 0, 1);
mtq_science    = clamp(mtq_efficiency * need_frac, 0.2, 0.6);
D.CR0100       = [0.75, mtq_science, 0.60, 0.30, 0.30];

D.CW0500       = [0.10, science_wheel_frac, 0.10, 0.90, 1.00];
D.CubeComputer = [1,1,1,1,1];

% ---- Output
disp('%% ---- Duty cycles (fraction of phase time active)');
disp('% Order: [LEOP, Science, Safe, Collision, Decommission]');
fprintf('D.CubeStar     = [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f];\n', D.CubeStar);
fprintf('D.CubeMag      = [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f];\n', D.CubeMag);
fprintf('D.CubeSense    = [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f];\n', D.CubeSense);
fprintf('D.CR0100       = [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f];\n', D.CR0100);
fprintf('D.CW0500       = [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f];\n', D.CW0500);
fprintf('D.CubeComputer = [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f];\n', D.CubeComputer);

% ===== helpers =====
function f_sunlit = sunlit_fraction_seasonal(h_km,N)
Re=6378.137e3; r=Re+h_km*1e3; psi=acos(Re/r);
beta=linspace(-23.5,23.5,N)*pi/180;
f_ecl=zeros(size(beta));
for i=1:numel(beta)
    b=beta(i); beta_crit=asin(Re/r);
    if abs(b)>=beta_crit, f_ecl(i)=0; else
        arg=cos(psi)/max(cos(b),1e-9); arg=min(max(arg,-1),1);
        f_ecl(i)=(1/pi)*acos(arg);
    end
end
f_sunlit = max(0, min(1, 1-mean(f_ecl)));
end
function y=clamp(x,a,b), y=min(max(x,a),b); end
