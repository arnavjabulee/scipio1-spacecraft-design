function report = validate_lqr_kalman
%VALIDATE_LQR_KALMAN Exercise the nonlinear LQR/KF example and its limits.
% Run from a fresh MATLAB session: these scripts clear the base workspace.
% Generated execution copies, metrics and plots go to ignored outputs/.
repository = fileparts(fileparts(mfilename('fullpath')));
sourcePath = fullfile(repository,'matlab','aocs','nominal_control', ...
    'sun_tracking_lqr_kalman_simulation.m');
source = fileread(sourcePath);
issues = checkcode(sourcePath,'-id');
assert(isempty(issues),'The LQR/KF script has Code Analyzer findings.');
destination = fullfile(repository,'outputs','lqr_kalman','validation');
if ~isfolder(destination), mkdir(destination); end
previousVisibility = get(groot,'DefaultFigureVisible');
restoreVisibility = onCleanup(@() set(groot,'DefaultFigureVisible',previousVisibility));
set(groot,'DefaultFigureVisible','off');

% Default nonlinear plant, noisy sensors and hard actuator saturation.
report.nominal = execute_case(source,'nominal',destination,true);
assert(report.nominal.max_attitude_error_deg<10, ...
    'The default example left the intended small-angle range.');
assert(all(report.nominal.attitude_estimation_rmse_deg < ...
    report.nominal.raw_attitude_measurement_rmse_deg), ...
    'Attitude filtering did not improve on the noisy measurements.');
assert(all(report.nominal.rate_estimation_rmse_rad_s < ...
    report.nominal.raw_rate_measurement_rmse_rad_s), ...
    'Rate filtering did not improve on the noisy measurements.');

% Isolate control signs, the quaternion reference, and convergence: remove
% truth disturbances and measurement noise, retain the initial perturbation.
noiseless = strrep(source,'measurement_std.*randn(6,1)','zeros(6,1)');
endOfTorque = '+cross(p.r_cp_srp_body,srp_force)+cross(p.m_res_body,B_body);';
assert(contains(noiseless,endOfTorque),'Disturbance test hook not found.');
unforced = strrep(noiseless,endOfTorque, ...
    [endOfTorque newline '    disturbance(:) = 0;']);
report.unforced = execute_case(unforced,'unforced',destination,false);
assert(report.unforced.final_sun_error_deg<0.01, ...
    'The undisturbed, noiseless perturbation did not converge.');
assert(report.unforced.last_orbit_rms_sun_error_deg<0.02, ...
    'The undisturbed last orbit did not settle.');

% Exact reference equilibrium should remain stationary without disturbances.
equilibrium = strrep(unforced,'theta_initial = deg2rad([1;-1.5;0.5]);', ...
    'theta_initial = zeros(3,1);');
equilibrium = strrep(equilibrium, ...
    'omega_initial = deg2rad([0.005;-0.005;0.0025]);', ...
    'omega_initial = zeros(3,1);');
report.equilibrium = execute_case(equilibrium,'equilibrium',destination,false);
assert(report.equilibrium.max_attitude_error_deg<0.001, ...
    'Reference equilibrium drifted; check kinematic/control conventions.');

% Halving the sample/integration step should give a nearby disturbed orbit.
% Use noiseless data to separate discretization changes from noise realizations.
report.step_two_seconds = execute_case(noiseless,'step_two',destination,false);
finer = strrep(noiseless,'sample_time_requested = 2;', ...
    'sample_time_requested = 1;');
report.step_one_second = execute_case(finer,'step_one',destination,false);
report.rms_step_difference_deg = abs( ...
    report.step_one_second.last_orbit_rms_sun_error_deg ...
    -report.step_two_seconds.last_orbit_rms_sun_error_deg);
assert(report.rms_step_difference_deg<0.1, ...
    'Last-orbit RMS changed by more than 0.1 degree with half the step.');
report.matlabVersion = version;
fileId = fopen(fullfile(destination,'validation.json'),'w','n','UTF-8');
assert(fileId>=0,'Cannot write validation summary.');
fwrite(fileId,jsonencode(report,PrettyPrint=true),'char'); fclose(fileId);
fprintf('LQR/KF validation passed: nominal, unforced, equilibrium, half-step.\n');
fprintf('Half-step last-orbit RMS difference: %.6f deg\n', ...
    report.rms_step_difference_deg);
end

function metrics = execute_case(source,name,destination,exportFigures)
path = fullfile(destination,[name '.m']);
fileId = fopen(path,'w','n','UTF-8');
assert(fileId>=0,'Cannot write validation execution copy.');
fwrite(fileId,source,'char'); fclose(fileId);
fprintf('Validation case: %s\n',name);
evalin('base',"run('"+replace(string(path),"'","''")+"');");
metrics = evalin('base','metrics');
if exportFigures
    figures = findall(groot,'Type','figure');
    [~,order] = sort([figures.Number]); figures = figures(order);
    for index = 1:numel(figures)
        exportgraphics(figures(index), ...
            fullfile(destination,sprintf('nominal_%02d.png',index)), ...
            'Resolution',150);
    end
end
end
