function report = benchmark_pd_lqr(seeds)
%BENCHMARK_PD_LQR Paired nominal Sun-pointing comparison on held-out noise.
% Original PD gains, same nonlinear plant, initial state, sensor noise, sample
% period, integration, two-orbit duration and magnetic actuator limits.
% Evaluate the second orbit. Default seeds 4:10 were not used for tuning.
% Execution copies and results are written under ignored outputs/.
arguments
    seeds (1,:) double {mustBeInteger,mustBeNonnegative} = 4:10
end
assert(~isempty(seeds) && numel(unique(seeds))==numel(seeds), ...
    'Supply distinct noise seeds.');
repository = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(repository,'matlab','aocs','nominal_control', ...
    'sun_tracking_lqr_kalman_simulation.m'));
pdSource = fileread(fullfile(repository,'matlab','aocs','nominal_control', ...
    'sun_tracking_pd_simulation.m'));
kpToken = regexp(pdSource,'Kp\s*=\s*([0-9eE.+-]+)','tokens','once');
kdToken = regexp(pdSource,'Kd\s*=\s*([0-9eE.+-]+)','tokens','once');
assert(~isempty(kpToken) && ~isempty(kdToken),'Cannot read original PD gains.');
Kp = str2double(kpToken{1}); Kd = str2double(kdToken{1});
plotStart = strfind(source,'%% Plots (');
localStart = strfind(source,'%% Local functions');
assert(isscalar(plotStart) && isscalar(localStart),'Cannot locate script blocks.');
source = [source(1:plotStart-1) source(localStart:end)];

% PD leaves roll unconstrained; its logging must work beyond the LQR's local
% attitude range. This quaternion chart is equivalent in the local LQR range.
chartStart = strfind(source,'function theta = attitude_error(');
chartEnd = strfind(source,'function C = quat_to_dcm(');
chart = sprintf([ ...
    'function theta = attitude_error(C_bi,C_ref)\n' ...
    '    dq = dcm_to_quat(C_bi*C_ref.'');\n' ...
    '    size_vector = norm(dq(2:4));\n' ...
    '    if size_vector<1e-12, theta=2*dq(2:4); else\n' ...
    '        theta=2*atan2(size_vector,dq(1))*dq(2:4)/size_vector;\n' ...
    '    end\nend\n\n']);
source = [source(1:chartStart-1) chart source(chartEnd:end)];
originalControl = 'm_requested = -K_lqr(:,:,phase)*x_posterior+m_feedforward(:,phase);';
assert(contains(source,originalControl),'Cannot locate the LQR feedback law.');
pdControl = sprintf([ ...
    'measured_attitude = expm(-skew(measurements(1:3,k)))*C_ref;\n' ...
    '    measured_sun = measured_attitude*p.s_eci_const;\n' ...
    '    measured_field = measured_attitude*C_ref.''*B_reference(:,phase);\n' ...
    '    pd_torque = -%.17g*cross(measured_sun,[1;0;0]) ...\n' ...
    '        -%.17g*measurements(4:6,k);\n' ...
    '    m_requested = cross(measured_field,pd_torque)/dot(measured_field,measured_field);'],Kp,Kd);

destination = fullfile(repository,'outputs','lqr_kalman','benchmark');
if ~isfolder(destination), mkdir(destination); end
warningState = warning('off','LQRKF:SmallAngleRange');
restoreWarnings = onCleanup(@() warning(warningState));
records = struct([]); traces = struct([]); index = 0;
controllers = {'pd','lqr'};
for controllerIndex = 1:numel(controllers)
    controller = controllers{controllerIndex};
    code = source;
    if strcmp(controller,'pd')
        % The KF/design code is not in PD's feedback path. PD receives the
        % raw noisy attitude/rate measurements, its original cross-product
        % sun error, original gains and original magnetic dipole allocation.
        code = strrep(code,originalControl,pdControl);
    end
    for seed = seeds
        executed = strrep(code,"rng(1,'twister');", ...
            sprintf("rng(%d,'twister');",seed));
        assert(~strcmp(executed,code) || seed==1,'Cannot set benchmark seed.');
        path = fullfile(destination,sprintf('%s_seed_%d.m',controller,seed));
        fileId = fopen(path,'w','n','UTF-8');
        assert(fileId>=0,'Cannot write benchmark execution copy.');
        fwrite(fileId,executed,'char'); fclose(fileId);
        fprintf('Paired benchmark: %s, seed %d\n',controller,seed);
        evalin('base',"run('"+replace(string(path),"'","''")+"');");
        data = evalin('base',[ ...
            'struct(''t'',t,''period'',Torb,''sun'',sun_err_deg,' ...
            '''dipole'',m_cmd_body,''saturation'',saturation,' ...
            '''q_initial'',q_hist(:,1),''w_initial'',w_hist(:,1),' ...
            '''max_attitude'',metrics.max_attitude_error_deg)']);
        mask = data.t>=data.period;
        times = data.t(mask); duration = times(end)-times(1);
        record.controller = controller; record.seed = seed;
        record.rms_sun_deg = sqrt(trapz(times,data.sun(mask).^2)/duration);
        record.peak_sun_deg = max(data.sun(mask));
        record.rms_dipole_Am2 = sqrt(trapz(times, ...
            sum(data.dipole(:,mask).^2,1))/duration);
        record.saturated_percent = 100*mean(any(data.saturation(:,mask),1));
        record.max_attitude_deg = data.max_attitude;
        record.initial_state = [data.q_initial;data.w_initial].';
        index = index+1;
        if index==1, records=record; else, records(index)=record; end
        if seed==seeds(1)
            trace.controller = controller;
            trace.time_minutes = data.t/60; trace.sun_error_deg = data.sun;
            if controllerIndex==1, traces=trace; else, traces(controllerIndex)=trace; end
        end
    end
end
pd = records(strcmp({records.controller},'pd'));
lqr = records(strcmp({records.controller},'lqr'));
assert(norm(vertcat(pd.initial_state)-vertcat(lqr.initial_state),'fro')<1e-12, ...
    'Paired tests used different physical initial states.');
report.seeds = seeds;
report.pd_gains = [Kp Kd]; report.matlabVersion = version;
report.records = records;
for controllerIndex = 1:numel(controllers)
    name = controllers{controllerIndex};
    rows = records(strcmp({records.controller},name));
    report.summary.(name).rms_sun_deg = mean([rows.rms_sun_deg]);
    report.summary.(name).peak_sun_deg = mean([rows.peak_sun_deg]);
    report.summary.(name).rms_dipole_Am2 = mean([rows.rms_dipole_Am2]);
end
report.rms_improvement_percent = 100*(1-report.summary.lqr.rms_sun_deg ...
    /report.summary.pd.rms_sun_deg);
report.peak_improvement_percent = 100*(1-report.summary.lqr.peak_sun_deg ...
    /report.summary.pd.peak_sun_deg);
report.dipole_ratio = report.summary.lqr.rms_dipole_Am2 ...
    /report.summary.pd.rms_dipole_Am2;
report.passed = all([lqr.rms_sun_deg]<[pd.rms_sun_deg]) ...
    && all([lqr.peak_sun_deg]<[pd.peak_sun_deg]) ...
    && all([lqr.rms_dipole_Am2]<=1.1*[pd.rms_dipole_Am2]) ...
    && all([lqr.max_attitude_deg]<10);
fileId = fopen(fullfile(destination,'benchmark.json'),'w','n','UTF-8');
assert(fileId>=0,'Cannot write benchmark summary.');
fwrite(fileId,jsonencode(report,PrettyPrint=true),'char'); fclose(fileId);
save(fullfile(destination,'traces.mat'),'traces');
fprintf('Held-out RMS improvement %.2f%%; peak improvement %.2f%%; dipole ratio %.4f\n', ...
    report.rms_improvement_percent,report.peak_improvement_percent,report.dipole_ratio);
assert(report.passed, ...
    'Tuned LQR must beat PD RMS/peak per seed, use <=110%% effort and stay below 10 deg.');
end
