function manifest = run_analysis(mode, outputDirectory)
%RUN_ANALYSIS Run the released analyses and record reproducible outputs.
%   run_analysis('smoke') checks all scripts with two Monte Carlo runs.
%   run_analysis('full') uses each script's original run count.
%   Results, executed source, RNG state, MATLAB/toolbox versions and figures
%   are saved under outputs/matlab/<mode>. Run from a fresh MATLAB session:
%   these original scripts clear the base workspace and close figures.
arguments
    mode (1,1) string {mustBeMember(mode,["smoke","full"])} = "smoke"
    outputDirectory (1,1) string = ""
end
repository = fileparts(fileparts(mfilename('fullpath')));
if outputDirectory == ""
    outputDirectory = fullfile(repository,'outputs','matlab',mode);
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
outputDirectory = string(java.io.File(outputDirectory).getCanonicalPath());
sourceFiles = dir(fullfile(repository,'matlab','**','*.m'));
sourceFiles = sourceFiles(~strcmp({sourceFiles.name},'run_analysis.m'));
manifest.mode = mode;
manifest.generatedUtc = char(datetime('now','TimeZone','UTC'));
manifest.matlabVersion = version;
manifest.products = ver;
manifest.runs = struct([]);
previousVisibility = get(groot,'DefaultFigureVisible');
restoreVisibility = onCleanup(@() set(groot,'DefaultFigureVisible',previousVisibility));
set(groot,'DefaultFigureVisible','off');
for index = 1:numel(sourceFiles)
    original = fullfile(sourceFiles(index).folder,sourceFiles(index).name);
    [~,name] = fileparts(original);
    destination = fullfile(outputDirectory,name);
    if ~isfolder(destination), mkdir(destination); end
    source = fileread(original);
    executed = source;
    if mode == "smoke"
        executed = regexprep(executed,'(?m)^N_runs\s*=\s*\d+\s*;', 'N_runs = 2;');
    end
    executable = fullfile(destination,[name '.m']);
    fileId = fopen(executable,'w','n','UTF-8');
    assert(fileId >= 0,'Unable to write executed source.');
    fwrite(fileId,executed,'char'); fclose(fileId);
    record.source = erase(string(original),string(repository)+filesep);
    record.executedSource = executable;
    record.status = "passed";
    record.error = "";
    record.figures = strings(0,1);
    record.codeAnalyzer = checkcode(original,'-id');
    record.sourceSha256 = sha256(source);
    record.executedSha256 = sha256(executed);
    rng(1,'twister');
    record.initialRng = rng;
    started = tic;
    fprintf('Checking %s (%s)\n',name,mode);
    try
        evalin('base',"run('"+replace(string(executable),"'","''")+"');");
        handles = findall(groot,'Type','figure');
        if ~isempty(handles)
            [~,order] = sort([handles.Number]); handles = handles(order);
        end
        for figureIndex = 1:numel(handles)
            filename = fullfile(destination,sprintf('figure_%02d.png',figureIndex));
            exportgraphics(handles(figureIndex),filename,'Resolution',150);
            record.figures(end+1,1) = filename;
        end
        % Keep scalar/array values and parameter structures, excluding
        % function handles and graphics, in a workspace file for inspection.
        evalin('base',"save('"+replace(fullfile(destination,'workspace.mat'),"'","''")+"');");
    catch exception
        record.status = "failed";
        record.error = string(getReport(exception,'extended','hyperlinks','off'));
        fprintf(2,'%s\n',record.error);
    end
    record.elapsedSeconds = toc(started);
    if index == 1
        manifest.runs = record;
    else
        manifest.runs(index) = record;
    end
    close all;
    fileId = fopen(fullfile(outputDirectory,'manifest.json'),'w','n','UTF-8');
    fwrite(fileId,jsonencode(manifest,PrettyPrint=true),'char'); fclose(fileId);
end
fprintf('%d/%d scripts passed. Manifest: %s\n', ...
    sum([manifest.runs.status] == "passed"),numel(manifest.runs), ...
    fullfile(outputDirectory,'manifest.json'));
end

function digest = sha256(value)
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(unicode2native(value,'UTF-8'));
bytes = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(bytes,2).',1,[]));
end
