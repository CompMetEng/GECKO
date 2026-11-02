function report = upgradeGECKOproject(projectPath, varargin)
% upgradeGECKOproject  Bring a GECKO 3 project to the GECKO 4 layout.
%
% Lists, and with 'dryRun', false carries out, the steps that bring a
% project made with GECKO 3 (or an earlier GECKO 4) to the current layout:
%
% 1. adapter: an adapter class without model_adapter.toml gets one
%    (convertAdapterToToml). Methods the class overrides are listed: they
%    need carrying over by hand. The class itself stays.
% 2. custom kcats: data/customKcats.tsv becomes data/kcats/custom.tsv
%    (convertCustomKcats).
% 3. stored predictions: data/DLKcat.tsv and stored OpenKineticsPredictor
%    results become kcat source files (kcatsFromDLKcat, kcatsFromOKP with
%    'useStored'), without running anything.
% Steps 2 and 3 need the full ecModel the custom kcats and predictions are
% for ('ecModel'); without it they are listed and skipped.
% 4. saved ecModels: models/*.yml ecModels made with GECKO 3 get their
%    protein usage and pool reactions turned forward (upgradeEcModel).
% 5. project_format = 4 is written to model_adapter.toml.
% 6. scripts: lines of the project's .m files that call renamed or replaced
%    functions are listed, not changed.
%
% A file that is replaced or converted is first moved to
% legacy/<yyyy-mm-dd>/ in the project, keeping its place in the folder
% tree. Steps that are already done are left out, so the function can run
% again. geckopy has the same as `geckopy upgrade` (upgrade_project).
%
% Parameters
% ----------
% projectPath : char
%     the project folder.
%
% Name-Value Arguments
% --------------------
% dryRun : logical
%     only list the steps (default true).
% ecModel : struct
%     the ecModel the stored predictions were made for (step 3).
%
% Returns
% -------
% report : table
%     one row per step: step, status ('to do', 'done', 'skipped', 'failed',
%     'by hand') and note; script lines are listed one per row.
%
% Examples
% --------
%     upgradeGECKOproject(projectPath)                    % see what it would do
%     upgradeGECKOproject(projectPath, 'dryRun', false, 'ecModel', ecModel)
%
% See also
% --------
% convertAdapterToToml, convertCustomKcats, upgradeEcModel, createGECKOproject

checkOptions(varargin, {'dryRun','ecModel'}, 'upgradeGECKOproject');
p = parseGECKOargs(varargin, {'dryRun', true; 'ecModel', []});
projectPath = char(projectPath);
if ~isfolder(projectPath)
    error('upgradeGECKOproject:path', 'No folder %s.', projectPath)
end
apply = ~p.dryRun;
legacy = fullfile(projectPath, 'legacy', char(datetime('now', 'Format', 'yyyy-MM-dd')));
toml = fullfile(projectPath, 'model_adapter.toml');
rows = cell(0, 3);
    function add(step, status, note)
        rows(end+1, :) = {step, status, note};
    end
    function keepOriginal(file)
        % Move a file to legacy/<date>/, at its place in the folder tree.
        rel = file(numel(projectPath)+2:end);
        target = fullfile(legacy, rel);
        if ~isfolder(fileparts(target))
            mkdir(fileparts(target));
        end
        movefile(file, target);
    end
status = 'to do';
if apply
    status = 'done';
end

% 1. Adapter class -> model_adapter.toml
if ~isfile(toml)
    classFile = adapterClass(projectPath);
    if isempty(classFile)
        add('adapter', 'by hand', 'no model_adapter.toml and no adapter class; see createGECKOproject');
    else
        note = sprintf('%s -> model_adapter.toml', relative(classFile, projectPath));
        if apply
            r = convertAdapterToToml(classFile, toml);
            if ~isempty(r.methods)
                add('adapter methods', 'by hand', ['carry over: ' strjoin(r.methods, ', ')]);
            end
        end
        add('adapter', status, note);
    end
end

% 2. customKcats.tsv -> data/kcats/custom.tsv
custom = fullfile(projectPath, 'data', 'customKcats.tsv');
customOut = fullfile(projectPath, 'data', 'kcats', 'custom.tsv');
if isfile(custom) && ~isfile(customOut)
    note = 'data/customKcats.tsv -> data/kcats/custom.tsv';
    if isempty(p.ecModel)
        add('custom kcats', 'skipped', [note '; give the ecModel they are for (''ecModel'')']);
    else
        if apply
            [~, r] = convertCustomKcats(p.ecModel, 'customKcats', custom, 'outFile', customOut, ...
                'modelAdapter', projectAdapter(projectPath));
            keepOriginal(custom);
            note = sprintf('%s (%d rows)', note, r.rows);
        end
        add('custom kcats', status, note);
    end
end

% 3. Stored predictions -> kcat source files
stored = {};
if isfile(fullfile(projectPath, 'data', 'DLKcat.tsv')) && ~isfile(fullfile(projectPath, 'data', 'kcats', 'dlkcat.tsv'))
    stored{end+1} = 'dlkcat';
end
okpFiles = dir(fullfile(projectPath, 'data', 'OKP_output*.csv'));
for k = 1:numel(okpFiles)
    method = okpMethod(projectPath, okpFiles(k).name);
    if ~isempty(method)
        name = ['okp_' lower(regexprep(regexprep(method, '[^A-Za-z0-9]+', '_'), '^_+|_+$', ''))];
        if ~isfile(fullfile(projectPath, 'data', 'kcats', [name '.tsv']))
            stored{end+1} = method; %#ok<AGROW>
        end
    end
end
stored = unique(stored, 'stable');
if ~isempty(stored)
    note = ['stored predictions -> data/kcats: ' strjoin(stored, ', ')];
    if isempty(p.ecModel)
        add('predictions', 'skipped', [note '; give the ecModel they were made for (''ecModel'')']);
    else
        if apply
            adapter = projectAdapter(projectPath);
            for k = 1:numel(stored)
                try
                    if strcmp(stored{k}, 'dlkcat')
                        kcatsFromDLKcat(p.ecModel, 'modelAdapter', adapter);
                    else
                        kcatsFromOKP(p.ecModel, 'method', stored{k}, 'useStored', true, 'modelAdapter', adapter);
                    end
                    add('predictions', status, ['stored predictions -> data/kcats: ' stored{k}]);
                catch err
                    add('predictions', 'failed', sprintf('%s: %s', stored{k}, err.message));
                end
            end
        else
            add('predictions', status, note);
        end
    end
end

% 4. Saved ecModels: protein reactions turned forward
models = dir(fullfile(projectPath, 'models', '*.yml'));
for k = 1:numel(models)
    file = fullfile(models(k).folder, models(k).name);
    if ~contains(fileread(file), 'prot_pool')
        continue
    end
    try
        model = readYAMLmodel(file);
    catch err
        add('saved ecModel', 'failed', sprintf('%s: %s', relative(file, projectPath), err.message));
        continue
    end
    if isfield(model.annotation, 'gecko_format') && str2double(model.annotation.gecko_format) > 4
        add('saved ecModel', 'failed', sprintf('%s: ecModel format %s is newer than this GECKO', ...
            relative(file, projectPath), model.annotation.gecko_format));
        continue
    end
    [model, changed] = upgradeEcModel(model);
    if isempty(changed)
        continue
    end
    if apply
        model.annotation.generator = ['GECKO ' geckoVersion()];
        keepOriginal(file);
        writeYAMLmodel(model, file);
    end
    add('saved ecModel', status, sprintf('%s: %d protein reactions turned forward', ...
        relative(file, projectPath), numel(changed)));
end

% 5. project_format in model_adapter.toml
adapterRow = strcmp(rows(:, 1), 'adapter') & ~strcmp(rows(:, 2), 'by hand');
if isfile(toml) || any(adapterRow)
    if ~isfile(toml) || isempty(regexp(fileread(toml), '(?m)^\s*project_format\s*=', 'once'))
        if apply
            text = fileread(toml);
            fid = fopen(toml, 'w', 'n', 'UTF-8');
            fprintf(fid, 'project_format = 4\n%s', text);
            fclose(fid);
        end
        add('project format', status, 'project_format = 4 in model_adapter.toml');
    end
end

% 6. Scripts that call renamed or replaced functions
renamed = { ...
    'enzymeUsage', 'getEnzymeUsage'; 'pfbaEnzymes', 'getPfbaEnzymes'; ...
    'selectKcatValue', 'assignKcatValues'; 'startGECKOproject', 'createGECKOproject'; ...
    'updateGECKOdoc', 'buildGECKOdoc'; 'copyECtoGEM', 'applyECcodes'; ...
    'addCarbonNum', 'fillCarbonNum'; 'updateprior', 'updatePrior'; ...
    'mergeDLKcatAndFuzzyKcats', 'mergeKcats'; ...
    'fuzzyKcatMatching', 'kcatsFromBrenda (or kcatsFromECTree), then assignKcats'; ...
    'readDLKcatOutput', 'kcatsFromDLKcat, then assignKcats'; ...
    'readOpenKineticsPredictorOutput', 'kcatsFromOKP, then assignKcats'; ...
    'writeOpenKineticsPredictorInput', 'kcatsFromOKP'; ...
    'applyCustomKcats', 'data/kcats/custom.tsv and assignKcats'};
scripts = dir(fullfile(projectPath, '**', '*.m'));
for k = 1:numel(scripts)
    file = fullfile(scripts(k).folder, scripts(k).name);
    if startsWith(file, fullfile(projectPath, 'legacy'))
        continue
    end
    lines = strsplit(fileread(file), {'\r\n', '\n'});
    for r = 1:size(renamed, 1)
        hit = find(~cellfun(@isempty, regexp(lines, ['(?<![\w.])' renamed{r, 1} '\s*\('], 'once')));
        for h = hit(:)'
            add('script', 'by hand', sprintf('%s:%d: %s -> %s', relative(file, projectPath), h, ...
                renamed{r, 1}, renamed{r, 2}));
        end
    end
end

report = cell2table(rows, 'VariableNames', {'step', 'status', 'note'});
if isempty(rows)
    fprintf('%s has the current GECKO 4 layout.\n', projectPath);
else
    disp(report)
    if ~apply && any(strcmp(rows(:, 2), 'to do'))
        fprintf('Run upgradeGECKOproject(projectPath, ''dryRun'', false) to carry out the steps; originals go to %s.\n', legacy);
    end
end
end

function adapter = projectAdapter(projectPath)
% The project's adapter: from model_adapter.toml, or its adapter class.
if isfile(fullfile(projectPath, 'model_adapter.toml'))
    adapter = TomlAdapter(projectPath);
else
    adapter = ModelAdapterManager.getAdapter(adapterClass(projectPath));
    adapter.params.path = projectPath;
end
end

function file = adapterClass(projectPath)
% The adapter class in the project folder, or '' when there is none.
file = '';
candidates = dir(fullfile(projectPath, '*.m'));
for k = 1:numel(candidates)
    f = fullfile(candidates(k).folder, candidates(k).name);
    if ~isempty(regexp(fileread(f), '(?m)^\s*classdef\s+\w+\s*<\s*ModelAdapter', 'once'))
        file = f;
        return
    end
end
end

function method = okpMethod(projectPath, outName)
% The OKP method of a stored result: from its own job file
% (OKP_job_<tag>.txt), or from OKP_job.txt for the shared OKP_output.csv.
method = '';
if strcmp(outName, 'OKP_output.csv')
    jobFile = fullfile(projectPath, 'data', 'OKP_job.txt');
else
    jobFile = fullfile(projectPath, 'data', ['OKP_job_' regexprep(outName, '^OKP_output_|\.csv$', '') '.txt']);
end
if ~isfile(jobFile)
    return
end
tok = regexp(fileread(jobFile), '(?m)^method:\s*(\S+)', 'tokens', 'once');
if ~isempty(tok)
    method = tok{1};
end
end

function r = relative(file, projectPath)
r = file(numel(projectPath)+2:end);
end

function checkOptions(args, names, functionName)
% Options are name-value pairs with known names, so a misspelled option is
% an error instead of being taken as a positional argument.
id = [functionName ':arguments'];
given = args(1:2:end);
if mod(numel(args), 2) ~= 0 || ~all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), given))
    error(id, 'Options of %s are name-value pairs. Options: %s.', functionName, strjoin(names, ', '))
end
unknown = setdiff(cellfun(@char, given, 'UniformOutput', false), names);
if ~isempty(unknown)
    error(id, 'Unknown option of %s: %s. Options: %s.', functionName, strjoin(unknown, ', '), ...
        strjoin(names, ', '))
end
end

function v = geckoVersion()
file = fullfile(findGECKOroot, 'version.txt');
if isfile(file)
    v = strtrim(fileread(file));
else
    v = 'unknown';
end
end
