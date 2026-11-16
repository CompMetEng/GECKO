function [kcatList, done] = kcatsFromDLKcat(model, varargin)
% kcatsFromDLKcat  Gather kcats from DLKcat and keep them as a source file.
%
% Predicts kcats with DLKcat and writes them to data/kcats/dlkcat.tsv in
% the project folder (see saveKcatSource). The option 'via' sets where the
% predictions come from:
%
% - 'okp': the DLKcat method of OpenKineticsPredictor (see kcatsFromOKP),
%   which needs an API key (see resolveOkpApiKey). OKP is a web service
%   with a job queue: the first call submits a job, a later call fetches
%   the result, and done is false while the job runs.
% - 'file': the output of an earlier run in data/DLKcat.tsv, for example
%   of a GECKO 3 project: a tab-separated file with a header and the
%   columns reaction, gene, substrate, SMILES, sequence and kcat. Rows
%   without a numeric kcat are left out, and rows with a substrate that is
%   not in model.metNames are left out with a warning.
% - 'docker': a local run of DLKcat in its Docker image, which needs
%   Docker. It writes the input to data/DLKcat.tsv (one row per
%   enzyme-substrate pair with a SMILES) and replaces it with the
%   predictions. This is the GECKO 3 method; 'okp' runs
%   the same predictor without Docker. A run records how many reactions it
%   covered in data/DLKcat_run.txt, so that the source file says so also
%   when the output is read later.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure), with metabolite SMILES
%     (findMetSmiles) and protein sequences.
%
% Name-Value Arguments
% --------------------
% via : char
%     'okp', 'file' or 'docker' (default: 'file' when data/DLKcat.tsv
%     exists, 'okp' otherwise).
% run : logical
%     true is the same as 'via', 'docker' (default false).
% wait : logical
%     with 'via', 'okp': wait until the job is done (default false).
% ecRxns : logical
%     for which model.ec.rxns DLKcat is run (default: all). Only used when
%     a prediction is started ('docker', or 'okp' when a job is submitted).
% save : logical
%     write the source file (default true).
% saveAs : char
%     name of the source file, without .tsv (default 'dlkcat').
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     the kcat values (see kcatsFromOKP for the fields of an OKP result);
%     empty while an OKP job runs.
% done : logical
%     whether the kcats were available.
%
% Examples
% --------
%     kcatsFromDLKcat(ecModel);                    % submits an OKP job
%     kcatsFromDLKcat(ecModel);                    % later: fetches and writes
%     kcatsFromDLKcat(ecModel, 'via', 'docker');   % local run with Docker
%
% See also
% --------
% kcatsFromOKP, extractEnzymeSubstratePairs, saveKcatSource

checkOptions(varargin, {'via','run','wait','ecRxns','save','saveAs','modelAdapter'}, 'kcatsFromDLKcat');
p = parseGECKOargs(varargin, { ...
    'via',          []; ...
    'run',          false; ...
    'wait',         false; ...
    'ecRxns',       []; ...
    'save',         true; ...
    'saveAs',       'dlkcat'; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
outFile = fullfile(modelAdapter.params.path, 'data', 'DLKcat.tsv');
if p.run
    if ~isempty(p.via) && ~strcmpi(p.via, 'docker')
        error('kcatsFromDLKcat:arguments', '''run'', true is the same as ''via'', ''docker''.')
    end
    via = 'docker';
elseif ~isempty(p.via)
    via = validatestring(char(p.via), {'okp', 'file', 'docker'});
elseif isfile(outFile)
    via = 'file';
else
    via = 'okp';
end
if strcmp(via, 'file') && ~isempty(p.ecRxns)
    error('kcatsFromDLKcat:ecRxns', ['ecRxns selects the reactions DLKcat is run for, so ' ...
        'it needs ''via'', ''okp'' or ''docker''.'])
end
if strcmp(via, 'okp')
    [kcatList, done] = kcatsFromOKP(model, 'method', 'DLKcat', 'wait', p.wait, ...
        'ecRxns', p.ecRxns, 'save', p.save, 'saveAs', p.saveAs, 'modelAdapter', modelAdapter);
    return
end
done = true;
if p.save
    file = fullfile(modelAdapter.params.path, 'data', 'kcats', [char(p.saveAs) '.tsv']);
    if strcmp(via, 'docker')
        warnSubset(p.ecRxns, numel(model.ec.rxns), file);
    end
end

runFile = fullfile(modelAdapter.params.path, 'data', 'DLKcat_run.txt');
if strcmp(via, 'docker')
    % As in runDLKcat below: on macOS, MATLAB started from the launcher may not
    % have Docker on its path.
    if ismac
        setenv('PATH', strcat('/usr/local/bin', ':', getenv("PATH")));
    end
    [status, ~] = system('docker --version');
    if status ~= 0
        error('kcatsFromDLKcat:noDocker', ['A local DLKcat run needs Docker, which was not ' ...
            'found. Install it from https://www.docker.com, or use ''via'', ''okp''.'])
    end
    writeDLKcatInput(model, 'ecRxns', p.ecRxns, 'modelAdapter', modelAdapter, 'overwrite', true);
    runDLKcat('modelAdapter', modelAdapter);
    fid = fopen(runFile, 'w');
    if ~isempty(p.ecRxns)
        fprintf(fid, 'ecRxns: %d\n', countRxns(p.ecRxns));
    end
    fclose(fid);
elseif ~isfile(outFile)
    error('kcatsFromDLKcat:noOutput', ['No DLKcat output at %s. Use ''via'', ''okp'' to ' ...
        'predict the kcats.'], outFile)
end
kcatList = readDLKcatOutput(model, 'modelAdapter', modelAdapter);

if p.save
    parameters = struct('run', strcmp(via, 'docker'));
    if isfile(runFile)
        tok = regexp(fileread(runFile), 'ecRxns:\s*(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            parameters.ec_rxns = str2double(tok{1});
        end
    end
    saveKcatSource(model, kcatList, p.saveAs, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromDLKcat', 'parameters', parameters);
    fprintf('DLKcat kcats for %d reactions written to %s\n', ...
        numel(unique(kcatList.rxns(kcatList.kcats > 0))), file);
end
end

function n = countRxns(ecRxns)
if islogical(ecRxns)
    n = sum(ecRxns);
else
    n = numel(unique(ecRxns));
end
end

function warnSubset(ecRxns, nRxns, file)
% A subset written under the default name replaces the full file.
if ~isempty(ecRxns) && countRxns(ecRxns) < nRxns && isfile(file)
    printOrange(sprintf(['%s is replaced by kcats for %d reactions only. Give saveAs ' ...
        'to keep it.\n'], file, countRxns(ecRxns)));
end
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

function writtenTable = writeDLKcatInput(model, varargin)
% Writes the DLKcat input (data/DLKcat.tsv by default): one row per
% enzyme-substrate pair, with the kcat column NA.
p = parseGECKOargs(varargin, { ...
    'ecRxns',         []; ...
    'modelAdapter',   []; ...
    'onlyWithSmiles', []; ...
    'filename',       []; ...
    'overwrite',      []});
ecRxnsSel      = p.ecRxns;
modelAdapter   = p.modelAdapter;
onlyWithSmiles = p.onlyWithSmiles;
filename       = p.filename;
overwrite      = p.overwrite;

if isempty(ecRxnsSel)
    ecRxnsSel = true(numel(model.ec.rxns),1);
elseif ~islogical(ecRxnsSel)
    error('ecRxns should be provided as logical vector')
elseif numel(ecRxnsSel)~=numel(model.ec.rxns)
    error('Length of ecRxns is not the same as model.ec.rxns')
end
ecRxnsSel = find(ecRxnsSel); % Change to indices, into the full model.ec.rxns

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if isempty(onlyWithSmiles)
    onlyWithSmiles=true;
end

if isempty(filename)
    filename = fullfile(params.path,'data','DLKcat.tsv');
elseif ~endsWith(filename,'.tsv')
    error('If filename is provided, it should include the .tsv extension.')
end

if isempty(overwrite) || ~overwrite % If is true
    if exist(filename,'file')
        error([filename ' already exists, either delete it first, or set the overwrite input argument as true'])
    end
end

pairs = extractEnzymeSubstratePairs(model, 'ecRxns', ismember((1:numel(model.ec.rxns))', ecRxnsSel), ...
    'modelAdapter', modelAdapter);

% Prepare output
if isempty(pairs.ecRxn)
    out = cell(6,0);
    writeRows(filename, out)
    writtenTable = out;
    return
end
out(1,:) = model.ec.rxns(pairs.ecRxn);
out(2,:) = model.ec.genes(pairs.protein);
out(3,:) = model.metNames(pairs.met);
if isfield(model,'metSmiles')
    out(4,:) = model.metSmiles(pairs.met);
else
    out(4,:) = cell(numel(pairs.met),1);
end

out(5,:) = model.ec.sequence(pairs.protein);
if onlyWithSmiles
    out(:,cellfun(@isempty,out(4,:))) = [];
else
    out(4,cellfun(@isempty,out(4,:))) = {'None'};
end
out(6,:) = cell(numel(out(1,:)),1);
out(6,:) = {'NA'};

writeRows(filename, out)
writtenTable = out;
end

function writeRows(filename, out)
% Write the rows (one per column of out); an empty table gives an empty file.
fID = fopen(filename,'w');
if ~isempty(out)
    fprintf(fID,'%s\t%s\t%s\t%s\t%s\t%s\n',out{:});
end
fclose(fID);
fprintf('Model-specific DLKcat input stored at %s\n',filename);
end

function runDLKcat(varargin)
% Runs DLKcat in its Docker image on data/DLKcat.tsv, which it overwrites
% with the predictions.
p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'filePath',     []});
modelAdapter = p.modelAdapter;
filePath     = p.filePath;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
% Make sure path is full, not relative
[~, params.path] = fileattrib(params.path);
params.path=params.path.Name;

if isempty(filePath)
    filePath = fullfile(params.path,'data','DLKcat.tsv');
elseif strcmp(filePath(end),{'\','/'})
    filePath = fullfile(filePath,'DLKcat.tsv');
end
filePath = checkFileExistence(filePath,1);

copyfile(filePath, fullfile(params.path,'data','tempDLKcat.tsv'));
    

%% Check and install requirements
% On macOS, Docker might not be properly loaded if MATLAB is started via
% launcher and not terminal.
if ismac
    setenv('PATH', strcat('/usr/local/bin', ':', getenv("PATH")));
end

% Check if Docker is installed
[checks.docker.status, checks.docker.out] = system('docker --version');
if checks.docker.status ~= 0
    error('Cannot find Docker, make sure it is installed. If it is, it might be required to start Matlab from the command-line instead of the launcher in order for Docker to be detected and used.')
end

disp('Running DLKcat prediction, this may take many minutes, especially the first time.')
status = system(['docker run --rm -v "' fullfile(params.path,'/data') '":/data ghcr.io/sysbiochalmers/dlkcat-gecko:0.1 /bin/bash -c "python DLKcat.py /data/tempDLKcat.tsv /data/tempDLKcatOutput.tsv"']);
delete(fullfile(params.path,'/data/tempDLKcat.tsv'));

if status == 0 && exist(fullfile(params.path,'data/tempDLKcatOutput.tsv'))
    movefile(fullfile(params.path,'/data/tempDLKcatOutput.tsv'), filePath);
    disp('DKLcat prediction completed.');
else    
    error('DLKcat encountered an error or it did not create any output file.')
end
end

function kcatList = readDLKcatOutput(model, varargin)
% Reads the DLKcat output (data/DLKcat.tsv by default) into a kcatList.
% Rows without a finite numeric kcat are left out; rows with a substrate
% that is not in model.metNames are left out with a warning
% (readDLKcatOutput:unknownSubstrate). A file without kcats raises
% readDLKcatOutput:noKcats, a reaction that is not in model.ec.rxns
% readDLKcatOutput:unknownReaction.
p = parseGECKOargs(varargin, { ...
    'outFile',      []; ...
    'modelAdapter', []});
outFile      = p.outFile;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if isempty(outFile)
    fID      = fopen(fullfile(params.path,'data','DLKcat.tsv'),'r');
else
    fID      = fopen(outFile);
end
DLKcatOutput = textscan(fID,'%s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fID);

[rxns, genes, subs, kcats] = deal(DLKcatOutput{[1,2,3,6]});

% A row has a kcat when its kcat field is a plain finite number.
isNumber = ~cellfun(@isempty, regexp(strtrim(kcats), '^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$', 'once'));
kcatValues = nan(numel(kcats), 1);
kcatValues(isNumber) = str2double(kcats(isNumber));
usable = isfinite(kcatValues);

% Check if it contains any kcat values
if ~any(usable)
    error('readDLKcatOutput:noKcats', 'DLKcat file does not contain any kcat values, please run DLKcat first.')
end

% Check that the substrates are in the model (case-insensitive: an SBML loader can
% normalise metabolite-name capitalisation differently than whatever produced the
% DLKcat input, and that must not read as a genuinely missing substrate). Rows with
% another substrate are left out with a warning.
matchMets = ismember(lower(subs),lower(model.metNames));
if ~all(matchMets)
    warning('readDLKcatOutput:unknownSubstrate', 'Rows with a substrate that is not in model.metNames are left out: %s', ...
        strjoin(unique(subs(~matchMets))', ', '))
    usable = usable & matchMets;
end

% Check that the reactions of the rows with a kcat are in model.ec.rxns
matchRxns = ismember(rxns,model.ec.rxns);
if ~all(matchRxns(usable))
    errorText = 'DLKcat was likely run with an input file that was generated from another ecModel, as the following reactions from DLKcat output cannot be found in model.ec.rxns:';
    error('readDLKcatOutput:unknownReaction', '%s', ravenList(errorText,rxns(usable & ~matchRxns),false))
end

% Keep the entries with a numeric value
kcats = kcatValues(usable);
rxns  = rxns(usable);
genes = genes(usable);
subs  = subs(usable);

% Make kcatList structure
kcatList.source     = 'DLKcat';
kcatList.rxns       = rxns;
kcatList.genes      = genes;
kcatList.substrates = subs;
kcatList.kcats      = kcats;
end
