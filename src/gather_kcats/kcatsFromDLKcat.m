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
% - 'file': the output of an earlier run in data/DLKcat.tsv
%   (readDLKcatOutput), for example of a GECKO 3 project.
% - 'docker': a local run of DLKcat in its Docker image (writeDLKcatInput,
%   runDLKcat), which needs Docker. This is the GECKO 3 method; 'okp' runs
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
%     the kcat values, as returned by readDLKcatOutput or
%     readOpenKineticsPredictorOutput; empty while an OKP job runs.
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
% kcatsFromOKP, readDLKcatOutput, runDLKcat, saveKcatSource

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
    % As runDLKcat: on macOS, MATLAB started from the launcher may not
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
