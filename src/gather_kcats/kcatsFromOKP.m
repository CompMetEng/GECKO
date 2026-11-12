function [kcatList, done] = kcatsFromOKP(model, varargin)
% kcatsFromOKP  Gather kcats from OpenKineticsPredictor and keep them as a source file.
%
% Predicts kcats with one OpenKineticsPredictor (OKP) method and writes
% them to data/kcats/okp_<method>.tsv in the project folder (see
% saveKcatSource). OKP is a web service with a job queue, so this function
% can be called repeatedly:
%
% - without a job for this method, it writes the input of the job
%   (data/OKP.csv: protein sequences and single-substrate SMILES) and
%   submits it;
% - with a job, it fetches the result. When the job is done the file is
%   written; otherwise it reports the job's status and writes nothing.
%
% Each method keeps its own job in data/OKP_job_<method>.txt and its own
% result in data/OKP_output_<method>.csv, so jobs of several methods can
% run at the same time. The job file records which job the stored result
% came from; useStored only uses a result of the current job.
%
% An API key is needed. Keys are free: open
% https://predictor.openkinetics.org/api-docs, click "Generate" in the "API
% key generator" section and copy the key (it is shown once). A key has a
% daily prediction quota. The key is a secret and does not go in the model
% adapter; resolveOkpApiKey says where it is read from.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure), with
%     metabolite SMILES (findMetSmiles) and protein sequences.
%
% Name-Value Arguments
% --------------------
% method : char
%     OKP method, e.g. 'CataPro', 'CatPred', 'CatRange', 'DLKcat',
%     'EITLEM', 'KinForm-H', 'KinForm-L', 'OmniESI', 'TurNup', 'UniKP'
%     (default: params.okp.method of the model adapter, or 'KinForm-H').
% wait : logical
%     wait until the job is done (default false).
% resubmit : logical
%     submit a new job even if one exists for this method; its earlier
%     result is removed (default false).
% useStored : logical
%     read the stored result of this method's job, without contacting OKP
%     (default false).
% ecRxns : logical
%     for which model.ec.rxns kcats are predicted (default: all). Only
%     used when a job is submitted.
% save : logical
%     write the source file (default true).
% saveAs : char
%     name of the source file, without .tsv (default 'okp_<method>', in
%     lower case, with other characters than letters and digits as '_').
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     the kcat values, with the fields source ('OpenKineticsPredictor'),
%     rxns, genes, substrates, kcats, similarity (the enzyme's highest
%     sequence similarity to the method's training data, in percent; NaN
%     without it) and kcatSource (the OKP source of each value, such as
%     'OKP-KinForm-H' or 'OKP-BRENDA'); empty while the job is not done.
%     A CatRange class ('10 to 100 s^-1') is read as its geometric
%     midpoint.
% done : logical
%     whether the result was available.
%
% Examples
% --------
%     kcatsFromOKP(ecModel, 'method', 'CataPro');          % submits
%     kcatsFromOKP(ecModel, 'method', 'CataPro');          % later: fetches and writes
%     kcatsFromOKP(ecModel, 'method', 'UniKP', 'wait', true);
%
% See also
% --------
% kcatsFromBrenda, kcatsFromDLKcat, resolveOkpApiKey, saveKcatSource

checkOptions(varargin, {'method','wait','resubmit','useStored','ecRxns','save','saveAs','modelAdapter'}, 'kcatsFromOKP');
p = parseGECKOargs(varargin, { ...
    'method',       []; ...
    'wait',         false; ...
    'resubmit',     false; ...
    'useStored',    false; ...
    'ecRxns',       []; ...
    'save',         true; ...
    'saveAs',       []; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
method = char(p.method);
if isempty(method)
    if isfield(params, 'okp') && isfield(params.okp, 'method') && ~isempty(params.okp.method)
        method = char(params.okp.method);
    else
        method = 'KinForm-H';
    end
end
if any(strcmp(method, {'CatPred', 'TurNup'}))
    error('kcatsFromOKP:method', ['%s reads whole reactions (its substrate set, and products ' ...
        'for TurNup); kcatsFromOKP submits enzyme-substrate pairs only.'], method)
end
tag = lower(regexprep(method, '[^A-Za-z0-9]+', '_'));
dataDir = fullfile(params.path, 'data');
jobFile = fullfile(dataDir, ['OKP_job_' tag '.txt']);
outFile = fullfile(dataDir, ['OKP_output_' tag '.csv']);
name = char(p.saveAs);
if isempty(name)
    name = ['okp_' tag];
end
job = readJobFile(jobFile);
submitting = ~p.useStored && (p.resubmit || ~isfield(job, 'jobId'));
if ~isempty(p.ecRxns) && ~submitting
    error('kcatsFromOKP:ecRxns', ['ecRxns selects the reactions of a new job, but this call ' ...
        'submits none (a job for %s exists; use ''resubmit'', true for a new one).'], method)
end
if p.save
    file = fullfile(modelAdapter.params.path, 'data', 'kcats', [char(name) '.tsv']);
    if submitting
        warnSubset(p.ecRxns, numel(model.ec.rxns), file);
    end
end

if p.useStored
    if ~isfield(job, 'jobId')
        adoptSharedResult(model, modelAdapter, dataDir, method, jobFile, outFile);
        job = readJobFile(jobFile);
    end
    if ~isfile(outFile) || ~isfield(job, 'resultOf') || ~strcmp(job.resultOf, job.jobId)
        error('kcatsFromOKP:noStored', ['No stored result of the current %s job (%s). Call ' ...
            'kcatsFromOKP without ''useStored'' to fetch it.'], method, jobFile)
    end
    kcatList = readOpenKineticsPredictorOutput(model, 'outFile', outFile, 'modelAdapter', modelAdapter);
    done = true;
else
    if submitting
        if isfile(outFile)
            delete(outFile);
        end
        jobId = submitOpenKineticsPredictor(model, method, p.ecRxns, modelAdapter);
        job = struct('jobId', jobId, 'method', method, 'submittedAt', ...
            char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
        if ~isempty(p.ecRxns)
            job.ecRxns = sprintf('%d', countRxns(p.ecRxns));
        end
        writeJobFile(jobFile, job);
    end
    [done, kcatList] = fetchOpenKineticsPredictor(model, job.jobId, p.wait, outFile, modelAdapter);
    if done
        job.resultOf = job.jobId;
        writeJobFile(jobFile, job);
    end
end

if ~done
    fprintf('Call kcatsFromOKP(model, ''method'', ''%s'') again later to collect the result.\n', method);
    return
end
if p.save
    parameters = struct('method', method);
    if isfield(job, 'ecRxns')
        parameters.ec_rxns = str2double(job.ecRxns);
    end
    if isfield(kcatList, 'similarity')
        kcatList.detail = repmat({''}, numel(kcatList.rxns), 1);
        has = ~isnan(kcatList.similarity);
        kcatList.detail(has) = cellstr("similarity=" + string(kcatList.similarity(has)));
    end
    saveKcatSource(model, kcatList, name, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromOKP', 'parameters', parameters, 'data', ['OKP job ' job.jobId]);
    fprintf('OKP (%s) kcats for %d reactions written to %s\n', method, ...
        numel(unique(kcatList.rxns(kcatList.kcats > 0))), file);
end
end

function adoptSharedResult(model, modelAdapter, dataDir, method, jobFile, outFile)
% A result that a GECKO 4 beta stored outside kcatsFromOKP, in
% data/OKP_output.csv with data/OKP_job.txt, is taken over for this method
% only when its job is for this method and its predictions are labelled
% with this method.
shared = fullfile(dataDir, 'OKP_output.csv');
sharedJob = readJobFile(fullfile(dataDir, 'OKP_job.txt'));
ok = isfile(shared) && isfield(sharedJob, 'jobId') && isfield(sharedJob, 'method') && ...
    strcmpi(sharedJob.method, method);
if ok
    list = readOpenKineticsPredictorOutput(model, 'outFile', shared, 'modelAdapter', modelAdapter);
    predicted = list.kcatSource(startsWith(list.kcatSource, 'OKP-') & ...
        ~ismember(lower(list.kcatSource), {'okp-brenda','okp-sabio-rk','okp-uniprot'}));
    ok = ~isempty(predicted) && all(strcmpi(predicted, ['OKP-' method]));
end
if ~ok
    error('kcatsFromOKP:noStored', ['No stored OKP result for %s (%s). Call kcatsFromOKP ' ...
        'without ''useStored'' to submit or fetch a job.'], method, outFile)
end
copyfile(shared, outFile);
sharedJob.method = method;
sharedJob.resultOf = sharedJob.jobId;
writeJobFile(jobFile, sharedJob);
end

function job = readJobFile(file)
% Fields of an OKP job file ('key: value' per line), as a struct.
job = struct();
if ~isfile(file)
    return
end
lines = strsplit(fileread(file), newline);
for i = 1:numel(lines)
    tok = regexp(lines{i}, '^([A-Za-z]\w*):\s*(.*?)\s*$', 'tokens', 'once');
    if ~isempty(tok)
        job.(tok{1}) = tok{2};
    end
end
end

function writeJobFile(file, job)
fid = fopen(file, 'w');
for f = {'jobId','method','submittedAt','ecRxns','resultOf'}
    if isfield(job, f{1})
        fprintf(fid, '%s: %s\n', f{1}, job.(f{1}));
    end
end
fclose(fid);
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

function jobId = submitOpenKineticsPredictor(model, method, ecRxns, modelAdapter)
% Writes the input of a job (data/OKP.csv) and submits it to the OKP REST API
% (https://predictor.openkinetics.org/). The input is written before the API
% key is resolved, so it can be checked without a key.
params = modelAdapter.params;
if isempty(ecRxns)
    ecRxns = true(numel(model.ec.rxns),1);
elseif ~islogical(ecRxns) || numel(ecRxns) ~= numel(model.ec.rxns)
    error('ecRxns should be a logical vector with length equal to model.ec.rxns')
end
% Request options: defaults, then params.okp, then the method of this job.
okp = struct('method','KinForm-H', ...
             'handleLongSequences','truncate', ...
             'includeSimilarityColumns',true,'canonicalizeSubstrates',true);
if isfield(params,'okp')
    fn = fieldnames(params.okp);
    for i=1:numel(fn); okp.(fn{i}) = params.okp.(fn{i}); end
end
okp.method = method;

filename = fullfile(params.path,'data','OKP.csv');
writeOpenKineticsPredictorInput(model, 'ecRxns', ecRxns, 'modelAdapter', modelAdapter, ...
    'onlyWithSmiles', true, 'overwrite', true);
apiKey = resolveOkpApiKey('', params.path);

%% Submit via the OKP REST API
% GECKO only ever predicts kcat.
url = 'https://predictor.openkinetics.org/api/v1/submit/';

targetsJson = '["kcat"]';
methodsJson = ['{"kcat":"' okp.method '"}'];

import matlab.net.http.*
import matlab.net.http.io.*
provider = MultipartFormProvider( ...
    'file',                     FileProvider(filename), ...
    'targets',                  targetsJson, ...
    'methods',                  methodsJson, ...
    'handleLongSequences',      okp.handleLongSequences, ...
    'includeSimilarityColumns', boolToStr(okp.includeSimilarityColumns), ...
    'canonicalizeSubstrates',   boolToStr(okp.canonicalizeSubstrates));
header = matlab.net.http.HeaderField('Authorization', ['Bearer ' apiKey]);
req    = RequestMessage('POST', header, provider);
resp   = req.send(url);

data = resp.Body.Data;
if resp.StatusCode ~= matlab.net.http.StatusCode.OK && ...
        resp.StatusCode ~= matlab.net.http.StatusCode.Created
    error('OKP submit failed (HTTP %d): %s', double(resp.StatusCode), okpErrorText(data));
end
if ~isstruct(data) || ~isfield(data,'jobId')
    error('OKP submit returned an unexpected response (no jobId).');
end
jobId = data.jobId;
fprintf('Submitted OKP job %s (method: %s).\n', jobId, okp.method);
end

function [done, kcatList] = fetchOpenKineticsPredictor(model, jobId, wait, outFile, modelAdapter)
% Checks the status of an OKP job and, when it has finished, stores the
% result in outFile and reads it. With wait, polls every 30 s until the job
% has finished.
params = modelAdapter.params;
pollInterval = 30;
apiKey = resolveOkpApiKey('', params.path);
base = 'https://predictor.openkinetics.org/api/v1';
opts = weboptions('HeaderFields', {'Authorization', ['Bearer ' apiKey]}, ...
                   'Timeout', 60, 'ContentType', 'json');

while true
    statusResp = webread([base '/status/' jobId '/'], opts);
    state = lower(string(statusResp.status));

    switch state
        case "completed"
            csvText = webread([base '/result/' jobId '/'], ...
                weboptions('HeaderFields', {'Authorization', ['Bearer ' apiKey]}, ...
                           'Timeout', 120, 'ContentType', 'text'));
            fID = fopen(outFile,'w');
            fwrite(fID, csvText);
            fclose(fID);
            fprintf('OKP job %s completed; result stored at %s\n', jobId, outFile);
            kcatList = readOpenKineticsPredictorOutput(model, 'outFile', outFile, 'modelAdapter', modelAdapter);
            done = true;
            return

        case "failed"
            error('OKP job %s failed. Check https://predictor.openkinetics.org/ for details.', jobId);

        otherwise % pending / running / queued
            pct = '';
            if isfield(statusResp,'progress') && isfield(statusResp.progress,'predictionsMade') ...
                    && isfield(statusResp.progress,'predictionsTotal') && statusResp.progress.predictionsTotal > 0
                pct = sprintf(' (%d/%d predictions)', statusResp.progress.predictionsMade, statusResp.progress.predictionsTotal);
            end
            if wait
                fprintf('OKP job %s status: %s%s. Waiting %d s...\n', jobId, char(statusResp.status), pct, pollInterval);
                pause(pollInterval);
            else
                fprintf('OKP job %s not finished (status: %s%s).\n', ...
                         jobId, char(statusResp.status), pct);
                done = false;
                kcatList = [];
                return
            end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = boolToStr(b)
if islogical(b) || isnumeric(b)
    if b; s = 'true'; else; s = 'false'; end
else
    s = char(b);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = okpErrorText(data)
if isstruct(data) && isfield(data,'error')
    txt = data.error;
elseif ischar(data) || isstring(data)
    txt = char(data);
else
    txt = 'unknown error';
end
end

function writeOpenKineticsPredictorInput(model, varargin)
% Writes data/OKP.csv, the input of a job: protein sequences and
% single-substrate SMILES.
[geckoPath, ~] = findGECKOroot();

%% Parse inputs
p = parseGECKOargs(varargin, { ...
    'ecRxns',         []; ...
    'modelAdapter',   []; ...
    'onlyWithSmiles', []; ...
    'overwrite',      []});
ecRxns         = p.ecRxns;
modelAdapter   = p.modelAdapter;
onlyWithSmiles = p.onlyWithSmiles;
overwrite      = p.overwrite;

if isempty(ecRxns)
    ecRxns = true(numel(model.ec.rxns),1);
elseif ~islogical(ecRxns) || numel(ecRxns) ~= numel(model.ec.rxns)
    error('ecRxns should be a logical vector with length equal to model.ec.rxns')
end

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either provide a modelAdapter or set the default in ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if isempty(onlyWithSmiles)
    onlyWithSmiles = true;
end

if isempty(overwrite)
    overwrite = false;
end

%% Output files
filename = fullfile(params.path,'data','OKP.csv');
if ~overwrite && exist(filename,'file')
    error('%s already exists. Set overwrite=true to replace it.', filename)
end

%% Map ec.rxns to model.rxns
if model.ec.geckoLight
   origRxns = extractAfter(model.ec.rxns, 4);
else
   origRxns = model.ec.rxns;
end

ecRxnsIdx = find(ecRxns);
[found, origRxnIdxs] = ismember(origRxns(ecRxnsIdx), model.rxns);
if ~all(found)
    error('Not all reactions in model.ec.rxns found in model.rxns')
end

%% Build reduced stoichiometric matrix (ignore certain metabolites)
% Load metabolites to ignore (reuses DLKcat files)
if exist(fullfile(params.path,'data','DLKcatIgnoreMets.tsv'),'file')
    ignoreFile = fullfile(params.path,'data','DLKcatIgnoreMets.tsv');
else
    ignoreFile = fullfile(geckoPath,'databases','DLKcatIgnoreMets.tsv');
end
fID = fopen(ignoreFile);
fileData = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[ignoreMets, ignoreSmiles] = deal(fileData{[1,2]});

% Normalize metabolite names (remove special characters, lowercase)
metsNorm        = lower(regexprep(model.metNames,'[^0-9a-zA-Z]+',''));
ignoreMetsNorm  = lower(regexprep(ignoreMets,'[^0-9a-zA-Z]+',''));
ignoreSmiles(cellfun(@isempty,ignoreSmiles)) = [];

% Mark metabolites to ignore
ignoreMetsIdx = ismember(metsNorm, ignoreMetsNorm);
if isfield(model,'metSmiles')
    ignoreMetsIdx = ignoreMetsIdx | ismember(model.metSmiles, ignoreSmiles);
end
ignoreMetsIdx = ignoreMetsIdx | startsWith(model.mets,'prot_');

reducedS = model.S;
reducedS(ignoreMetsIdx,:) = 0;

%% Remove currency metabolite pairs (reuses DLKcat files)
if exist(fullfile(params.path,'data','DLKcatCurrencyMets.tsv'),'file')
    currencyFile = fullfile(params.path,'data','DLKcatCurrencyMets.tsv');
else
    currencyFile = fullfile(geckoPath,'databases','DLKcatCurrencyMets.tsv');
end
fID = fopen(currencyFile);
fileData = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[currencyMets(:,1), currencyMets(:,2)] = deal(fileData{[1,2]});
currencyMets = lower(regexprep(currencyMets,'[^0-9a-zA-Z]+',''));
for i=1:size(currencyMets,1)
    subs = strcmp(currencyMets(i,1),metsNorm);
    prod = strcmp(currencyMets(i,2),metsNorm);
    [~,subsRxns]=find(reducedS(subs,:));
    [~,prodRxns]=find(reducedS(prod,:));
    pairRxns = intersect(subsRxns,prodRxns);
    tempRedS=reducedS;
    tempRedS([find(subs);find(prod)],pairRxns) = 0;
    % Do not remove currency mets if no substrate remains
    rxnsWithRemainingSubstrates = any(tempRedS(:,pairRxns) < 0,1);
    reducedS([find(subs);find(prod)],intersect(pairRxns,pairRxns(rxnsWithRemainingSubstrates))) = 0;
end


%% Extract substrates for selected reactions only
clearedS = reducedS(:, origRxnIdxs);

% Find all substrates (negative stoichiometry)
[substrateIdxs, reactionIdxs] = find(clearedS < 0);

% Find proteins catalyzing each reaction. reactionIdxs indexes clearedS's
% columns, i.e. positions within the requested ecRxnsIdx subset -- map back
% through ecRxnsIdx to the absolute model.ec.rxnEnzMat row before indexing it.
[proteinIdxs, rxnEnzIdxs] = find(model.ec.rxnEnzMat(ecRxnsIdx(reactionIdxs),:)');

%% Build output: sequence and SMILES
sequences = model.ec.sequence(proteinIdxs);

if isfield(model,'metSmiles')
    smiles = model.metSmiles(substrateIdxs(rxnEnzIdxs));
else
    smiles = cell(size(sequences));
end

%% Filter empty entries
if onlyWithSmiles
    valid = ~cellfun(@isempty, sequences) & ~cellfun(@isempty, smiles);
else
    valid = ~cellfun(@isempty, sequences);
    smiles(cellfun(@isempty, smiles)) = {'None'};
end

sequences = sequences(valid);
smiles = smiles(valid);

%% Write output file
% Sequence and SMILES
outTable = [sequences(:)'; smiles(:)'];

fID = fopen(filename,'w');
fprintf(fID,'Protein Sequence,Substrate\n');
fprintf(fID,'%s,%s\n', outTable{:});
fclose(fID);
fprintf('OpenKineticsPredictor input stored at %s\n', filename);
end

function kcatList = readOpenKineticsPredictorOutput(model, varargin)
% Reads an OKP result file (kcat (1/s), Source kcat, Extra Info kcat,
% Protein Sequence, Substrate) into a kcatList. A CatRange class ('10 to
% 100 s^-1') is read as its geometric midpoint, or half the upper bound for
% the class from 0. kcatSource keeps the OKP source of each value, with
% the prefix 'OKP-'.
p = parseGECKOargs(varargin, { ...
    'outFile',      []; ...
    'modelAdapter', []});
outFile      = p.outFile;

%% Read OKP output file using readcell to correctly read quoted fields
% that may contain commas in the 'Extra Info kcat' column.
if ~exist(outFile,'file')
    error('Cannot find file: %s', outFile)
end
raw = readcell(outFile, 'Delimiter', ',', 'TextType', 'char');

if size(raw,2) < 5
    error('OKP output file does not have the expected 5 columns (kcat, Source kcat, Extra Info kcat, Protein Sequence, Substrate): %s', outFile)
end

% Drop header row; the 'max similarity to <method> training data' column,
% when present, gives the similarity.
header = cellfun(@(h) char(string(h)), raw(1,:), 'UniformOutput', false);
simCol = find(startsWith(header, 'max similarity to'), 1);
raw = raw(2:end, :);
similarity = nan(size(raw, 1), 1);
if ~isempty(simCol)
    for i = 1:size(raw, 1)
        v = raw{i, simCol};
        if isnumeric(v) && ~isempty(v)
            similarity(i) = v;
        elseif ischar(v) || isstring(v)
            similarity(i) = str2double(v);
        end
    end
end

kcatsCol      = raw(:,1);
sourceKcatCol = raw(:,2);
sequences     = raw(:,4);
smiles        = raw(:,5);

% Normalize types: readcell returns doubles for numeric cells, char for
% strings, and `missing` for empty fields.
kcats = nan(size(kcatsCol));
for i = 1:numel(kcatsCol)
    v = kcatsCol{i};
    if isnumeric(v)
        kcats(i) = v;
    elseif ischar(v) || isstring(v)
        kcats(i) = kcatFromText(char(v));
    end
end

sourceKcatCol = normalizeCellStrings(sourceKcatCol);
sequences     = normalizeCellStrings(sequences);
smiles        = normalizeCellStrings(smiles);

%% Filter out entries with no numeric kcat value (NaN). Zero values are
% retained here and filtered downstream by assignKcatValues.
valid = ~isnan(kcats);
kcats         = kcats(valid);
sequences     = sequences(valid);
smiles        = smiles(valid);
sourceKcatCol = sourceKcatCol(valid);
similarity    = similarity(valid);

if isempty(kcats)
    error('OKP file does not contain any numeric kcat values. Run the job at https://predictor.openkinetics.org/ first.')
end

if ~isfield(model,'metSmiles')
    error('model.metSmiles is required to map SMILES back to metabolites.')
end

%% Build ec.rxns -> model.rxns index mapping (inverse of the mapping done
%  in writeOpenKineticsPredictorInput)
if model.ec.geckoLight
    origRxns = extractAfter(model.ec.rxns, 4);
else
    origRxns = model.ec.rxns;
end
[foundEcRxn, ecRxnToModelRxn] = ismember(origRxns, model.rxns);
if ~all(foundEcRxn)
    error('Not all entries in model.ec.rxns could be matched to model.rxns.')
end

%% Resolve each row to matching ec.rxns indices
rxnsOut       = cell(0,1);
genesOut      = cell(0,1);
substratesOut = cell(0,1);
kcatsOut      = zeros(0,1);
sourceOut     = cell(0,1);
similarityOut = zeros(0,1);
unmatched     = 0;

for i = 1:numel(sequences)
    % Map protein sequence -> ec.genes indices (columns of rxnEnzMat)
    proteinIdxs = find(strcmp(model.ec.sequence, sequences{i}));
    % Map SMILES -> metabolite indices (rows of model.S)
    metIdxs = find(strcmp(model.metSmiles, smiles{i}));

    if isempty(proteinIdxs) || isempty(metIdxs)
        unmatched = unmatched + 1;
        continue
    end

    % model.rxns columns where any of these mets are substrates (S < 0)
    modelRxnsUsingSub = any(model.S(metIdxs, :) < 0, 1);

    % Lift to ec.rxns: does the underlying model.rxn consume this metabolite?
    ecRxnUsesSub = modelRxnsUsingSub(ecRxnToModelRxn).';

    % ec.rxns catalyzed by any of these protein(s)
    ecRxnCatalyzed = any(model.ec.rxnEnzMat(:, proteinIdxs) ~= 0, 2);

    candidateEcRxns = find(ecRxnUsesSub & ecRxnCatalyzed);

    if isempty(candidateEcRxns)
        unmatched = unmatched + 1;
        continue
    end

    for k = 1:numel(candidateEcRxns)
        ecIdx = candidateEcRxns(k);

        % Pick the gene whose protein both matches this sequence and
        % catalyzes this ec.rxn
        catalysts       = find(model.ec.rxnEnzMat(ecIdx, :) ~= 0);
        matchingProtein = intersect(proteinIdxs, catalysts, 'stable');
        if ~isempty(matchingProtein)
            geneName = model.ec.genes{matchingProtein(1)};
        else
            geneName = '';
        end

        % Pick a metabolite name whose SMILES matches and that is a
        % substrate of the underlying model.rxn
        origRxnIdx       = ecRxnToModelRxn(ecIdx);
        metIsSub         = model.S(metIdxs, origRxnIdx) < 0;
        substrateMetIdxs = metIdxs(metIsSub);
        if ~isempty(substrateMetIdxs)
            subName = model.metNames{substrateMetIdxs(1)};
        else
            subName = '';
        end

        rxnsOut{end+1,1}       = model.ec.rxns{ecIdx};
        genesOut{end+1,1}      = geneName;
        substratesOut{end+1,1} = subName;
        kcatsOut(end+1,1)      = kcats(i);
        sourceOut{end+1,1}     = sourceKcatCol{i};
        similarityOut(end+1,1) = similarity(i);
    end
end

if unmatched > 0
    warning('%d of %d OKP entries could not be matched to an ec.rxn (protein sequence and/or SMILES not found in model, or no enzyme+substrate pair matches).', ...
        unmatched, numel(sequences));
end

if isempty(rxnsOut)
    error('No OKP entries could be mapped to model.ec.rxns.')
end

%% Build per-entry source labels
% Strip leading 'Prediction from ' to keep labels compact, then prefix all
% with 'OKP-' so the overall pipeline stage is still identifiable in
% model.ec.source.
kcatSource = regexprep(sourceOut, '^Prediction from\s+', '');
kcatSource = strcat('OKP-', kcatSource);

%% Make kcatList structure
kcatList.source     = 'OpenKineticsPredictor';
kcatList.rxns       = rxnsOut;
kcatList.genes      = genesOut;
kcatList.substrates = substratesOut;
kcatList.kcats      = kcatsOut;
kcatList.kcatSource = kcatSource;
kcatList.similarity = similarityOut;
end

function c = normalizeCellStrings(c)
% Convert each cell entry to a char row vector; missing/empty -> ''.
for i = 1:numel(c)
    v = c{i};
    if ismissing(v)
        c{i} = '';
    elseif isstring(v)
        c{i} = char(v);
    elseif isnumeric(v)
        if isnan(v)
            c{i} = '';
        else
            c{i} = num2str(v);
        end
    elseif ~ischar(v)
        c{i} = char(string(v));
    end
end
end

function kcat = kcatFromText(v)
% A number, or a CatRange class range '<low> to <high> s^-1' as its
% geometric midpoint (half the upper bound when the class starts at 0).
kcat = str2double(v);
if ~isnan(kcat)
    return
end
bounds = regexp(strtrim(v), '^([0-9.eE+-]+)\s+to\s+([0-9.eE+-]+)(\s*s\^-1)?$', 'tokens', 'once');
if isempty(bounds)
    return
end
low = str2double(bounds{1});
high = str2double(bounds{2});
if isnan(low) || isnan(high)
    return
elseif low > 0
    kcat = sqrt(low * high);
else
    kcat = high / 2;
end
end
