function [model,noSMILES] = findMetSmiles(model, varargin)
% findMetSmiles  Query PubChem by metabolite names to obtain SMILES.
%
% Matches are cached in <params.path>/data/smilesDB.tsv, which is checked
% first on subsequent runs so repeated queries are avoided. A name for which
% PubChem answers that it has no compound with that name (HTTP status 404, or
% 400 for a name it cannot parse) is cached with an empty SMILES. A request
% that gets no answer (server error, timeout, no connection) is retried up to
% 10 times; if it still fails, a warning (findMetSmiles:requestFailed) is
% issued, nothing is cached for that name and the name is returned in noSMILES,
% so that a later call asks again. Metabolites named 'prot_...' are skipped
% and never cached. If the model
% already has a metSmiles field, non-empty entries are left unchanged. When
% the cache holds a name more than once, the last row is used. A metabolite
% without a name is skipped with a warning (findMetSmiles:emptyName) and
% gets no SMILES.
%
% Parameters
% ----------
% model : struct
%     a model whose metNames field is used to find the relevant SMILES.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% verbose : logical
%     whether progress should be reported (default true).
%
% Returns
% -------
% model : struct
%     model with model.metSmiles specified.
% noSMILES : cell
%     metabolite names for which no SMILES could be found.
%
p = parseRAVENargs(varargin, { ...
    'modelAdapter', []; ...
    'verbose',      true});
modelAdapter = p.modelAdapter;
verbose      = p.verbose;
if isempty(verbose); verbose = true; end

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

nEmpty = sum(cellfun(@isempty,model.metNames));
if nEmpty > 0
    warning('findMetSmiles:emptyName','%d metabolite(s) have no name and get no SMILES.',nEmpty)
end
[uniqueNames, ~, uniqueIdx] = unique(regexprep(model.metNames,'^prot_.*',''));
uniqueSmiles(1:numel(uniqueNames),1) = {''};
metMatch = false(length(uniqueNames),1);
metMatch(strcmp(uniqueNames,'')) = 1; % No need trying to match empty fields
if verbose; fprintf('Check for local SMILES database... '); end
smilesDBfile = (fullfile(params.path,'data','smilesDB.tsv'));
if exist(smilesDBfile,'file')==2
    fID = fopen(smilesDBfile,'r');
    raw = textscan(fID,'%s %s','Delimiter','\t','HeaderLines',0);
    fclose(fID);
    %The last row of a name that occurs more than once is used
    smilesDB.names = flipud(raw{1});
    smilesDB.smile = flipud(raw{2});
    [inDB, metIdx] = ismember(uniqueNames,smilesDB.names);
    uniqueSmiles(inDB) = smilesDB.smile(metIdx(inDB));
    metMatch = metMatch | inDB;
    if verbose; fprintf('done.\n'); end
else
    if verbose; fprintf('not found.\n'); end
end

if any(~metMatch)
    PB = progressReport(numel(uniqueNames), 'Querying PubChem for SMILES by metabolite names');
    webOptions = weboptions('Timeout', 30);
    failedNames = cell(0,1);
    for i = 1:numel(uniqueNames)
        PB.count;
        if metMatch(i)
            continue;
        end
        smileResult = '';
        failed = false;
        for attempt = 1:10
            try
                smileResult = webread(['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' uniqueNames{i} '/property/CanonicalSMILES/TXT'], webOptions);
                %Sometimes multiple lines are given, with alternative SMILES. Only
                %keep the first suggestion.
                smileResult = strtrim(strtok(smileResult, newline));
                failed = false;
                break
            catch exception
                %PubChem answers 404 for a name it does not know and 400 for a
                %name it cannot parse: both are final answers (no match). Any
                %other failure (server busy, error 500, no connection) is retried.
                if strcmp(exception.identifier, 'MATLAB:webservices:HTTP404StatusCodeError') || ...
                        strcmp(exception.identifier, 'MATLAB:webservices:HTTP400StatusCodeError')
                    smileResult = '';
                    failed = false;
                    break
                end
                smileResult = '';
                failed = true;
                pause(0.5*attempt)
            end
        end
        pause(0.2) % PubChem allows at most 5 requests per second
        if failed
            %The request did not get an answer: nothing is cached, so a later
            %call asks again.
            failedNames{end+1,1} = uniqueNames{i}; %#ok<AGROW>
            continue
        end
        uniqueSmiles{i,1} = smileResult;
        % Append one line each time, in case internet connection is lost
        % halfway. Open & close file each time to avoid leaving the file
        % open when breaking the function. A name without a match gets an
        % empty SMILES, so that it is not asked for again.
        out = [uniqueNames(i), uniqueSmiles(i)];
        fID = fopen(smilesDBfile,'a');
        fprintf(fID,'%s\t%s\n',out{:});
        fclose(fID);
    end
    PB.done;
    if ~isempty(failedNames)
        warning('findMetSmiles:requestFailed','%d request(s) to PubChem failed and are not cached; call again to retry.',numel(failedNames))
    end
    if verbose
        fprintf('Model-specific SMILES database stored at %s\n',smilesDBfile);
    end
end
newSmiles = uniqueSmiles(uniqueIdx);
noSMILES = cellfun(@isempty,uniqueSmiles);
successRatio = 1-(numel(find(noSMILES))/numel(uniqueSmiles));
fprintf('SMILES could be found for %s%% of the unique metabolite names.\n',num2str(successRatio*100,'%.0f'))
noSMILES = uniqueNames(noSMILES);

if ~isfield(model,'metSmiles') || all(cellfun(@isempty,model.metSmiles))
    model.metSmiles = newSmiles;
else
    emptySmiles = cellfun(@isempty,model.metSmiles);
    model.metSmiles(emptySmiles) = newSmiles(emptySmiles);
end
end
