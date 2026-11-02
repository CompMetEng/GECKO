function fluxData = loadFluxData(varargin)
% loadFluxData  Load total protein measurements and exchange flux data.
%
% Loads total protein content and exchange flux data (e.g. carbon source,
% O2, CO2) from a tab-delimited file.
%
% Name-Value Arguments
% --------------------
% fluxDataFile : char
%     path to file with flux data (default reads data/fluxData.tsv as
%     specified in modelAdapter).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% fluxData : struct
%     structure with flux data.
%
% Notes
% -----
% The function raises an error when the file does not exist
% (loadFluxData:fileNotFound), when a line has a different number of cells
% than the header (loadFluxData:raggedRow) or when the file has no exchange
% flux columns (loadFluxData:noExchangeColumns). An exchange column header has
% the form '<metabolite> (<reaction id>)'; the reaction id may contain
% parentheses.
%
% The fluxData structure has the following fields:
%
% - conds : sampling condition.
% - Ptot : total protein (g/gDCW).
% - grRate : growth rate (1/h).
% - exchFluxes : exchange fluxes (mmol/gDCWh).
% - exchMets : exchanged metabolites, matching exchFluxes.
% - exchRxnIDs : exchange reaction IDs, matching exchMets.
% - evotuneRMSEweight : if column existed in fluxData.tsv, weights for RMSE calculation in evotune kcat tuning.
% - source : if column existed in fluxData.tsv, description where the data comes from.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     fluxData = loadFluxData(fluxDataFile, modelAdapter);
%     fluxData = loadFluxData('fluxDataFile', fluxDataFile);
%
% See also
% --------
% loadProtData

p = parseGECKOargs(varargin, { ...
    'fluxDataFile', []; ...
    'modelAdapter', []});
fluxDataFile = p.fluxDataFile;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(fluxDataFile)
    fluxDataFile = fullfile(params.path,'data','fluxData.tsv');
end

%Load total protein content and flux data
if ~isfile(fluxDataFile)
    error('loadFluxData:fileNotFound', 'The flux data file %s does not exist.', fluxDataFile);
end
fID       = fopen(fluxDataFile);
formatStr = '%s';
data      = textscan(fID,formatStr,'Delimiter','\n');
fclose(fID);
fluxDataRaw  = [];
for i=1:length(data{1})
    row      = data{1}(i);
    row      = strsplit(row{1},'\t');
    row      = row(1:end);
    if i > 1 && numel(row) ~= size(fluxDataRaw,2)
        error('loadFluxData:raggedRow', ...
            'Line %d of %s has %d cells, the header has %d.', ...
            i, fluxDataFile, numel(row), size(fluxDataRaw,2));
    end
    fluxDataRaw = [fluxDataRaw; row]; 
end
if size(fluxDataRaw,2) == 0
    error('The input file %s is likely not tab-delimited.', fluxDataFile)
end
fluxData            = [];

%Find additional fields
extraFieldNames     = {'evotuneRMSEweight','source'};
extraFieldNumeric   = [true,false];      
[logicalExtra, extraFields]    = ismember(extraFieldNames,fluxDataRaw(1,:));
extraFields(~logicalExtra) = [];
extraFieldNames(~logicalExtra) = [];
extraFieldNumeric(~logicalExtra) = [];
if any(extraFields)
    for i=1:numel(extraFields)
        if extraFieldNumeric(i)
            fluxData.(extraFieldNames{i}) = str2double(fluxDataRaw(2:end,extraFields(i)));
        else
            fluxData.(extraFieldNames{i}) = fluxDataRaw(2:end,extraFields(i));
        end
    end
    fluxDataRaw(:,extraFields) = [];
end

%Extract observed byProduct names from file
if size(fluxDataRaw,2) < 4
    error('loadFluxData:noExchangeColumns', ...
        '%s has no exchange flux columns after the columns condition, Ptot and grRate.', fluxDataFile);
end
exchRxns = fluxDataRaw(1,4:end);
exchMets = exchRxns;
for i=1:numel(exchRxns)
    [exchMets{i}, exchRxns{i}] = splitExchHeader(exchRxns{i});
end

fluxData.conds      = fluxDataRaw(2:end,1);
fluxData.Ptot       = str2double(fluxDataRaw(2:end,2));
fluxData.grRate     = str2double(fluxDataRaw(2:end,3));
fluxData.exchFluxes = str2double(fluxDataRaw(2:end,4:end));
%fluxData.exchFluxes(isnan(fluxData.exchFluxes)) = 0;
fluxData.exchMets   = exchMets;
fluxData.exchRxnIDs = exchRxns;
end

function [met, rxn] = splitExchHeader(header)
% Split a column header '<metabolite> (<reaction id>)' at the parenthesis that
% matches the final closing one, so a reaction id may contain parentheses. A
% header of another form gives the whole header as the reaction id.
met = strtrim(header);
rxn = header;
if ~endsWith(header,')')
    return
end
depth = 0;
for k=numel(header):-1:1
    if header(k) == ')'
        depth = depth + 1;
    elseif header(k) == '('
        depth = depth - 1;
        if depth == 0
            if k < numel(header)-1
                met = strtrim(header(1:k-1));
                rxn = header(k+1:end-1);
            end
            return
        end
    end
end
end
