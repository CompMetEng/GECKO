function kcatList = readDLKcatOutput(model, varargin)
% readDLKcatOutput  Read the DLKcat output file into a kcatList structure.
%
% Reads the DLKcat output file and constructs a kcatList structure, that can
% be used by assignKcatValues() to populate the ecModel with kcat values.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% outFile : char
%     name and path of the DLKcat output file (default: data/DLKcat.tsv from
%     the obj.params.path folder specified in the modelAdapter).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% strict : logical
%     raise an error for a substrate that is not in model.metNames
%     (readDLKcatOutput:unknownSubstrate). By default the rows with such a
%     substrate are left out with a warning (default false).
%
% Rows without a numeric kcat ('NA', 'None', empty, text, 'inf') are left out
% before the reactions are checked. A file without any such row raises
% readDLKcatOutput:noKcats. A reaction that is not in model.ec.rxns raises
% readDLKcatOutput:unknownReaction.
%
% Returns
% -------
% kcatList : struct
%     structure array with list of DLKcat derived kcat values, with separate
%     entries for each kcat value.
%
% Notes
% -----
% The kcatList structure has the following fields:
%
% - source : 'DLKcat'.
% - rxns : reaction identifiers.
% - genes : gene identifiers.
% - substrates : substrate names.
% - kcats : predicted kcat value in /sec.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     kcatList = readDLKcatOutput(model);
%     kcatList = readDLKcatOutput(model, 'outFile', 'myDLKcat.tsv');
%
% See also
% --------
% runDLKcat, writeDLKcatInput, assignKcatValues

p = parseGECKOargs(varargin, { ...
    'outFile',      []; ...
    'modelAdapter', []; ...
    'strict',       false});
outFile      = p.outFile;
modelAdapter = p.modelAdapter;
strict       = p.strict;

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
    error('readDLKcatOutput:noKcats', 'DLKcat file does not contain any kcat values, please run runDLKcat() first.')
end

% Check that the substrates are in the model (case-insensitive: an SBML loader can
% normalise metabolite-name capitalisation differently than whatever produced the
% DLKcat input, and that must not read as a genuinely missing substrate). Rows with
% another substrate are left out with a warning, or raise an error with 'strict'.
matchMets = ismember(lower(subs),lower(model.metNames));
if ~all(matchMets)
    errorText = 'DLKcat was likely run with an input file that was generated from another ecModel, as the following substrates from DLKcat output cannot be found in model.metNames:';
    if strict
        error('readDLKcatOutput:unknownSubstrate', '%s', ravenList(errorText,subs(~matchMets),false))
    end
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
