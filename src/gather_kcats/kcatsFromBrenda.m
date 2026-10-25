function kcatList = kcatsFromBrenda(model, varargin)
% kcatsFromBrenda  Gather kcats from BRENDA and keep them as a source file.
%
% Matches the model's EC numbers and substrates to BRENDA by fuzzy
% matching (see fuzzyKcatMatching), and writes the result to
% data/kcats/brenda.tsv in the project folder (see saveKcatSource). The
% matching runs locally on the BRENDA files of GECKO.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure), with EC
%     numbers (getECfromDatabase or getECfromGEM).
%
% Name-Value Arguments
% --------------------
% aggregation : char
%     which BRENDA value per EC number, substrate and organism is used:
%     'max' or 'median' (default: params.kcatAggregation of the model
%     adapter, or 'max').
% forceWClvl : double
%     force a minimum wildcard level in the EC numbers (default 0).
% ecRxns : logical or double
%     for which model.ec.rxns kcats are gathered (default: all). A logical
%     vector or indices of model.ec.rxns, see fuzzyKcatMatching; [] selects
%     all reactions and any other empty selection is an error.
% save : logical
%     write the source file (default true).
% saveAs : char
%     name of the source file, without .tsv (default 'brenda'). With
%     ecRxns, the file then only holds those reactions.
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     the kcat values, as returned by fuzzyKcatMatching.
%
% Examples
% --------
%     kcatsFromBrenda(ecModel);
%     kcatList = kcatsFromBrenda(ecModel, 'aggregation', 'median');
%
% See also
% --------
% fuzzyKcatMatching, kcatsFromDLKcat, kcatsFromOKP, saveKcatSource, loadKcatSource

checkOptions(varargin, {'aggregation','forceWClvl','ecRxns','save','saveAs','modelAdapter'}, 'kcatsFromBrenda');
p = parseGECKOargs(varargin, { ...
    'aggregation',  []; ...
    'forceWClvl',   0; ...
    'ecRxns',       []; ...
    'save',         true; ...
    'saveAs',       'brenda'; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
aggregation = resolveKcatAggregation(p.aggregation, modelAdapter);
if isSubset(p.ecRxns) && countSelected(p.ecRxns) == 0
    error('kcatsFromBrenda:emptySubset', ...
        'ecRxns selects no reactions; leave it out to gather for all reactions.');
end
if p.save
    file = fullfile(modelAdapter.params.path, 'data', 'kcats', [char(p.saveAs) '.tsv']);
    warnSubset(p.ecRxns, numel(model.ec.rxns), file);
end

kcatList = fuzzyKcatMatching(model, 'ecRxns', p.ecRxns, 'modelAdapter', modelAdapter, ...
    'forceWClvl', p.forceWClvl, 'kcatAggregation', aggregation);

if p.save
    parameters = struct('aggregation', aggregation, 'force_wildcard_level', p.forceWClvl);
    if isSubset(p.ecRxns)
        parameters.ec_rxns = numel(kcatList.rxns);
    end
    saveKcatSource(model, kcatList, p.saveAs, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromBrenda', 'parameters', parameters, ...
        'data', brendaRelease(modelAdapter));
    fprintf('BRENDA kcats for %d reactions written to %s\n', sum(kcatList.kcats > 0), file);
end
end

function release = brendaRelease(modelAdapter)
% The release line at the top of kcat.tsv, e.g. 'BRENDA release 2026.1'.
release = '';
file = fullfile(modelAdapter.getBrendaDBFolder(), 'kcat.tsv');
fid = fopen(file, 'r');
if fid < 0
    return
end
line = fgetl(fid);
fclose(fid);
tok = regexp(line, '^#\s*(BRENDA release [^ ]+)', 'tokens', 'once');
if ~isempty(tok)
    release = tok{1};
end
end

function tf = isSubset(ecRxns)
% [] selects all reactions; any other value selects the reactions it names.
tf = ~(isnumeric(ecRxns) && isequal(size(ecRxns),[0 0]));
end

function n = countSelected(ecRxns)
if islogical(ecRxns)
    n = sum(ecRxns);
else
    n = numel(ecRxns);
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
if isSubset(ecRxns) && countRxns(ecRxns) < nRxns && isfile(file)
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
