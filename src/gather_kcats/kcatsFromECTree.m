function kcatList = kcatsFromECTree(model, varargin)
% kcatsFromECTree  Estimate kcats from the EC tree and keep them as a source file.
%
% Each reaction gets the median kcat of all wild-type BRENDA measurements
% of its EC number, from any organism and substrate, shrunk toward the
% median of the EC's 3-digit and 2-digit classes with weight n / (n + k)
% for n measurements (see ecTreeEstimate). A reaction with several EC
% numbers uses the one with most measurements. Reactions without an EC
% number, or whose EC classes have no measurements, get no value. The
% result is written to data/kcats/ec_tree.tsv in the project folder (see
% saveKcatSource).
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure), with EC
%     numbers (getECfromDatabase or getECfromGEM).
%
% Name-Value Arguments
% --------------------
% k : double
%     shrinkage constant (default 5): the number of measurements at which
%     an EC's own median and its class count equally.
% ecRxns : logical or double
%     for which model.ec.rxns kcats are estimated, as a logical vector or as
%     indices (default: all). A reaction named more than once gives one
%     row. An empty ecRxns given explicitly selects no reaction.
% measurements : struct
%     BRENDA measurement table (default: loadBrendaMeasurements with the
%     model adapter).
% save : logical
%     write the source file (default true).
% saveAs : char
%     name of the source file, without .tsv (default 'ec_tree'). With
%     ecRxns, the file then only holds those reactions.
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     the kcat values, with fields source, rxns, eccodes, kcats and detail
%     (the number of measurements of the EC, n, and the deepest EC level
%     with measurements, depth).
%
% Examples
% --------
%     kcatsFromECTree(ecModel);
%     kcatList = kcatsFromECTree(ecModel, 'k', 10, 'save', false);
%
% See also
% --------
% ecTreeEstimate, loadBrendaMeasurements, kcatsFromBrenda, saveKcatSource, assignKcats

checkOptions(varargin, {'k','ecRxns','measurements','save','saveAs','modelAdapter'}, 'kcatsFromECTree');
p = parseGECKOargs(varargin, { ...
    'k',            5; ...
    'ecRxns',       []; ...
    'measurements', []; ...
    'save',         true; ...
    'saveAs',       'ec_tree'; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
hasSubset = any(strcmp(cellfun(@char, varargin(1:2:end), 'UniformOutput', false), 'ecRxns'));
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
if ~isfield(model.ec, 'eccodes')
    error('No EC codes defined in model.ec.eccodes. First run getECfromGEM() and/or getECfromDatabase().')
end
if p.save
    file = fullfile(modelAdapter.params.path, 'data', 'kcats', [char(p.saveAs) '.tsv']);
    warnSubset(p.ecRxns, hasSubset, numel(model.ec.rxns), file);
end
measurements = p.measurements;
if isempty(measurements)
    measurements = loadBrendaMeasurements('modelAdapter', modelAdapter);
end

idx = (1:numel(model.ec.rxns))';
if hasSubset
    selected = false(numel(idx), 1);
    selected(p.ecRxns(:)) = true;
    idx = idx(selected);
end
eccodes = model.ec.eccodes(idx);
eccodes(cellfun(@(e) ~ischar(e) && ~isstring(e), eccodes)) = {''};
[values, n, depth, ecUsed] = ecTreeEstimate(eccodes, measurements, p.k);
found = depth > 0;
idx = idx(found);
kcatList.source  = 'ec_tree';
kcatList.rxns    = model.ec.rxns(idx);
kcatList.eccodes = ecUsed(found);
kcatList.kcats   = 10 .^ values(found);
kcatList.detail  = cellstr("n=" + string(n(found)) + "; depth=" + string(depth(found)));

if p.save
    parameters = struct('k', p.k);
    if hasSubset
        parameters.ec_rxns = countRxns(p.ecRxns);
    end
    saveKcatSource(model, kcatList, p.saveAs, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromECTree', 'parameters', parameters, ...
        'data', measurements.release);
    fprintf('EC-tree kcats for %d reactions written to %s\n', numel(kcatList.rxns), file);
end
end

function n = countRxns(ecRxns)
if islogical(ecRxns)
    n = sum(ecRxns);
else
    n = numel(unique(ecRxns));
end
end

function warnSubset(ecRxns, hasSubset, nRxns, file)
% A subset written under the default name replaces the full file.
if hasSubset && countRxns(ecRxns) < nRxns && isfile(file)
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
