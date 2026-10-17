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
% ecRxns : logical
%     for which model.ec.rxns kcats are estimated (default: all).
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

checkNameValue(varargin, {'k','ecRxns','measurements','save','saveAs','modelAdapter'}, 'kcatsFromECTree');
p = parseGECKOargs(varargin, { ...
    'k',            5; ...
    'ecRxns',       []; ...
    'measurements', []; ...
    'save',         true; ...
    'saveAs',       'ec_tree'; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
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
    file = kcatSourceFile(p.saveAs, [], modelAdapter); % fails before estimating, not after
    warnSubset(p.ecRxns, numel(model.ec.rxns), file);
end
measurements = p.measurements;
if isempty(measurements)
    measurements = loadBrendaMeasurements('modelAdapter', modelAdapter);
end

idx = (1:numel(model.ec.rxns))';
if ~isempty(p.ecRxns)
    idx = idx(p.ecRxns(:));
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
    if ~isempty(p.ecRxns)
        parameters.ec_rxns = countRxns(p.ecRxns);
    end
    saveKcatSource(model, kcatList, p.saveAs, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromECTree', 'parameters', parameters, ...
        'data', measurements.release);
    fprintf('EC-tree kcats for %d reactions written to %s\n', numel(kcatList.rxns), file);
end
end
