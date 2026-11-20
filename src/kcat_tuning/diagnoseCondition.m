function conflicts = diagnoseCondition(model, fluxData, varargin)
% diagnoseCondition  Which measured values of a condition conflict with the model.
%
% In the model without enzyme limits (the protein pool and the usage of every
% enzyme unlimited), after the anaerobic switch and the protein content of the
% condition, the bands of the measured exchange fluxes and the growth rate are
% allowed to stretch, at a cost of the stretch relative to the measured value,
% and the total cost is minimized. The measurements that have to move for the
% condition to become feasible are the ones that conflict with the model, with
% its constraints or with each other. A minimal total stretch is not unique in
% general, so the table names one way to make the condition feasible, not
% every way.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure).
% fluxData : struct
%     measured growth rates and exchange fluxes, as read by loadFluxData.
%
% Name-Value Arguments
% --------------------
% condition : double or char
%     index or name of the condition in fluxData.conds (default 1).
% desiredGrowthRate : double
%     growth rate that replaces the measured one (default []: the measured
%     growth rate).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% fluxTolerance : double
%     a measured flux v is met within [v*(1-t), v*(1+t)], widened by
%     fluxAbsTolerance on both sides (default 0.1).
% fluxAbsTolerance : double
%     absolute widening of the band in mmol/gDCW/h (default 1e-3).
% fluxTolerances : cell
%     a two-column cell, reaction ID and tolerance, that replaces
%     fluxTolerance for those reactions (default {}).
%
% Returns
% -------
% conflicts : table
%     one row per measured exchange and one for growth, with the variables
%     reaction, measured, bandLow, bandHigh, achieved (the flux of the least
%     stretched solution), stretch (the distance from the achieved flux to the
%     band, 0 inside it) and relativeStretch (stretch divided by the size of
%     the measured value, where a value below 0.05 counts as 0.05), sorted with
%     the largest relativeStretch first.
%
% Examples
% --------
%     fluxData = loadFluxData();
%     conflicts = diagnoseCondition(ecModel, fluxData, 'condition', 'Std');

p = parseGECKOargs(varargin, { ...
    'condition',         1; ...
    'desiredGrowthRate', []; ...
    'modelAdapter',      []; ...
    'fluxTolerance',     0.1; ...
    'fluxAbsTolerance',  1e-3; ...
    'fluxTolerances',    {}});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
condIdx = conditionIndex(fluxData, p.condition);
bioIdx = find(strcmp(model.rxns, modelAdapter.params.bioRxn));
target = fluxData.grRate(condIdx);
if ~isempty(p.desiredGrowthRate)
    target = p.desiredGrowthRate;
end
sizeFloor = 0.05;

m = setupCondition(model, fluxData, condIdx, modelAdapter);
m.ub(strcmp(m.rxns, 'prot_pool_exchange')) = 1000;
if ~model.ec.geckoLight
    usageRxns = startsWith(m.rxns, 'usage_prot_');
    m.lb(usageRxns) = 0;
    m.ub(usageRxns) = 1000;
end

%the banded reactions: every measured exchange, and growth
idx = []; lo = []; hi = []; measured = [];
for i = 1:numel(fluxData.exchRxnIDs)
    v = fluxData.exchFluxes(condIdx,i);
    if isnan(v) || abs(v) == 1000
        continue;
    end
    k = find(strcmp(m.rxns, fluxData.exchRxnIDs{i}));
    if isempty(k)
        error('diagnoseCondition:missingExchange', ...
            'The exchange reaction %s of the condition is not in the model.', fluxData.exchRxnIDs{i});
    end
    t = p.fluxTolerance;
    if ~isempty(p.fluxTolerances)
        j = find(strcmp(p.fluxTolerances(:,1), fluxData.exchRxnIDs{i}));
        if ~isempty(j)
            t = p.fluxTolerances{j,2};
        end
    end
    band = sort([v*(1-t), v*(1+t)]);
    idx(end+1,1) = k; %#ok<AGROW>
    lo(end+1,1) = band(1) - p.fluxAbsTolerance; %#ok<AGROW>
    hi(end+1,1) = band(2) + p.fluxAbsTolerance; %#ok<AGROW>
    measured(end+1,1) = v; %#ok<AGROW>
end
idx = [idx; bioIdx];
lo = [lo; target];
hi = [hi; Inf];
measured = [measured; target];
sizes = max(abs(measured), sizeFloor);
n = numel(idx);

%v = u + dHigh - dLow, with u inside the band and the stretches d >= 0
metIds = arrayfun(@(k) sprintf('stretchMet_%d',k), (1:n)', 'UniformOutput', false);
metsToAdd = struct();
metsToAdd.mets = metIds;
metsToAdd.metNames = metIds;
metsToAdd.compartments = repmat(m.comps(1), n, 1);
m = addMets(m, metsToAdd);
firstMet = numel(m.mets) - n;
m.c(:) = 0;
rxnIds = {}; rxnMets = {}; rxnCoeffs = {}; rxnLb = []; rxnUb = []; rxnC = [];
for k = 1:n
    m.S(firstMet+k, idx(k)) = 1;
    m.lb(idx(k)) = -1000;
    m.ub(idx(k)) = 1000;
    if idx(k) == bioIdx
        m.lb(idx(k)) = 0;
    end
    upper = min(hi(k), 1000);
    rxnIds{end+1,1} = sprintf('stretchBand_%d',k); rxnMets{end+1,1} = metIds(k); %#ok<AGROW>
    rxnCoeffs{end+1,1} = -1; rxnLb(end+1,1) = lo(k); rxnUb(end+1,1) = upper; rxnC(end+1,1) = 0; %#ok<AGROW>
    if isfinite(hi(k))
        rxnIds{end+1,1} = sprintf('stretchHigh_%d',k); rxnMets{end+1,1} = metIds(k); %#ok<AGROW>
        rxnCoeffs{end+1,1} = -1; rxnLb(end+1,1) = 0; rxnUb(end+1,1) = 1000; rxnC(end+1,1) = -1/sizes(k); %#ok<AGROW>
    end
    rxnIds{end+1,1} = sprintf('stretchLow_%d',k); rxnMets{end+1,1} = metIds(k); %#ok<AGROW>
    rxnCoeffs{end+1,1} = 1; rxnLb(end+1,1) = 0; rxnUb(end+1,1) = 1000; rxnC(end+1,1) = -1/sizes(k); %#ok<AGROW>
end
rxnsToAdd = struct();
rxnsToAdd.rxns = rxnIds;
rxnsToAdd.rxnNames = rxnIds;
rxnsToAdd.mets = rxnMets;
rxnsToAdd.stoichCoeffs = rxnCoeffs;
rxnsToAdd.lb = rxnLb;
rxnsToAdd.ub = rxnUb;
rxnsToAdd.c = rxnC;
m = addRxns(m, rxnsToAdd, 1);

res = solveLP(m, 0);
if isempty(res.f)
    error('diagnoseCondition:blocked', ...
        'The condition is infeasible even with every band stretched; the model itself blocks it.');
end
achieved = res.x(idx);
stretch = max([lo - achieved, zeros(n,1), achieved - hi], [], 2);
reaction = m.rxns(idx);
relativeStretch = stretch ./ sizes;
conflicts = table(reaction, measured, lo, hi, achieved, stretch, relativeStretch, ...
    'VariableNames', {'reaction','measured','bandLow','bandHigh','achieved','stretch','relativeStretch'});
conflicts = sortrows(conflicts, 'relativeStretch', 'descend');
end

function idx = conditionIndex(fluxData, condition)
% Index of a condition given by its name or its row.
if ischar(condition) || isstring(condition)
    idx = find(strcmp(fluxData.conds, char(condition)));
    if isempty(idx)
        error('diagnoseCondition:unknownCondition', ...
            'Condition %s not found in fluxData.conds.', char(condition));
    end
else
    idx = condition;
    if idx < 1 || idx > numel(fluxData.conds)
        error('diagnoseCondition:unknownCondition', ...
            'Condition index %d is out of range (%d conditions).', idx, numel(fluxData.conds));
    end
end
end

function m = setupCondition(m, fluxData, condIdx, modelAdapter)
% Apply the anaerobic switch and the protein content of a condition, through
% the functions of the model adapter.
o2 = find(strcmp(fluxData.exchMets, 'oxygen'));
if ~isempty(o2) && fluxData.exchFluxes(condIdx,o2(1)) == 0
    m = modelAdapter.makeModelAnaerobic(m);
end
ptot = fluxData.Ptot(condIdx);
if ~(isnan(ptot) || ptot == 0)
    m = modelAdapter.changeProteinBiomass(m, ptot);
end
end
