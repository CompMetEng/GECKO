function names = infeasibleConditions(model, fluxData, varargin)
% infeasibleConditions  Conditions that are infeasible without enzyme limits.
%
% Checks every condition of fluxData in the model with the protein pool and
% the usage of every enzyme unlimited, which is the model without enzyme
% constraints: growth at the measured rate (or desiredGrowthRate) with every
% measured exchange flux within its band, after the anaerobic switch and the
% protein content of the condition (the makeModelAnaerobic and
% changeProteinBiomass functions of the model adapter). This is the check
% that sensitivityTuning makes first for a condition. No kcat can fix such a
% condition. Give these conditions weight 0 in the flux data, so that no
% tuning is asked to fit them. The model is not changed.
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
% desiredGrowthRate : double
%     growth rate that replaces the measured one in every condition (default
%     []: the measured growth rate).
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
% names : cell
%     the names (fluxData.conds) of the infeasible conditions.
%
% Examples
% --------
%     fluxData = loadFluxData();
%     names = infeasibleConditions(ecModel, fluxData);

p = parseGECKOargs(varargin, { ...
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
params = modelAdapter.params;
bioIdx  = strcmp(model.rxns, params.bioRxn);
poolEx  = strcmp(model.rxns, 'prot_pool_exchange');
if model.ec.geckoLight
    usageRxns = false(numel(model.rxns),1);
else
    usageRxns = startsWith(model.rxns, 'usage_prot_');
end
names = {};
for cond = 1:numel(fluxData.conds)
    [rxnIdx, lb, ub] = conditionBounds(model, fluxData, cond, p.fluxTolerance, p.fluxAbsTolerance, p.fluxTolerances);
    target = fluxData.grRate(cond);
    if ~isempty(p.desiredGrowthRate)
        target = p.desiredGrowthRate;
    end
    mc = setupCondition(model, fluxData, cond, modelAdapter);
    mc.lb(rxnIdx) = lb;
    mc.ub(rxnIdx) = ub;
    mc.lb(bioIdx) = target;
    mc.ub(poolEx) = 1000;
    mc.lb(usageRxns) = 0;
    mc.ub(usageRxns) = 1000;
    mc.c = double(bioIdx);
    res = solveLP(mc,0);
    if isempty(res.f)
        names{end+1,1} = fluxData.conds{cond}; %#ok<AGROW>
    end
end
end

function [rxnIdx, lb, ub] = conditionBounds(m, fluxData, condIdx, tol, absTol, perFlux)
% Bounds that hold every measured exchange of a condition within its band.
rxnIdx = []; lb = []; ub = [];
for i = 1:numel(fluxData.exchRxnIDs)
    v = fluxData.exchFluxes(condIdx,i);
    if isnan(v)
        continue;
    end
    idx = find(strcmp(m.rxns, fluxData.exchRxnIDs{i}));
    if isempty(idx)
        error('infeasibleConditions:missingExchange', ...
            'The exchange reaction %s of the condition is not in the model.', fluxData.exchRxnIDs{i});
    end
    if abs(v) == 1000
        if v < 0
            lo = -1000; hi = 0;
        else
            lo = 0; hi = 1000;
        end
    else
        t = tol;
        if ~isempty(perFlux)
            k = find(strcmp(perFlux(:,1), fluxData.exchRxnIDs{i}));
            if ~isempty(k)
                t = perFlux{k,2};
            end
        end
        band = sort([v*(1-t), v*(1+t)]);
        lo = band(1) - absTol;
        hi = band(2) + absTol;
    end
    rxnIdx(end+1,1) = idx; %#ok<AGROW>
    lb(end+1,1) = lo; %#ok<AGROW>
    ub(end+1,1) = hi; %#ok<AGROW>
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
