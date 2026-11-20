function [model, tunedKcats] = sensitivityTuning(model, varargin)
% sensitivityTuning  Relax the most limiting kcats to reach a growth rate.
%
% Iteratively increases the kcat of whichever enzyme currently limits growth
% the most, until the model reaches the desired growth rate, updating the
% tuned values in model.ec.kcat. The loop also stops, with a warning, when an
% iteration raises growth by less than 1e-9.
%
% With fluxData, the target is a measured condition instead of a growth rate:
% growth at the measured rate, and every measured exchange flux within
% fluxTolerance of its value. Each iteration checks the condition; when it is
% infeasible, the protein pool is lifted as a limit, the protein drawn from it
% is minimized and the kcat of the enzyme that draws most is raised. A
% condition that is infeasible without any enzyme limit, which is the model
% without enzyme constraints, is not something a kcat can fix (a measurement
% that conflicts with the model): the tuning stops before any change with
% status 'infeasibleInGem', and the condition is better left out of the tuning
% (give it weight 0 in the flux data, see infeasibleConditions). A condition
% that is feasible without enzyme limits but not with the pool lifted is
% stopped by a limit of single enzymes, status 'unfixable'.
%
% This is the standard kcat tuning: it needs one maximum growth rate.
% cmaesKcatTuning is the advanced tuning, which fits kcats to maximum
% growth rates and exchange fluxes, measured in one or more conditions.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% desiredGrowthRate : double
%     kcats will be relaxed until this growth rate is reached (default: the
%     experimental growth rate gr_exp from the model adapter).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% foldChange : double
%     kcat values will be increased by this fold-change (default 10).
% protToIgnore : cell
%     vector of protein ids to be ignored in tuned kcats, e.g. {'P38122',
%     'Q99271'} (default []).
% verbose : logical
%     whether progress should be reported (default true).
% fluxData : struct
%     measured growth rate and exchange fluxes, as read by loadFluxData. The
%     model's own medium is kept apart from the measured exchanges, which
%     take the bounds of the condition (an exchange of -1000 or 1000 stays
%     unconstrained in that direction, and NaN is not measured).
% condition : double or char
%     index or name of the condition in fluxData.conds (default 1). The
%     measured growth rate is the target unless desiredGrowthRate is given.
% fluxTolerance : double
%     a measured flux v is met within [v*(1-t), v*(1+t)], widened by
%     fluxAbsTolerance on both sides (default 0.1).
% fluxAbsTolerance : double
%     absolute widening of the band in mmol/gDCW/h, which fluxes near zero
%     need (default 1e-3).
% fluxTolerances : cell
%     a two-column cell, reaction ID and tolerance, that replaces
%     fluxTolerance for those reactions (default {}).
%
% Returns
% -------
% model : struct
%     ecModel with updated model.ec.kcat.
% tunedKcats : struct
%     structure with information on tuned kcat values.
%
% Notes
% -----
% The tunedKcats structure has the following fields:
%
% - rxns : identifiers of reactions with tuned kcat values.
% - rxnNames : names of the reactions in tunedKcats.rxns.
% - enzymes : enzymes that catalyze the reactions in tunedKcats.rxns, whose
%   kcat value has been tuned.
% - oldKcat : kcat values in the input model.
% - newKcat : kcat values in the output model, after tuning.
% - source : source label of each kcat, taken from the input model before
%   tuning.
% - status : why the tuning stopped: 'reached', 'infeasibleInGem', 'unfixable',
%   'stalled', 'infeasible', 'noCandidate' or 'maxIterations'.
%
% model.ec.source of a tuned kcat keeps its origin and detail and gains the
% tuning status, e.g. 'brenda (wc=0, origin=1) | tuned=sensitivity
% (prior=0.05)', where prior is the kcat before the first tuning (see
% formatKcatSource). An origin that is missing is written as 'unknown'.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [model, tunedKcats] = sensitivityTuning(model);
%     [model, tunedKcats] = sensitivityTuning(model, 'foldChange', 10);
%     fluxData = loadFluxData();
%     [model, tunedKcats] = sensitivityTuning(model, 'fluxData', fluxData, 'condition', 'Std');
%     [model, tunedKcats] = sensitivityTuning(model, desiredGrowthRate, modelAdapter, foldChange, protToIgnore, verbose);

p = parseGECKOargs(varargin, { ...
    'desiredGrowthRate', []; ...
    'modelAdapter',      []; ...
    'foldChange',        10; ...
    'protToIgnore',      {}; ...
    'verbose',           true; ...
    'fluxData',          []; ...
    'condition',         1; ...
    'fluxTolerance',     0.1; ...
    'fluxAbsTolerance',  1e-3; ...
    'fluxTolerances',    {}});
desiredGrowthRate = p.desiredGrowthRate;
modelAdapter      = p.modelAdapter;
foldChange        = p.foldChange;
protToIgnore      = p.protToIgnore;
verbose           = p.verbose;
fluxData          = p.fluxData;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
if ~isempty(fluxData)
    condIdx = conditionIndex(fluxData, p.condition);
    if isempty(desiredGrowthRate)
        desiredGrowthRate = fluxData.grRate(condIdx);
    end
end
if isempty(desiredGrowthRate)
    if ~isfield(params, 'gR_exp') || isempty(params.gR_exp)
        error('GECKO:grExpNotSet', 'gr_exp is not set in model_adapter.toml')
    end
    desiredGrowthRate = params.gR_exp;
end

kcatList = [];
m = model;
m.c = double(strcmp(m.rxns, params.bioRxn));% Make sure that growth is maximized

[res,hs] = solveLP(m);
if isempty(res.x)
    error('FBA of input model gives no valid result. Reduce protein pool constraint with setProtPoolSize and/or check if exchange constraints are correctly defined.')
end
% The loop stops when growth rises by less than growthStallTol over the
% previous iteration (the first iteration has nothing to compare with).
lastGrowth = NaN;
growthStallTol = 1e-9;
status = 'maxIterations';
if ~isempty(fluxData)
    [m, kcatList, status] = tuneToCondition(m, fluxData, condIdx, desiredGrowthRate, ...
        p.fluxTolerance, p.fluxAbsTolerance, p.fluxTolerances, foldChange, protToIgnore, ...
        params, verbose);
elseif ~m.ec.geckoLight
    %for the full model, we first find the draw reaction with the most flux
    drawRxns = startsWith(m.rxns, 'usage_prot_');
    idxToIgnore = cellfun(@(x) find(strcmpi(model.rxns, ['usage_prot_' x])), protToIgnore);
    iteration = 1;
    while true
        [res,hs] = solveLP(m,0,[],hs); %skip parsimonius, only takes time
        if isempty(res.f)
            printOrange('WARNING: The model has no feasible solution with the increased kcats, tuning stops.\n');
            status = 'infeasible';
            break;
        end
        if ~isnan(lastGrowth) && res.f <= lastGrowth + growthStallTol
            printOrange('WARNING: No growth increase from increased kcats - check if the constraints on the uptake reactions are too tight.\n');
            status = 'stalled';
            break;
        end
        lastGrowth = res.f;
        if verbose; disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]); end
        if (lastGrowth >= desiredGrowthRate)
            status = 'reached';
            break;
        end
        %If you get an error here, it is likely due to numerical issues in the solver
        %The trick where we don't allow low kcats is to fix that, but maybe
        %it is not enough.
        iteration            = iteration + 1;
        %find the highest draw_prot rxn flux
        drawFluxes           = zeros(length(drawRxns),1);
        drawFluxes(drawRxns) = res.x(drawRxns);
        % Remove from the list user defined proteins
        drawFluxes(idxToIgnore) = 0;
        if ~any(drawFluxes > 0)
            printOrange('WARNING: No protein usage reaction carries flux, tuning cannot continue.\n');
            status = 'noCandidate';
            break;
        end
        [~,sel]              = max(drawFluxes); % since bounds 0 to 1000
        %Now get the metabolite
        metSel               = m.S(:,sel) > 0; % positive coeff (forward usage)
        %now find the reaction with the largest consumption of this protein
        protFluxes           = m.S(metSel,:).' .* res.x; %negative
        [~,rxnSel]           = min(protFluxes);
        kcatList             = [kcatList, rxnSel];
        rxn                  = m.rxns(rxnSel);
        targetSubRxn         = strcmp(m.ec.rxns, rxn);
        m                    = markTuned(m,targetSubRxn);
        m.ec.kcat(targetSubRxn)     = m.ec.kcat(targetSubRxn) .* foldChange;
        m                    = applyKcatConstraints(m,targetSubRxn);
    end

else
    origRxns = extractAfter(m.ec.rxns,4);
    %find the reactions involved in proteins to be ignored
    idxToIgnore = cellfun(@(x) find(m.ec.rxnEnzMat(:, strcmpi(m.ec.enzymes, x))), protToIgnore, 'UniformOutput', false);
    %create an unique vector
    idxToIgnore = unique(cat(1, idxToIgnore{:}));
    %get the correct idx in model.rxns
    idxToIgnore = cellfun(@(x) find(strcmpi(m.rxns, x)), origRxns(idxToIgnore));
    iteration = 1;
    while true
        res = solveLP(m,0); %skip parsimonius, only takes time
        if isempty(res.f)
            printOrange('WARNING: The model has no feasible solution with the increased kcats, tuning stops.\n');
            status = 'infeasible';
            break;
        end
        if ~isnan(lastGrowth) && res.f <= lastGrowth + growthStallTol
            printOrange('WARNING: No growth increase from increased kcats - check if the constraints on the uptake reactions are too tight.\n');
            status = 'stalled';
            break;
        end
        lastGrowth = res.f;
        if verbose; disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]); end
        if (lastGrowth >= desiredGrowthRate)
            status = 'reached';
            break;
        end
        %If you get an error here, it is likely due to numerical issues in the solver
        %The trick where we don't allow low kcats is to fix that, but maybe
        %it is not enough.
        iteration       = iteration + 1;
        %find the highest protein usage flux
        protPoolStoich  = m.S(strcmp(m.mets, 'prot_pool'),:).';
        protPoolStoich(idxToIgnore) = 0;
        if ~any(res.x .* protPoolStoich < 0)
            printOrange('WARNING: No reaction draws on the protein pool, tuning cannot continue.\n');
            status = 'noCandidate';
            break;
        end
        [~,sel]         = min(res.x .* protPoolStoich); %max consumption
        kcatList        = [kcatList, sel];
        rxn             = m.rxns(sel.');
        targetSubRxns   = strcmp(origRxns, rxn);
        m               = markTuned(m,targetSubRxns);
        m.ec.kcat(targetSubRxns) = m.ec.kcat(targetSubRxns) .* foldChange;
        m               = applyKcatConstraints(m,rxn);
    end
end

kcatList        = unique(kcatList);
tunedKcats.rxns = m.rxns(kcatList);
tunedKcats.rxnNames = m.rxnNames(kcatList);
if ~model.ec.geckoLight
    [~, rxnIdx]     = ismember(tunedKcats.rxns,m.ec.rxns);
else
    [~, rxnIdx]     = ismember(tunedKcats.rxns,origRxns);
end
tunedKcats.enzymes  = cell(numel(kcatList),1);
for i=1:numel(rxnIdx)
    [~, metIdx]     = find(m.ec.rxnEnzMat(rxnIdx(i),:));
    tunedKcats.enzymes{i}   = strjoin(m.ec.enzymes(metIdx),';');
end
tunedKcats.oldKcat  = model.ec.kcat(rxnIdx);
tunedKcats.newKcat  = m.ec.kcat(rxnIdx);
tunedKcats.source   = model.ec.source(rxnIdx);
tunedKcats.status   = status;

model = m;
end

function m = markTuned(m,rows)
% Add the tuning method sensitivity to the source of the ec rows, with the
% kcat before the first tuning as prior. Called before the kcat changes.
for i = find(rows(:))'
    info = parseKcatSource(m.ec.source{i});
    origin = info.origin;
    if isempty(origin)
        origin = 'unknown';
    end
    tuned = info.tuned;
    if ~any(strcmp(tuned,'sensitivity'))
        tuned{end+1} = 'sensitivity';
    end
    prior = info.prior;
    if isnan(prior)
        prior = m.ec.kcat(i);
    end
    m.ec.source{i} = formatKcatSource(origin,info.detail,tuned,prior);
end
end

function idx = conditionIndex(fluxData, condition)
% Index of a condition given by its name or its row.
if ischar(condition) || isstring(condition)
    idx = find(strcmp(fluxData.conds, char(condition)));
    if isempty(idx)
        error('sensitivityTuning:unknownCondition', ...
            'Condition %s not found in fluxData.conds.', char(condition));
    end
else
    idx = condition;
    if idx < 1 || idx > numel(fluxData.conds)
        error('sensitivityTuning:unknownCondition', ...
            'Condition index %d is out of range (%d conditions).', idx, numel(fluxData.conds));
    end
end
end

function [m, kcatList, status] = tuneToCondition(m, fluxData, condIdx, targetGrowth, tol, absTol, perFlux, foldChange, protToIgnore, params, verbose)
% Raise kcats until the model meets growth and the measured fluxes of a condition.
kcatList = [];
status = 'maxIterations';
maxIterations = 1000;
bioIdx  = strcmp(m.rxns, params.bioRxn);
poolEx  = strcmp(m.rxns, 'prot_pool_exchange');
poolMet = strcmp(m.mets, 'prot_pool');
if ~any(poolEx)
    error('sensitivityTuning:noPoolExchange', ...
        'The model has no prot_pool_exchange reaction, which the tuning to a condition needs.');
end
[rxnIdx, lb, ub] = conditionBounds(m, fluxData, condIdx, tol, absTol, perFlux);
light = m.ec.geckoLight;
if ~light
    drawRxns = startsWith(m.rxns, 'usage_prot_');
    idxToIgnore = cellfun(@(x) find(strcmpi(m.rxns, ['usage_prot_' x])), protToIgnore);
else
    origRxns = extractAfter(m.ec.rxns,4);
    idxToIgnore = cellfun(@(x) find(m.ec.rxnEnzMat(:, strcmpi(m.ec.enzymes, x))), protToIgnore, 'UniformOutput', false);
    idxToIgnore = unique(cat(1, idxToIgnore{:}));
    idxToIgnore = cellfun(@(x) find(strcmpi(m.rxns, x)), origRxns(idxToIgnore));
end
if ~light
    usageRxns = startsWith(m.rxns, 'usage_prot_');
else
    usageRxns = false(numel(m.rxns),1);
end
mc = m;
mc.lb(rxnIdx) = lb;
mc.ub(rxnIdx) = ub;
mc.lb(bioIdx) = targetGrowth;
mc.ub(poolEx) = 1000;
mc.lb(usageRxns) = 0;
mc.ub(usageRxns) = 1000;
mc.c = double(bioIdx);
res = solveLP(mc,0);
if isempty(res.f)
    printOrange('WARNING: The condition is infeasible without any enzyme limit, which is the model without enzyme constraints, so no kcat can meet it; check the measured fluxes against the model.\n');
    status = 'infeasibleInGem';
    return;
end
lastPool = NaN;
for iteration = 1:maxIterations
    mc = m;
    mc.lb(rxnIdx) = lb;
    mc.ub(rxnIdx) = ub;
    mc.lb(bioIdx) = targetGrowth;
    mc.c = double(bioIdx);
    res = solveLP(mc,0);
    if ~isempty(res.f)
        if verbose; disp(['The condition is met after ' num2str(iteration) ' iteration(s).']); end
        status = 'reached';
        return;
    end
    mc.ub(poolEx) = 1000;
    mc.c = -double(poolEx);
    res = solveLP(mc,0);
    if isempty(res.f)
        printOrange('WARNING: The condition is feasible without enzyme limits but infeasible with the protein pool lifted, so a limit of single enzymes stops it.\n');
        status = 'unfixable';
        return;
    end
    poolUsed = -res.f;
    if ~isnan(lastPool) && poolUsed >= lastPool - 1e-9
        printOrange('WARNING: The protein needed for the condition did not fall, tuning stops.\n');
        status = 'stalled';
        return;
    end
    lastPool = poolUsed;
    if verbose; disp(['Iteration ' num2str(iteration) ': the condition needs ' num2str(poolUsed) ' of protein']); end
    if ~light
        draw = zeros(numel(m.rxns),1);
        draw(drawRxns) = res.x(drawRxns) .* abs(full(m.S(poolMet, drawRxns)).');
        draw(idxToIgnore) = 0;
        if ~any(draw > 0)
            printOrange('WARNING: No protein usage reaction carries flux, tuning cannot continue.\n');
            status = 'noCandidate';
            return;
        end
        [~,sel]      = max(draw);
        metSel       = m.S(:,sel) > 0;
        protFluxes   = m.S(metSel,:).' .* res.x;
        [~,rxnSel]   = min(protFluxes);
        kcatList     = [kcatList, rxnSel];
        targetSubRxn = strcmp(m.ec.rxns, m.rxns(rxnSel));
        m            = markTuned(m,targetSubRxn);
        m.ec.kcat(targetSubRxn) = m.ec.kcat(targetSubRxn) .* foldChange;
        m            = applyKcatConstraints(m,targetSubRxn);
    else
        protPoolStoich = m.S(poolMet,:).';
        protPoolStoich(idxToIgnore) = 0;
        if ~any(res.x .* protPoolStoich < 0)
            printOrange('WARNING: No reaction draws on the protein pool, tuning cannot continue.\n');
            status = 'noCandidate';
            return;
        end
        [~,sel]         = min(res.x .* protPoolStoich);
        kcatList        = [kcatList, sel];
        rxn             = m.rxns(sel.');
        targetSubRxns   = strcmp(origRxns, rxn);
        m               = markTuned(m,targetSubRxns);
        m.ec.kcat(targetSubRxns) = m.ec.kcat(targetSubRxns) .* foldChange;
        m               = applyKcatConstraints(m,rxn);
    end
end
printOrange('WARNING: The condition was not met within the maximum number of iterations, tuning stops.\n');
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
        error('sensitivityTuning:missingExchange', ...
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
