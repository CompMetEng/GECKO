function [model, flexEnz] = relaxEnzConcs(model, varargin)
% relaxEnzConcs  Relax enzyme concentrations by the least total amount.
%
% Relaxes the upper bounds of the usage reactions of the measured enzymes
% (and, where needed, of the protein pool exchange) by the least total
% amount that lets the model reach the target growth rate, in one linear
% program. flexibilizeEnzConcs and relaxProteomicsGreedy relax one enzyme at
% a time, so which enzymes end up relaxed depends on how the most limiting
% enzyme is ranked; relaxEnzConcs does not depend on such a ranking.
%
% The upper bound of the usage reaction of every measured enzyme becomes
% base * (1 + x), where base is its current upper bound (at most 1000), and
% the sum of the folds x is minimised with the growth rate at least
% expGrowth. The protein pool exchange may be raised by a fold as well, at
% poolWeight times the price of the fold of an enzyme, so that it is raised
% only where raising enzymes does not suffice. A current upper bound is
% never lowered. The concentrations in model.ec.concs are not changed.
%
% Parameters
% ----------
% model : struct
%     a full ecModel with model.ec.concs filled and the concentrations
%     applied as upper bounds (fillEnzConcs and constrainEnzConcs).
%
% Name-Value Arguments
% --------------------
% expGrowth : double
%     the growth rate to reach (default: gR_exp of the model adapter).
% bioRxn : char
%     the biomass reaction (default: bioRxn of the model adapter).
% perComplex : logical
%     whether the measured subunits of a complex (the subunits it owns, see
%     getComplexFamilies) share one fold, so that concentrations levelled
%     by levelComplexConcs stay levelled. A subunit with an upper bound of 0
%     gets its own fold (default false).
% weight : char
%     'relative' minimises the sum of the folds, so that an enzyme measured
%     at almost nothing is not relaxed freely because its absolute cost is
%     small; 'absolute' minimises the mg/gDCW added (default 'relative').
% poolWeight : double
%     price of one fold of the protein pool relative to one fold of an
%     enzyme, or per mg with weight 'absolute' (default 1000).
% floor : double
%     a base below this value is taken as this value when scaling its fold,
%     so that an enzyme measured at 0 can still be relaxed, at a high
%     relative price (default 1e-9).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% model : struct
%     the ecModel with the relaxed upper bounds of the usage reactions and
%     of prot_pool_exchange. The upper bound of the biomass reaction is
%     raised to expGrowth when it was lower.
% flexEnz : struct
%     the relaxed enzymes, with the fields of flexibilizeEnzConcs:
%
%     - uniprotIDs : enzymes whose usage upper bound was relaxed, and
%       'prot_pool' when the protein pool exchange was raised.
%     - oldConcs : the upper bounds before relaxation.
%     - flexConcs : the relaxed upper bounds.
%     - ratioIncr : flexConcs ./ oldConcs; the entries are sorted by this
%       field, highest first.
%     - frequence : 1 for every entry.
%
% Notes
% -----
% The error relaxEnzConcs:noMeasured is raised when no enzyme has a
% concentration, relaxEnzConcs:infeasible when no relaxation of the
% measured enzymes and the protein pool reaches expGrowth (the growth rate
% is then not limited by protein), and relaxEnzConcs:arguments for an
% unknown weight. A fold of at most 1e-9 counts as no relaxation. geckopy
% has the same function (relax_enz_concs).
%
% Examples
% --------
%     [model, flexEnz] = relaxEnzConcs(model);
%     [model, flexEnz] = relaxEnzConcs(model, 'expGrowth', 0.3, 'perComplex', true);
%
% See also
% --------
% flexibilizeEnzConcs, relaxProteomicsGreedy, constrainEnzConcs, levelComplexConcs

p = parseGECKOargs(varargin, { ...
    'expGrowth',    []; ...
    'bioRxn',       []; ...
    'perComplex',   false; ...
    'weight',       'relative'; ...
    'poolWeight',   1000; ...
    'floor',        1e-9; ...
    'modelAdapter', []});
weight     = char(p.weight);
poolWeight = p.poolWeight;
floorConc  = p.floor;
if ~any(strcmp(weight, {'relative', 'absolute'}))
    error('relaxEnzConcs:arguments', 'weight must be ''relative'' or ''absolute'', not ''%s''.', weight)
end
expGrowth = p.expGrowth;
bioRxn    = p.bioRxn;
if isempty(expGrowth) || isempty(bioRxn)
    modelAdapter = p.modelAdapter;
    if isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end
    params = modelAdapter.getParameters();
    if isempty(expGrowth)
        expGrowth = params.gR_exp;
    end
    if isempty(bioRxn)
        bioRxn = params.bioRxn;
    end
end
bioIdx  = getIndexes(model, bioRxn, 'rxns');
poolIdx = find(strcmp(model.rxns, 'prot_pool_exchange'));
if isempty(poolIdx)
    error('The model has no prot_pool_exchange reaction; makeEcModel adds it.')
end

measured = find(isfinite(model.ec.concs(:)));
if isempty(measured)
    error('relaxEnzConcs:noMeasured', ['No enzyme concentrations are measured. ' ...
        'Run fillEnzConcs and constrainEnzConcs first.'])
end
[found, usageIdx] = ismember(strcat('usage_prot_', model.ec.enzymes(measured)), model.rxns);
if ~all(found)
    error('Usage reactions are not defined for all enzymes. This is done by makeEcModel.')
end
oldCaps = model.ub(usageIdx);
% Folds count from the current upper bound, so a bound raised before is never lowered
base = min(oldCaps, 1000);
scale = max(base, floorConc);

% Units share one fold: single enzymes, or the measured subunits of a complex
nMeasured = numel(measured);
unitOf = zeros(nMeasured, 1);
nUnits = 0;
if p.perComplex
    families = getComplexFamilies(model);
end
for i = 1:nMeasured
    if unitOf(i) > 0
        continue
    end
    nUnits = nUnits + 1;
    members = i;
    if p.perComplex && base(i) > 0 && families.owner(measured(i)) > 0
        cplx = families.complexes{families.owner(measured(i))};
        [~, owned] = ismember(cplx.enzymes(cplx.owned), model.ec.enzymes(measured));
        owned = owned(owned > 0);
        members = owned(unitOf(owned) == 0 & base(owned) > 0);
    end
    unitOf(members) = nUnits;
end

% The linear program: a copy of the model with one fold variable per unit
% and one for the pool, and one row per relaxed upper bound:
% usage - scale * x <= base, and pool - poolCap * xPool <= poolCap
nRxns = numel(model.rxns);
nMets = numel(model.mets);
poolCap = model.ub(poolIdx);
lp = model;
lp.ub(usageIdx) = 1000;
lp.ub(poolIdx) = 1000;
lp.lb(bioIdx) = expGrowth;
lp.ub(bioIdx) = max(model.ub(bioIdx), expGrowth);
foldCols = nRxns + (1:nUnits)';
poolCol = nRxns + nUnits + 1;
rows = sparse([1:nMeasured, 1:nMeasured, nMeasured+1, nMeasured+1], ...
    [usageIdx', foldCols(unitOf)', poolIdx, poolCol], ...
    [ones(1, nMeasured), -scale', 1, -poolCap], nMeasured + 1, poolCol);
lp.S = [lp.S, sparse(nMets, nUnits + 1); rows];
lp.lb = [lp.lb; zeros(nUnits + 1, 1)];
lp.ub = [lp.ub; inf(nUnits + 1, 1)];
lp.rxns = [lp.rxns; arrayfun(@(k) sprintf('relax_unit_%d', k), (1:nUnits)', 'UniformOutput', false); ...
    {'relax_pool'}];
lp.mets = [lp.mets; strcat('relax_cap_', model.ec.enzymes(measured)); {'relax_cap_pool'}];
if size(lp.b, 2) == 1
    lp.b = [lp.b lp.b];
end
lp.b = [lp.b; -inf(nMeasured + 1, 1), [base; poolCap]];
% solveLP maximises, so the costs are negative
lp.c = zeros(poolCol, 1);
if strcmp(weight, 'relative')
    lp.c(foldCols) = -1;
    lp.c(poolCol) = -poolWeight;
else
    lp.c(foldCols) = -accumarray(unitOf, scale, [nUnits 1]);
    lp.c(poolCol) = -poolWeight * poolCap;
end

sol = solveLP(lp);
if sol.stat ~= 1
    error('relaxEnzConcs:infeasible', ['No relaxation of the measured enzymes or the protein ' ...
        'pool reaches growth %g; the target is not limited by protein.'], expGrowth)
end

folds = sol.x(foldCols);
relaxed = folds(unitOf) > 1e-9;
newCaps = base + scale .* folds(unitOf);
model.ub(usageIdx(relaxed)) = newCaps(relaxed);
model.ub(bioIdx) = max(model.ub(bioIdx), expGrowth);

flexEnz.uniprotIDs = model.ec.enzymes(measured(relaxed));
flexEnz.oldConcs   = base(relaxed);
flexEnz.flexConcs  = newCaps(relaxed);
poolFold = sol.x(poolCol);
if poolFold > 1e-9
    newPool = poolCap * (1 + poolFold);
    model.ub(poolIdx) = newPool;
    flexEnz.uniprotIDs(end+1, 1) = {'prot_pool'};
    flexEnz.oldConcs(end+1, 1)   = poolCap;
    flexEnz.flexConcs(end+1, 1)  = newPool;
    fprintf('%d enzymes relaxed, protein pool %g -> %g.\n', sum(relaxed), poolCap, newPool);
else
    fprintf('%d enzymes relaxed.\n', sum(relaxed));
end
flexEnz.frequence = ones(numel(flexEnz.uniprotIDs), 1);
flexEnz.ratioIncr = flexEnz.flexConcs ./ flexEnz.oldConcs;
[~, order] = sort(flexEnz.ratioIncr, 'descend');
flexEnz.uniprotIDs = flexEnz.uniprotIDs(order);
flexEnz.oldConcs   = flexEnz.oldConcs(order);
flexEnz.flexConcs  = flexEnz.flexConcs(order);
flexEnz.ratioIncr  = flexEnz.ratioIncr(order);
flexEnz.frequence  = flexEnz.frequence(order);
end
