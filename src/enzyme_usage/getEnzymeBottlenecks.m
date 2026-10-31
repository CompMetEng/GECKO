function bottlenecks = getEnzymeBottlenecks(model, varargin)
% getEnzymeBottlenecks  Top-N enzymes limiting the current objective.
%
% Solves the model and ranks every enzyme by the change of the objective
% after the upper bound of its usage_prot_<id> reaction is raised by
% foldChange and the model is solved again. The objective is maximised, so
% the gain is the rise of the maximised objective. Only enzymes whose usage runs at
% its upper bound are solved again; raising a bound that does not bind
% cannot change the optimum, so every other enzyme has an objectiveGain of 0
% without a solve.
%
% The shadow price of the prot_<id> mass-balance constraint is reported for
% reference and does not enter the order: at a degenerate optimum, which is
% common in enzyme-constrained models, several constraints are active and
% the solver assigns the dual to one of them depending on the basis.
% Enzymes with equal objectiveGain (rounded to 9 decimals) keep their order
% in model.ec.enzymes.
%
% Parameters
% ----------
% model : struct
%     a full ecModel in GECKO 3 format (with ecModel.ec structure). Not
%     supported for gecko-light models, which have no per-enzyme
%     prot_<id> / usage_prot_<id> machinery.
%
% Name-Value Arguments
% --------------------
% top : double
%     number of enzymes to return (default 10).
% foldChange : double
%     a value how much to increase the usage upper bound of an enzyme that
%     runs at its upper bound; must be greater than 1 (default 2).
%
% Returns
% -------
% bottlenecks : table
%     one row per returned enzyme, sorted by descending objectiveGain, with
%     columns uniprot, gene, objectiveGain (objective after raising the
%     bound minus the objective of the unchanged model; 0 for enzymes not at
%     their upper bound, for gains below 1e-10 and for solves that are not
%     optimal), shadowPrice, flux, capUsage (abs(flux)/upperBound, NaN if
%     upperBound is 0), upperBound and multipleActive. An enzyme is at its
%     upper bound when upperBound is above 0 and capUsage is at least
%     1 - 1e-6. multipleActive is true when at least one other constraint
%     is active: another enzyme at its upper bound, or prot_pool_exchange at
%     its upper bound (within 1e-6 relative). Where it is true, the dual is
%     not informative about a single enzyme. The table property
%     UserData.nResolved holds the number of enzymes that were solved again
%     and UserData.activeConstraints the number of active constraints.
%
% Examples
% --------
%     bottlenecks = getEnzymeBottlenecks(ecModel);
%     bottlenecks = getEnzymeBottlenecks(ecModel, 'top', 5, 'foldChange', 3);

p = parseGECKOargs(varargin, {'top', []; 'foldChange', []});
top = p.top;
foldChange = p.foldChange;
if isempty(top)
    top = 10;
end
if isempty(foldChange)
    foldChange = 2;
end
if ~isnumeric(foldChange) || ~isscalar(foldChange) || ~(foldChange > 1)
    error('getEnzymeBottlenecks:foldChange', 'foldChange must be greater than 1')
end

if isfield(model.ec, 'geckoLight') && model.ec.geckoLight
    error(['getEnzymeBottlenecks requires a full ecModel, not gecko-light: ' ...
        'light models have no per-enzyme prot_<id> / usage_prot_<id> machinery.'])
end

[sol, hs] = solveLP(model);
if sol.stat ~= 1
    error('getEnzymeBottlenecks: solver status %d; cannot rank bottlenecks.', sol.stat)
end

uniprot = model.ec.enzymes(:);
gene    = model.ec.genes(:);

protMets  = strcat('prot_', uniprot);
[found, metIdx] = ismember(protMets, model.mets);
if ~all(found)
    error('Cannot find prot_<id> metabolites for all enzymes. This is done by makeEcModel.')
end
usageRxns = strcat('usage_prot_', uniprot);
[found, rxnIdx] = ismember(usageRxns, model.rxns);
if ~all(found)
    error('Usage reactions are not defined for all enzymes. This is done by makeEcModel.')
end

shadowPrice = sol.sPrice(metIdx);
flux        = sol.x(rxnIdx);
upperBound  = model.ub(rxnIdx);
capUsage    = abs(flux) ./ upperBound;
capUsage(upperBound == 0) = NaN;

nEnz = numel(uniprot);
atCap = upperBound > 0 & capUsage >= 1 - 1e-6;
base = sol.f;
objectiveGain = zeros(nEnz, 1);
for i = find(atCap)'
    tempModel = model;
    tempModel.ub(rxnIdx(i)) = upperBound(i) * foldChange;
    [tempSol, hs] = solveLP(tempModel, 0, [], hs);
    if tempSol.stat == 1
        gain = tempSol.f - base;
        if gain >= 1e-10
            objectiveGain(i) = gain;
        end
    end
end

activeConstraints = nnz(atCap);
poolIdx = find(strcmp(model.rxns, 'prot_pool_exchange'), 1);
if ~isempty(poolIdx) && isfinite(model.ub(poolIdx)) && model.ub(poolIdx) > 0 ...
        && abs(sol.x(poolIdx)) >= model.ub(poolIdx) * (1 - 1e-6)
    activeConstraints = activeConstraints + 1;
end
multipleActive = (activeConstraints - double(atCap)) >= 1;

bottlenecks = table(uniprot, gene, objectiveGain, shadowPrice, flux, capUsage, upperBound, multipleActive);
[~, order] = sortrows([-round(objectiveGain, 9), (1:nEnz)']);
bottlenecks = bottlenecks(order, :);
bottlenecks = bottlenecks(1:min(top, height(bottlenecks)), :);
bottlenecks.Properties.UserData.nResolved = nnz(atCap);
bottlenecks.Properties.UserData.activeConstraints = activeConstraints;
end
