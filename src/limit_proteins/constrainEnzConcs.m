function model = constrainEnzConcs(model, varargin)
% constrainEnzConcs  Constrain enzyme usages by their concentration.
%
% Constrain enzyme usages by their concentration as provided in
% model.ec.concs. For enzymes with non-NaN entries in model.ec.concs, the
% upper bound of their enzyme usage reaction becomes the measured
% concentration. Their usage reactions keep drawing from the protein pool.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure) with enzyme
%     concentrations in model.ec.concs.
%
% Name-Value Arguments
% --------------------
% removeConstraints : logical
%     whether enzyme concentration constraints should be removed
%     (model.ec.concs will remain unchanged) (default false).
% restrictTo : cell
%     enzymes (Uniprot IDs) to process; the usage reactions of the other
%     enzymes are left unchanged. IDs that are not in model.ec.enzymes are
%     ignored (default: all enzymes in model.ec.enzymes).
%
% Returns
% -------
% model : struct
%     an ecModel constraint with available enzyme concentrations.
%
% Notes
% -----
% A negative concentration in model.ec.concs is contradictory data and
% raises the error constrainEnzConcs:negativeConcentration. To populate model.ec.concs you should run fillEnzConcs.
%
% A NaN entry means that the enzyme was not measured, and its usage draws
% from the protein pool. An entry of 0 means that the enzyme was measured as
% absent: its usage upper bound is 0, which removes every flux through the
% enzyme. flexibilizeEnzConcs multiplies the concentration and cannot raise a
% bound of 0; relaxProteomicsGreedy sets the bound back to its default.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     model = constrainEnzConcs(model);
%     model = constrainEnzConcs(model, 'removeConstraints', true);
%     model = constrainEnzConcs(model, 'restrictTo', {'P00330'});
%
% See also
% --------
% fillEnzConcs

%Enzyme with NaN entry in model.ec.concs => draw from prot_pool
%Enzyme with numeric entry in model.ec.concs => exchange reaction with
%enzyme level as UB

p = parseGECKOargs(varargin, { ...
    'removeConstraints', false; ...
    'restrictTo',        []});
removeConstraints = p.removeConstraints;
if isempty(p.restrictTo) && isnumeric(p.restrictTo)
    selected = true(numel(model.ec.enzymes), 1);
else
    selected = ismember(model.ec.enzymes, cellstr(p.restrictTo));
end

%Get indices of usage reactions 
usageRxns = strcat('usage_prot_',model.ec.enzymes);
[~, usageRxnsIdx] = ismember(usageRxns, model.rxns);

if any(usageRxnsIdx == 0)
    error('Usage reactions are not defined for all enzymes. This is done by makeEcModel.')
end
%Get index of protein pool metabolite
protPoolIdx = find(ismember(model.mets,'prot_pool'));
if ~any(protPoolIdx)
    error('Cannot find protein pool pseudometabolite.')
end

%Protein that should be constraint by UB
if removeConstraints
    protCons = false(numel(model.ec.enzymes), 1);
else
    protCons = selected & ~isnan(model.ec.concs);
    if any(model.ec.concs(protCons) < 0)
        error('constrainEnzConcs:negativeConcentration', ...
            'model.ec.concs has a negative concentration for %s.', ...
            strjoin(model.ec.enzymes(protCons & model.ec.concs < 0), ', '))
    end
end

%Set the selected reactions to draw from prot_pool (forward: prot_pool -> prot_<id>)
model.S(protPoolIdx, usageRxnsIdx(selected)) = -1;
model.lb(usageRxnsIdx(selected)) = 0;
model.ub(usageRxnsIdx(selected)) = 1000;

%If non-NaN in model.ec.concs, then constrain by UB
if any(protCons)
%    model.S(protPoolIdx, usageRxnsIdx(protCons)) = 0; % Since GECKO 3.2.0
    model.ub(usageRxnsIdx(protCons)) = model.ec.concs(protCons);
end
end
