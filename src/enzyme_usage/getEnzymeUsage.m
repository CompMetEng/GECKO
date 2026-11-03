function usageData = getEnzymeUsage(ecModel,fluxes,varargin)
% getEnzymeUsage  Compute enzyme usages from a provided flux distribution.
%
% Gives enzyme usages based on a provided flux distribution, as obtained
% from a full GECKO model. It can give:
%
% 1. absolute usage: the specific enzyme usage in mg/gDCW, which can be
%    given for enzymes with- and without concentration information;
% 2. capacity usage: the ratio of available enzyme that is used, calculated
%    by (absUsage/UB) (note that capacity usage is 0 if an enzyme
%    concentration was not constrained in the model);
% 3. UB: the upper bound of each enzyme exchange reaction, which may not be
%    the same as the enzyme concentration, if it has been flexibilized;
% 4. protID: the protein identifiers for each enzyme, in the same order as
%    ecModel.ec.enzymes.
%
% Parameters
% ----------
% ecModel : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% fluxes : double
%     vector of fluxes, for instance sol.x.
%
% Name-Value Arguments
% --------------------
% zero : logical
%     whether also enzymes with zero absolute usage should be included
%     (default true).
%
% Returns
% -------
% usageData : struct
%     structure with enzyme usage data.
%
% Notes
% -----
% The usageData structure has the following fields:
%
% - capUsage : vector of enzyme capacity usages.
% - absUsage : vector of absolute enzyme usages.
% - UB : vector of enzyme exchange reaction upper bounds.
% - protID : string array of matching protein IDs. An enzyme without a
%   usage_prot_<id> reaction in the model is left out, with the warning
%   getEnzymeUsage:missingUsageReaction.
% - fluxes : vector of fluxes, copy of input fluxes.
% - measuredProtein : sum of the measured concentrations in model.ec.concs
%   (mg/gDCW), over the enzymes with a non-NaN entry.
% - demandedProtein : sum of absUsage over the same enzymes (mg/gDCW).
% - demandRatio : demandedProtein / measuredProtein. NaN when no
%   concentration is measured; Inf when the measured concentrations sum to 0
%   and the demand is positive. A ratio above 1 means that the solution uses
%   more of the measured enzymes than was measured, which happens when the
%   upper bounds of the usage reactions are above the measurements (for
%   instance after flexibilizeEnzConcs). A ratio at or below 1 does not
%   exclude that one enzyme exceeds its measurement, since the sum is over
%   all measured enzymes; capUsage gives the ratio per enzyme. Enzymes
%   without a measured concentration draw from the protein pool and are not
%   counted.
% - measuredProteinActive, demandRatioActive : the same two values restricted
%   to the enzymes with an absUsage above 1e-7 mg/gDCW. Measured enzymes
%   that are inactive in the solution raise measuredProtein and can hide a
%   conflict among the active ones.
%   The five fields are computed over all enzymes with a usage reaction,
%   whatever the value of zero.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     usageData = getEnzymeUsage(ecModel,fluxes,zero);
%     usageData = getEnzymeUsage(ecModel,fluxes,'zero',false);
%
% See also
% --------
% reportEnzymeUsage

p = parseGECKOargs(varargin, { ...
    'zero', true});
zero = p.zero;
if ecModel.ec.geckoLight
    error('This function does not work on GECKO light models.')
end
usageData.protID      = ecModel.ec.enzymes;
concs                 = nan(numel(ecModel.ec.enzymes),1);
nConcs                = min(numel(ecModel.ec.concs),numel(concs));
concs(1:nConcs)       = ecModel.ec.concs(1:nConcs);
[~,rxnIdx] = ismember(strcat('usage_prot_',ecModel.ec.enzymes),ecModel.rxns);
noUsageRxn = rxnIdx == 0;
if any(noUsageRxn)
    warning('getEnzymeUsage:missingUsageReaction', ...
        'The model has no usage reaction for the enzymes %s; they are left out.', ...
        strjoin(usageData.protID(noUsageRxn), ', '))
    usageData.protID = usageData.protID(~noUsageRxn);
    concs = concs(~noUsageRxn);
    rxnIdx = rxnIdx(~noUsageRxn);
end

usageData.UB          = ecModel.ub(rxnIdx);
usageData.absUsage    = abs(fluxes(rxnIdx));
usageData.capUsage    = abs(usageData.absUsage./usageData.UB);
usageData.fluxes      = fluxes;

measured = ~isnan(concs);
active   = usageData.absUsage > 1e-7;
usageData.measuredProtein       = sum(concs(measured));
usageData.measuredProteinActive = sum(concs(measured & active));
usageData.demandedProtein       = sum(usageData.absUsage(measured));
usageData.demandRatio           = usageData.demandedProtein / usageData.measuredProtein;
usageData.demandRatioActive     = sum(usageData.absUsage(measured & active)) / usageData.measuredProteinActive;

if ~zero
    nonzero               = usageData.absUsage>0;
    usageData.absUsage    = usageData.absUsage(nonzero);
    usageData.capUsage    = usageData.capUsage(nonzero);
    usageData.UB          = usageData.UB(nonzero);
    usageData.protID      = usageData.protID(nonzero);
end
end
