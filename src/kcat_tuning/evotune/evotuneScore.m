function [objective,rmse,rmseDetail,carbonWeights] = evotuneScore(ecModel,experimentData,varargin)
% evotuneScore  Score one ecModel's current kcats against the experimental data.
%
% Internal helper used by cmaesKcatTuning, screenKcatLeverage and
% tunePriorPenaltyWeight to score one ecModel (with a given set of kcats
% already applied via applyKcatConstraints) against the measured growth
% rates and exchange fluxes. The flux conditions and the maximum-growth
% conditions are simulated with RAVEN's simulateConditions and combined into
% one distance with experimentDistance (see there for the definitions);
% this function supplies the model-specific parts: the biomass reaction,
% the anaerobic and protein-content functions of the model adapter, and the
% penalty on the kcats.
%
% Parameters
% ----------
% ecModel : struct
%     ecModel with the current candidate kcats already applied (via
%     applyKcatConstraints). Not mutated.
% experimentData : struct
%     structure with experimental data to be used, as returned by
%     loadEvotuneData.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% maxGrowthWeight : double
%     weight of the maximum-growth conditions against the flux conditions in
%     the distance: (rmseFlux + w*rmseMaxGrowth) / (w+1) (default 1.0, which
%     averages the two; 2 counts the maximum-growth conditions twice as
%     much).
% priorPenaltyWeight : double
%     weight on a prior term added to the objective (default 0, i.e. the
%     objective equals the distance). See kcatVec/kcat0/sigma0log below.
% kcatVec, kcat0, sigma0log : double
%     only used when priorPenaltyWeight > 0: the tunable kcat subset's
%     current values, prior values, and per-kcat prior log-space std dev,
%     all the same length and in the same order. The penalty added to the
%     objective is priorPenaltyWeight * mean((log(kcatVec./kcat0)./sigma0log).^2).
% carbonWeights : containers.Map
%     carbon weights of the reactions, as computed by computeCarbonWeights
%     (default: computed here).
%
% Returns
% -------
% objective : double
%     rmse, plus the prior penalty term when priorPenaltyWeight > 0. This
%     is what a search should minimise.
% rmse : double
%     the distance between the model and the measurements (carbon-weighted
%     RMSE in Cmmol/gDW/h, without the prior penalty), always comparable
%     across runs regardless of priorPenaltyWeight.
% rmseDetail : cell
%     per-condition RMSE values, labelled fluxData_<i> for the i-th flux
%     condition and maxGrowth_<i> for the i-th maximum-growth condition. A
%     condition with weight 0 is left out of the simulation and has NaN.
% carbonWeights : containers.Map
%     the carbon weights that were used. The weights depend on the
%     metabolite formulas and not on the kcats, so callers scoring many
%     candidates in a loop should pass this output back in as carbonWeights
%     for the next candidate, rather than recomputing it (an expensive
%     metabolite-formula scan) on every call.
%
% See also
% --------
% cmaesKcatTuning, screenKcatLeverage, loadEvotuneData, experimentDistance

p = parseRAVENargs(varargin, { ...
    'modelAdapter', []; ...
    'maxGrowthWeight', 1.0; ...
    'priorPenaltyWeight', 0.0; ...
    'kcatVec', []; ...
    'kcat0', []; ...
    'sigma0log', []; ...
    'carbonWeights', []});
modelAdapter        = p.modelAdapter;
maxGrowthWeight     = p.maxGrowthWeight;
priorPenaltyWeight  = p.priorPenaltyWeight;
carbonWeights       = p.carbonWeights;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end
bioRxnId = modelAdapter.params.bioRxn;

if isempty(carbonWeights)
    rxnIds = experimentData.zeroFlux(:);
    if ~isempty(experimentData.fluxData)
        rxnIds = [rxnIds; experimentData.fluxData.exchRxnIDs(:)];
    end
    if ~isempty(experimentData.maxGrowth)
        rxnIds = [rxnIds; experimentData.maxGrowth.exchRxnIDs(:)];
    end
    carbonWeights = computeCarbonWeights(ecModel, unique([rxnIds; {bioRxnId}]), bioRxnId);
end

makeAnaerobic = @(model) modelAdapter.makeModelAnaerobic(model);
changeProteinBiomass = @(model,ptot) modelAdapter.changeProteinBiomass(model,ptot);

fluxSims = [];
maxGrowthSims = [];
if ~isempty(experimentData.fluxData)
    fluxSims = simulateConditions(ecModel, experimentData.fluxData, true, ...
        experimentData.zeroFlux, bioRxnId, 'makeAnaerobic', makeAnaerobic, ...
        'changeProteinBiomass', changeProteinBiomass);
end
if ~isempty(experimentData.maxGrowth)
    maxGrowthSims = simulateConditions(ecModel, experimentData.maxGrowth, false, ...
        {}, bioRxnId, 'makeAnaerobic', makeAnaerobic, ...
        'changeProteinBiomass', changeProteinBiomass);
end
[rmse, detail] = experimentDistance(experimentData, fluxSims, maxGrowthSims, ...
    carbonWeights, bioRxnId, 'maxGrowthWeight', maxGrowthWeight);

rmseDetail = [];
if isfield(detail,'fluxData')
    rmseDetail = [rmseDetail; [cellstr("fluxData_" + string(1:numel(detail.fluxData)))', num2cell(detail.fluxData(:))]];
end
if isfield(detail,'maxGrowth')
    rmseDetail = [rmseDetail; [cellstr("maxGrowth_" + string(1:numel(detail.maxGrowth)))', num2cell(detail.maxGrowth(:))]];
end

objective = rmse;
if priorPenaltyWeight > 0 && ~isempty(p.kcatVec) && ~isempty(p.kcat0) && ~isempty(p.sigma0log)
    dev = log(p.kcatVec(:) ./ p.kcat0(:)) ./ p.sigma0log(:);
    objective = rmse + priorPenaltyWeight * mean(dev.^2);
end
end
