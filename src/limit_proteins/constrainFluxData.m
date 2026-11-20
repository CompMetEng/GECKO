function model = constrainFluxData(model, varargin)
% constrainFluxData  Constrain an ecModel to one measured flux condition.
%
% Use it to reproduce one measured condition before an analysis, for
% example a flux variability analysis or a proteomics integration. The
% biomass reaction is bounded by the measured growth rate and each exchange
% reaction of the data by its measured flux in the chosen condition. The
% data are the flux conditions of data/fluxData.tsv, read by RAVEN's
% loadFluxData. The preferred carbon source of the adapter is blocked first,
% so that the data alone decide which carbon source is used. The bounds are
% set by RAVEN's applyFluxDataConstraints, where the bounds of each setting
% of looseStrictFlux are listed; this function takes the biomass reaction
% and the carbon source from the model adapter.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% fluxData : struct
%     structure with flux data, with the fields:
%
%     - conds : sampling condition.
%     - Ptot : total protein (g/gDCW).
%     - grRate : growth rate (1/h).
%     - exchFluxes : exchange fluxes (mmol/gDCWh).
%     - exchMets : exchanged metabolites, matching exchFluxes.
%     - exchRxnIDs : exchange reaction IDs, matching exchMets.
%
%     Default: the data read from data/fluxData.tsv of the model adapter.
% condition : double or char
%     either index number or name of the sample condition in fluxData.conds
%     (default 1).
% maxMinGrowth : char
%     'max' if the provided growth rate should be set as maximum growth rate
%     (= upper bound), or 'min' if it should be set as minimum growth rate
%     (= lower bound). The latter option is suitable if minimization of
%     prot_pool_exchange is used as objective function (default 'max'). Any
%     other value raises constrainFluxData:invalidMaxMinGrowth.
% looseStrictFlux : char or double
%     how strictly constrained the exchange fluxes should be (default
%     'loose'):
%
%     - 'loose' : the exchange fluxes are constraint only by the "outer
%       bounds". If exchFluxes(i) > 0, LB = 0 and UB = exchFluxes(i). If
%       exchFluxes(i) < 0, LB = exchFluxes(i) and UB = 0.
%     - 0-100 : LB and UB constraints are set with a specified percentage of
%       variance around exchFluxes. If 10 is specified, LB = exchFluxes*0.95
%       and UB = exchFluxes*1.05. This allows for 10% variance around the
%       exchFluxes values, but strictly forces a flux through the exchRxns.
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% model : struct
%     an ecModel where fluxes are constraint.
%
% Notes
% -----
% An exchange reaction in fluxData.exchRxnIDs that is not in the model
% raises applyFluxDataConstraints:missingReaction, and a condition that is not
% in fluxData.conds raises applyFluxDataConstraints:unknownCondition. The
% preferred carbon source of the adapter (c_source) is blocked before the
% data is applied.
%
% If a provided constraint is either -1000 or 1000, then the function will
% update the reaction lower and upper bound to either allow uptake or
% excretion, irrespective of what option is given as the looseStrictFlux
% parameter.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     model = constrainFluxData(model);
%     model = constrainFluxData(model, 'condition', 2);
%
% See also
% --------
% loadFluxData, applyFluxDataConstraints

p = parseRAVENargs(varargin, { ...
    'fluxData',        []; ...
    'condition',       []; ...
    'maxMinGrowth',    []; ...
    'looseStrictFlux', []; ...
    'modelAdapter',    []});
fluxData        = p.fluxData;
modelAdapter    = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(fluxData)
    fluxData = loadFluxData(fullfile(params.path,'data','fluxData.tsv'));
end

options = {'cSource', params.c_source};
if ~isempty(p.condition)
    options = [options, {'condition', p.condition}];
end
if ~isempty(p.maxMinGrowth)
    options = [options, {'maxMinGrowth', p.maxMinGrowth}];
end
if ~isempty(p.looseStrictFlux)
    options = [options, {'looseStrictFlux', p.looseStrictFlux}];
end
model = applyFluxDataConstraints(model, fluxData, params.bioRxn, options{:});
end
