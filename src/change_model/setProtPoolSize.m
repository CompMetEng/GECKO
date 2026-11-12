function model = setProtPoolSize(model, varargin)
% setProtPoolSize  Set the limit of the total protein usage in the model.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% Ptot : double
%     total cellular protein content in g/gDCW. If not specified, the value
%     will be read from the model adapter. If not specified in model
%     adapter, 0.5 g/gDCW is assumed.
% f : double
%     estimated fraction of enzymes in the model. If not specified, the
%     value will be read from the model adapter. If not specified in model
%     adapter, 0.5 is assumed.
% sigma : double
%     estimated saturation factor. If not specified, the value will be read
%     from the model adapter. If not specified in model adapter, 0.5 is
%     assumed.
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% model : struct
%     ecModel with protein pool constraint set.
%
% Notes
% -----
% A model without a prot_pool_exchange reaction raises the error
% setProtPoolSize:noPoolExchange.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     model = setProtPoolSize(model, Ptot, f, sigma, modelAdapter);
%     model = setProtPoolSize(model, 'Ptot', Ptot);
%
% See also
% --------
% makeEcModel

p = parseGECKOargs(varargin, { ...
    'Ptot',         []; ...
    'f',            []; ...
    'sigma',        []; ...
    'modelAdapter', []});
Ptot         = p.Ptot;
f            = p.f;
sigma        = p.sigma;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter) && (isempty(Ptot) || isempty(f) || isempty(sigma))
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

if isempty(sigma)
    sigma = modelAdapter.getParameters().sigma;
end
if isempty(f)
    f = modelAdapter.getParameters().f;
end
if isempty(Ptot)
    Ptot = modelAdapter.getParameters().Ptot;
end

poolIdx = strcmp(model.rxns, 'prot_pool_exchange');
if ~any(poolIdx)
    error('setProtPoolSize:noPoolExchange', ...
        'The model has no prot_pool_exchange reaction; makeEcModel adds it.')
end
model.ub(poolIdx) = Ptot*f*sigma*1000;
end

