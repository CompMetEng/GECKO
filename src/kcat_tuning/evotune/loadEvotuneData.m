function experimentData = loadEvotuneData(varargin)
% loadEvotuneData  Load the measured conditions of the model adapter's data folder.
%
% Reads the three data files that cmaesKcatTuning scores a model against
% from the data subfolder of the model directory defined by the model
% adapter, with RAVEN's loadExperimentData. The files list conditions in
% which growth rates and exchange fluxes were measured: the flux conditions
% of fluxData.tsv (the uptake of the carbon source is limited to the
% measured rate, and growth and the measured fluxes are compared) and the
% maximum-growth conditions of maxGrowth.tsv (the uptake is unlimited, and
% only the maximal growth rate is compared), and the reactions of
% zeroFlux.tsv that carry no flux. See loadExperimentData for the layout of
% the files. Each file is optional: a missing file yields [] (fluxData,
% maxGrowth) or an empty cell array (zeroFlux) rather than an error.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% experimentData : struct
%     structure with the loaded experimental data.
%
% Notes
% -----
% The experimentData structure has the following fields:
%
% - fluxData : the flux conditions of fluxData.tsv. Every condition in the
%   file is used by the tuning; a weight of 0 in the column weight leaves a
%   condition out of the simulation and the distance. [] if the file isn't
%   present.
% - maxGrowth : the maximum-growth conditions of maxGrowth.tsv. [] if the
%   file isn't present.
% - zeroFlux : reaction IDs assumed to carry zero flux, loaded from
%   zeroFlux.tsv. Empty cell array if the file isn't present.
%
% See also
% --------
% cmaesKcatTuning, evotuneScore, loadExperimentData
p = parseRAVENargs(varargin, { ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default ecModel adapter in the ModelAdapterManager.')
    end
end

experimentData = loadExperimentData(fullfile(modelAdapter.params.path,'data'));
end
