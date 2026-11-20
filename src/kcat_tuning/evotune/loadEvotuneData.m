function experimentData = loadEvotuneData(varargin)
% loadEvotuneData  Load experimental data for evotune (CMA-ES) kcat tuning.
%
% Reads the three experimental data files used by cmaesKcatTuning from the
% data subfolder of the model directory defined by the model adapter, with
% RAVEN's loadExperimentData. Each file is optional: a missing file yields []
% (fluxData, maxGrowth) or an empty cell array (zeroFlux) rather than an
% error.
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
% - fluxData : flux data loaded from fluxData.tsv. Every condition in the
%   file is used by the tuning; a weight of 0 in the column weight leaves a
%   condition out. [] if the file isn't present.
% - maxGrowth : maximum growth-rate data loaded from maxGrowth.tsv. [] if
%   the file isn't present.
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
