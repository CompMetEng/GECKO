function varargout = bayesianSensitivityTuning(varargin) %#ok<STOUT>
% bayesianSensitivityTuning  Removed: Bayesian (ABC-SMC) kcat tuning of GECKO 3.
%
% GECKO 4 tunes kcats against experimental data with cmaesKcatTuning (the
% advanced kcat tuning), which reads the same kinds of data from
% data/evotuneMaxGrowth.tsv, data/evotuneFluxData.tsv and
% data/evotuneZeroExch.tsv (see loadEvotuneData). sensitivityTuning is the
% standard kcat tuning. This function raises an error that names them.
%
% See also
% --------
% cmaesKcatTuning, screenKcatLeverage, sensitivityTuning

error('GECKO:removedFunction', ['bayesianSensitivityTuning was removed in GECKO 4. ' ...
    'Use cmaesKcatTuning (advanced kcat tuning, data in data/evotune*.tsv) or ' ...
    'sensitivityTuning (standard kcat tuning).']);
end
