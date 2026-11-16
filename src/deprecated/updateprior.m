function varargout = updateprior(varargin) %#ok<STOUT>
% updateprior  Removed: helper of the Bayesian (ABC-SMC) kcat tuning of GECKO 3.
%
% This function raises an error that names its GECKO 4 replacement.
%
% See also
% --------
% cmaesKcatTuning, sensitivityTuning

error('GECKO:removedFunction', ['updateprior was removed in GECKO 4 with bayesianSensitivityTuning. ' ...
    'Use cmaesKcatTuning (advanced kcat tuning) or sensitivityTuning ' ...
    '(standard kcat tuning).']);
end
