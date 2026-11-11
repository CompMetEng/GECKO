function runDLKcat(varargin)
% runDLKcat  Deprecated: use kcatsFromDLKcat(model, 'via', 'docker').
%
% kcatsFromDLKcat writes the DLKcat input of the ecModel, runs DLKcat in its
% Docker image and reads the predictions. Calling runDLKcat raises an
% error.
%
% See also
% --------
% kcatsFromDLKcat

warning('GECKO:deprecatedFunction', ['runDLKcat is deprecated; use ' ...
    'kcatsFromDLKcat(model, ''via'', ''docker'').'])
error('GECKO:deprecatedFunction', ['runDLKcat no longer runs DLKcat on its own. Use ' ...
    'kcatsFromDLKcat(ecModel, ''via'', ''docker''), which also writes the input.'])
end
