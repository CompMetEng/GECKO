function varargout = writeDLKcatInput(varargin) %#ok<STOUT>
% writeDLKcatInput  Deprecated: use kcatsFromDLKcat(model, 'via', 'docker').
%
% kcatsFromDLKcat writes the DLKcat input of the ecModel, runs DLKcat in its
% Docker image and reads the predictions. Calling writeDLKcatInput raises
% an error. extractEnzymeSubstratePairs lists the enzyme-substrate pairs
% that the input holds.
%
% See also
% --------
% kcatsFromDLKcat, extractEnzymeSubstratePairs

warning('GECKO:deprecatedFunction', ['writeDLKcatInput is deprecated; use ' ...
    'kcatsFromDLKcat(model, ''via'', ''docker'').'])
error('GECKO:deprecatedFunction', ['writeDLKcatInput no longer writes the DLKcat input on ' ...
    'its own. Use kcatsFromDLKcat(ecModel, ''via'', ''docker''), which also runs DLKcat.'])
end
