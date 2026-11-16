function varargout = getcarbonnum(varargin) %#ok<STOUT>
% getcarbonnum  Removed: carbon numbers of exchange reactions, a helper of the Bayesian kcat tuning of GECKO 3.
%
% This function raises an error that names its GECKO 4 replacement.
%
% See also
% --------
% fillCarbonNum

error('GECKO:removedFunction', ['getcarbonnum was removed in GECKO 4. Use fillCarbonNum, which writes' ...
    ' the carbon number of each exchange reaction to model.excarbon.']);
end
