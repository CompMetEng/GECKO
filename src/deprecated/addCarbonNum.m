function varargout = addCarbonNum(varargin) %#ok<STOUT>
% addCarbonNum  Removed: carbon numbers of exchange reactions of GECKO 3.
%
% This function raises an error that names its replacement, RAVEN's
% computeCarbonWeights, which returns the carbon weight of each reaction as a
% containers.Map instead of writing the field model.excarbon.
%
% See also
% --------
% computeCarbonWeights

error('GECKO:removedFunction', ['addCarbonNum was removed in GECKO 4. Use computeCarbonWeights (RAVEN),' ...
    ' which returns the carbon weight of each reaction as a containers.Map.']);
end
