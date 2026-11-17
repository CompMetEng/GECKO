function model = getKcatAcrossIsozymes(model)
% getKcatAcrossIsozymes  Deprecated alias for fillKcatsFromIsozymes.
%
% Kept for backward compatibility with earlier GECKO 4 code. Will be removed
% in a future release; switch to fillKcatsFromIsozymes.
%
% See also
% --------
% fillKcatsFromIsozymes

warning('GECKO:deprecatedName', ...
    'getKcatAcrossIsozymes is deprecated; use fillKcatsFromIsozymes instead. The old name will be removed in a future release.');
model = fillKcatsFromIsozymes(model);
end
