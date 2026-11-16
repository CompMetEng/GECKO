function [model, changed] = zeroUnmeasuredIsozymes(model, varargin)
% zeroUnmeasuredIsozymes  Limit unmeasured isozymes that have a measured sibling.
%
% For a reaction with a measured and an unmeasured isozyme (see
% getIsozymeCoverage), the usage upper bound of the enzymes of the unmeasured
% isozyme is lowered, instead of leaving them to draw from the protein pool
% without a limit. The function lowers the usage upper bound in the model
% and does not change model.ec.concs. An enzyme that also catalyses a reaction
% without a measured sibling is left unchanged, because limiting it would
% limit that reaction as well.
%
% Not detecting a protein is not the same as the protein being absent: a
% protein below the detection limit of the measurement is missing from the
% data as well. Method 'zero' treats the unmeasured isozyme as absent and
% sets its usage upper bound to 0, which removes every flux through the
% enzyme and can lower the growth rate or make the model infeasible. Method
% 'detectionLimit' treats it as present at or below the detection limit and
% sets the upper bound to the detection limit.
%
% An enzyme lowered by this function has no concentration in model.ec.concs
% and is a candidate in relaxProteomicsGreedy like any other enzyme with an
% upper bound below the default.
%
% Run after constrainEnzConcs, which resets the usage upper bounds of all
% enzymes. The function is not called by other functions of GECKO.
%
% Parameters
% ----------
% model : struct
%     a full ecModel (with ecModel.ec structure), not
%     GECKO light, with concentrations in model.ec.concs, constrained by
%     constrainEnzConcs.
%
% Name-Value Arguments
% --------------------
% method : char
%     'zero' or 'detectionLimit' (default 'zero').
% limit : double
%     the detection limit in mg/gDCW for method 'detectionLimit' (default
%     the smallest concentration above 0 in model.ec.concs).
%
% Returns
% -------
% model : struct
%     an ecModel with lowered usage upper bounds. An upper bound that is
%     already below the new value is not raised.
% changed : table
%     one row per enzyme that the function acted on, with the variables
%     uniprot and newUpperBound.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     model = zeroUnmeasuredIsozymes(model);
%     model = zeroUnmeasuredIsozymes(model, 'method', 'detectionLimit');
%
% See also
% --------
% getIsozymeCoverage, constrainEnzConcs, fillEnzConcs

p = parseGECKOargs(varargin, { ...
    'method', 'zero'; ...
    'limit',  []});
method = p.method;
limit  = p.limit;
if isempty(method)
    method = 'zero';
end
if ~ismember(method, {'zero','detectionLimit'})
    error('zeroUnmeasuredIsozymes:method', ...
        'method must be ''zero'' or ''detectionLimit''.')
end

[~, zeroable] = getIsozymeCoverage(model);

if strcmp(method, 'zero')
    newUb = zeros(numel(zeroable), 1);
else
    if isempty(limit)
        positive = model.ec.concs(model.ec.concs > 0);
        if isempty(positive)
            limit = NaN;
        else
            limit = min(positive);
        end
    end
    if ~isscalar(limit) || ~isnumeric(limit) || isnan(limit) || limit < 0
        error('zeroUnmeasuredIsozymes:limit', ...
            'The detection limit must be a nonnegative number.')
    end
    newUb = repmat(double(limit), numel(zeroable), 1);
end

[found, usageIdx] = ismember(strcat('usage_prot_', zeroable), model.rxns);
if any(~found)
    error('zeroUnmeasuredIsozymes:noUsageReaction', ...
        'The model has no usage reaction for %s.', strjoin(zeroable(~found), ', '))
end
newUb = min(newUb(:), model.ub(usageIdx(:)));
usageIdx = usageIdx(:);
model.ub(usageIdx) = newUb;
changed = table(zeroable(:), newUb, 'VariableNames', {'uniprot','newUpperBound'});
end
