function [enz, controlCoeffs] = getConcControlCoeffs(model, varargin)
% getConcControlCoeffs  Calculate control coefficients of protein usage.
%
% For each candidate protein whose usage reaction is running above the
% given limit fraction of its upper bound, temporarily increases that
% protein's usage upper bound by foldChange, re-solves, and expresses the
% resulting growth-rate increase per unit of concentration increase as a
% control coefficient. Proteins below the limit, or whose increase does
% not measurably raise growth, get a coefficient of zero.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% proteins : cell
%     a list of proteins to calculate the coefficients (default
%     model.ec.enzymes).
% foldChange : double
%     a value how much to increase the protein concentration (default 2).
% limit : double
%     threshold on the usage/concentration ratio (how saturated a
%     protein's usage upper bound already is); only proteins with a ratio
%     above this value are evaluated for a control coefficient (default
%     0, i.e. any nonzero usage is evaluated).
% method : char
%     'finite_difference' (default) or 'shadow_price'. See Notes.
%
% Returns
% -------
% enz : logical
%     flags which of the input proteins had a usage/concentration ratio
%     above limit and were evaluated for a control coefficient.
% controlCoeffs : double
%     control coefficient per input protein (0 for proteins not evaluated
%     or whose increase did not measurably raise growth).
%
% Notes
% -----
% The default method 'finite_difference' raises the upper bound of one usage
% reaction by foldChange, re-solves and divides the rise of the objective by
% the rise of the bound, one solve per evaluated protein. A rise of at most
% 1e-10 counts as no sensitivity.
%
% The method 'shadow_price' reads, from a single solve, the absolute
% reduced cost of each evaluated usage reaction: the local derivative of the
% objective to the upper bound of that reaction. A nonzero reduced cost is
% kept only when raising that bound by 0.01 percent raises the objective by
% at least half of the rise the reduced cost predicts, and by more than
% 1e-12 times max(1, |objective|) (one extra solve per nonzero value), as
% the reduced cost can overestimate the effect at a degenerate optimum. A
% usage reaction held at a positive lower bound gets 0. When the solver
% gives no reduced costs (as with SCIP), the 'finite_difference' method is
% used instead. The two methods agree where raising the bound does not
% change the set of active constraints.
%
% When the model has no feasible solution, the error
% getConcControlCoeffs:infeasible is raised.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [enz, controlCoeffs] = getConcControlCoeffs(model);
%     [enz, controlCoeffs] = getConcControlCoeffs(model, 'foldChange', 3);
%     [enz, controlCoeffs] = getConcControlCoeffs(model, 'method', 'shadow_price');
%
% See also
% --------
% flexibilizeEnzConcs

p = parseGECKOargs(varargin, { ...
    'proteins',   []; ...
    'foldChange', []; ...
    'limit',      []; ...
    'method',     'finite_difference'});
proteins   = p.proteins;
foldChange = p.foldChange;
limit      = p.limit;
method     = char(p.method);

if isempty(limit)
    limit = 0;
end

if isempty(foldChange)
    foldChange = 2;
end

if isempty(proteins)
    proteins = model.ec.enzymes;
end
if ~any(strcmp(method, {'finite_difference', 'shadow_price'}))
    error('getConcControlCoeffs:method', ...
        'method must be ''finite_difference'' or ''shadow_price'', not ''%s''.', method)
end

enz = false(length(proteins),1);
controlCoeffs = zeros(length(proteins),1);

[sol,hs] = solveLP(model);
if sol.stat ~= 1
    error('getConcControlCoeffs:infeasible', ...
        'The model has no optimal solution, so no control coefficients can be calculated.')
end
initialGrowth = sol.f;

% Get enzyme index
[~, protIdx] = ismember(proteins, model.ec.enzymes);

% Get the protein usage reactions
protUsageRxns = strcat('usage_prot_', model.ec.enzymes(protIdx));
[~, protUsageRxnIdx] = ismember(protUsageRxns, model.rxns);

if strcmp(method, 'shadow_price') && ~isempty(sol.rCost)
    for i = 1:numel(proteins)
        idx = protUsageRxnIdx(i);
        prevConc = model.ub(idx);
        if ~((sol.x(idx)/prevConc) > limit)
            continue
        end
        enz(i) = 1;
        rc = abs(sol.rCost(idx));
        heldAtLower = model.lb(idx) > 0 && sol.x(idx) <= model.lb(idx)*(1 + 1e-9);
        if rc == 0 || heldAtLower
            continue
        end
        % Confirm the reduced cost with a small raise of the bound (0.01 percent)
        step = 1e-4;
        tempModel = model;
        tempModel.ub(idx) = prevConc*(1 + step);
        [tempSol,hs] = solveLP(tempModel,0,[],hs);
        if tempSol.stat ~= 1
            continue
        end
        gain = tempSol.f - initialGrowth;
        if gain >= max(0.5*rc*prevConc*step, 1e-12*max(1, abs(initialGrowth)))
            controlCoeffs(i) = rc;
        end
    end
    return
end

for i = 1:numel(proteins)
    % Get the previous concentration
    prevConc = model.ub(protUsageRxnIdx(i));

    % Only consider those with a usage close the UB
    if (sol.x(protUsageRxnIdx(i))/prevConc) > limit
        
        % Update the logical vector
        enz(i) = 1;

        % Create a temporal model since coeff will be calculated one enzyme at
        % the time, without other change
        tempModel = model;
        % Increase the concentration by flexfactor
        newConc = prevConc*(foldChange);
        tempModel.ub(protUsageRxnIdx(i)) = newConc;

        % Get the new growth rate after the adjustment
        [tempSol,hs] = solveLP(tempModel,0,[],hs);
        tempGrowth = tempSol.f;
        
        % Calculate the coeff only if new growth rate is significantly
        % higher than initial value
        if (tempGrowth-initialGrowth)>1e-10
            controlCoeffs(i) = (tempGrowth-initialGrowth)/(newConc-prevConc);
        end
    end

end
end
