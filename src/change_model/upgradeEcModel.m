function [model, changed] = upgradeEcModel(model)
% upgradeEcModel  Bring an ecModel made with GECKO 3 to the GECKO 4 format.
%
% GECKO 3 defined the protein usage reactions (usage_prot_<id>, prot_<id>
% -> prot_pool) and the protein pool exchange (prot_pool_exchange,
% prot_pool ->) in the reverse direction, used with negative flux. GECKO 4
% defines them forward (prot_pool -> prot_<id>, -> prot_pool), so that an
% enzyme's usage is the positive flux of its usage reaction. This turns
% each reaction that is still reverse: its column of S, its bounds and its
% objective coefficient change sign. Reactions already forward are left as
% they are, so calling it again changes nothing. Flux values of the model
% (fluxes of usage reactions and the pool exchange) change sign with it.
% geckopy has the same (upgrade_ec_model).
%
% Parameters
% ----------
% model : struct
%     an ecModel.
%
% Returns
% -------
% model : struct
%     the ecModel in the GECKO 4 format.
% changed : cell
%     the reactions that were turned.
%
% Examples
% --------
%     ecModel = upgradeEcModel(loadEcModel('ecYeastGEM_GECKO3.yml'));
%
% See also
% --------
% upgradeGECKOproject, makeEcModel, setProtPoolSize

pool = find(strcmp(model.mets, 'prot_pool'));
changed = {};
if isempty(pool)
    return
end
usage = find(startsWith(model.rxns, 'usage_prot_'));
reverse = usage(full(model.S(pool, usage)) > 0);
exchange = find(strcmp(model.rxns, 'prot_pool_exchange'));
if ~isempty(exchange) && full(model.S(pool, exchange)) < 0
    reverse(end+1, 1) = exchange;
end
if isempty(reverse)
    return
end
model.S(:, reverse) = -model.S(:, reverse);
lb = model.lb(reverse);
model.lb(reverse) = -model.ub(reverse);
model.ub(reverse) = -lb;
if isfield(model, 'c')
    model.c(reverse) = -model.c(reverse);
end
changed = model.rxns(reverse);
end
