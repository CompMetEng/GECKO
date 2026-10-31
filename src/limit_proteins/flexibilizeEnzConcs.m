function [model, flexEnz] = flexibilizeEnzConcs(model, varargin)
% flexibilizeEnzConcs  Flexibilize enzyme concentrations of an ecModel.
%
% Flexibilize enzyme concentration of an ecModel constrained with proteomics
% data. The upper bound of the protein usage reaction is changed, while the
% concentrations in ecModel.ec.concs remain unchanged.
%
% If no (more) limiting enzyme concentrations can be found, it might be the
% protein pool exchange that is limiting growth. In that case, an attempt
% will be made to relax the protein pool exchange reaction, and if the growth
% rate indeed increases, it is suggested to set the protein pool exchange
% unconstrained (ub=1000) before again running flexibilizeEnzConcs. Such
% situations, where a proteomics integrated ecModel is overconstrained, may
% occur if the ecModel should be able to simulate maximum growth rate (from
% e.g. batch cultivation).
%
% If relaxing the protein pool exchange does not increase the growth rate,
% then this is not due to enzyme constraints, but rather an issue with the
% metabolic network itself, or the set nutrient exchange is not sufficient.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% expGrowth : double
%     estimated experimental growth rate. If not specified, the value will be
%     read from the model adapter.
% foldChange : double
%     a value how much to increase the enzyme concentration (default 2).
% iterPerEnzyme : double
%     the number of iterations that an enzyme can be increased. A zero number
%     can be defined; if zero is defined no limit will be set, and it will
%     increase the enzyme concentration until it reaches the defined growth
%     rate (default 5).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% verbose : logical
%     whether progress should be reported (default true).
% perComplex : logical
%     whether every complex (see getComplexFamilies) is flexibilized as one
%     enzyme (default false). Enzymes and complexes are then ranked by the
%     growth gained per foldChange of raising all their usage upper bounds
%     together (a levelled complex binds on all subunits at once, so no
%     single subunit shows a gain), every measured subunit of the picked
%     complex is raised by the same factor, iterations are counted per
%     complex, and one factor per complex is kept in the final pass, so that
%     concentrations levelled by levelComplexConcs stay levelled. An upper
%     bound is never lowered below its current value. The frequence of a
%     complex is recorded on its first measured subunit.
%
% Returns
% -------
% model : struct
%     ecModel where the constraint of measured enzyme concentrations have
%     been flexibilized (=relaxed) in the S-matrix, to allow the model to
%     reach the growth rate defined in expGrowth, while the values in
%     ecModel.ec.concs have remained untouched.
% flexEnz : struct
%     array with information about flexibilized proteins, with the fields:
%
%     - uniprotIDs : enzymes whose usage UB was flexibilized.
%     - oldConcs : original concentrations, from model.ec.concs.
%     - flexConcs : flexibilized concentrations, new UB in model.
%     - ratioIncr : ratio by which the concentration increased; the enzymes
%       will be sorted by this field.
%     - frequence : numeric how often the enzyme has been step-wise
%       flexibilized.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [model, flexEnz] = flexibilizeEnzConcs(model);
%     [model, flexEnz] = flexibilizeEnzConcs(model, 'foldChange', 3);
%
% See also
% --------
% constrainEnzConcs, getConcControlCoeffs, levelComplexConcs, getComplexFamilies

p = parseGECKOargs(varargin, { ...
    'expGrowth',     []; ...
    'foldChange',    []; ...
    'iterPerEnzyme', []; ...
    'modelAdapter',  []; ...
    'verbose',       []; ...
    'perComplex',    []});
expGrowth     = p.expGrowth;
foldChange    = p.foldChange;
iterPerEnzyme = p.iterPerEnzyme;
modelAdapter  = p.modelAdapter;
verbose       = p.verbose;
perComplex    = p.perComplex;

if isempty(verbose)
    verbose = true;
end
if isempty(perComplex)
    perComplex = false;
end
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(iterPerEnzyme)
    iterPerEnzyme = 5;
end

if isempty(foldChange)
    foldChange = 2;
end

if isempty(expGrowth)
    expGrowth = modelAdapter.getParameters().gR_exp;
end

% If a zero value is defined, not iteration limit will be set
if iterPerEnzyme == 0
    iterPerEnzyme = inf;
end

bioRxnIdx = getIndexes(model,params.bioRxn,'rxns');
if model.ub(bioRxnIdx) < expGrowth
    fprintf(['The upper bound of the biomass reaction was %s, while expGrowth '...
             'is %s. The upper bound has been set to expGrowth.'],...
             num2str(model.ub(bioRxnIdx)), num2str(expGrowth))
    model.ub(bioRxnIdx) = expGrowth;
end

% Keep track if protein pool will be flexibilized
changedProtPool = false;

% In case the model have not been protein constrained
% model = constrainEnzConcs(model);

[sol,hs] = solveLP(model);
predGrowth = sol.f;

% Get those proteins with a concentration defined
protConcs = find(~isnan(model.ec.concs));

% Get enzymes names with a concentration value
proteins = model.ec.enzymes(protConcs);

% Store how many time is predicted a enzyme with the highest coeff
frequence = zeros(numel(proteins),1);

if perComplex
    [model, flexEnz] = flexibilizePerComplex(model, proteins, protConcs, ...
        predGrowth, expGrowth, foldChange, iterPerEnzyme, params, verbose);
    return
end

flexBreak=false;
if any(protConcs)
    while predGrowth < expGrowth

        % Get the control coefficients
        [~, controlCoeffs] = getConcControlCoeffs(model, proteins, foldChange, 0.75);

        % Stop looking for limiting proteins all coeffs are zero
        if any(controlCoeffs)

            % Find the maximum coefficient index
            [~,maxIdx] = max(controlCoeffs);

            % Increase +1 the frequence
            frequence(maxIdx) = frequence(maxIdx) + 1;

            % Allow to increase the UB maximum five times
            if frequence(maxIdx) <= iterPerEnzyme

                % Set how much will increase the UB
                increase = foldChange*frequence(maxIdx);

                % Get usage rxn for the highest coeff
                protUsageIdx = strcmpi(model.rxns, strcat('usage_prot_', proteins(maxIdx)));

                % replace the UB
                model.ub(protUsageIdx) = model.ec.concs(protConcs(maxIdx)) * (1+increase);

                % Get the new growth rate
                sol = solveLP(model,0,[],hs);
                predGrowth = sol.f;

                if verbose
                    disp(['Protein ' proteins{maxIdx} ' UB adjusted. Grow: ' num2str(predGrowth)])
                end
            else
                printOrange(['WARNING: Limit has been reached. Protein '  proteins{maxIdx} ' seems to be problematic. Consider changing the kcat.\n']);
                flexBreak=true;
                break
            end
        else
            disp('No (more) limiting enzymes have been found. Attempting to increase protein pool exchange...')
            protPoolIdx = getIndexes(model,'prot_pool_exchange','rxns');
            oldProtPool = model.ub(protPoolIdx);
            tempModel = setParam(model,'ub',protPoolIdx,1000);
            tempModel = setParam(tempModel,'ub',bioRxnIdx,expGrowth);
            sol = solveLP(tempModel);
            if (sol.f-predGrowth)>1e-3 % There is improvement in growth rate
                % Find new protein pool constraint
                predGrowth = sol.f;
                tempModel = setParam(tempModel,'lb',bioRxnIdx,predGrowth);
                % minimise pool usage (solveLP maximises the objective, so use -1)
                tempModel = setParam(tempModel,'obj',protPoolIdx,-1);
                sol = solveLP(tempModel);
                newProtPool = sol.x(protPoolIdx);
                model.ub(protPoolIdx) = newProtPool;
                fprintf(['Changing the upper bound of protein pool exchange from %s ',...
                         'to %s enabled a\ngrowth rate of %s. It can be helpful to set ', ...
                         'the upper bound of protein\npool exchange to 1000 before ',...
                         'running flexibilizeEnzConcs.\n'],num2str(oldProtPool), ...
                         num2str(newProtPool),num2str(predGrowth));
                changedProtPool = true;
            else
                fprintf(['Protein pool exchange was also not limiting. Inability to reach growth ',...
                        'rate is not related to\nenzyme constraints. Maximum growth rate ',...
                        'is %s.\n'], num2str(predGrowth))
            end
            % To return an usable model with these conditions, set the
            % highest growth rate reached as the experimental value
            expGrowth = predGrowth;
            break
        end
    end
else
    error('Protein concentrations have not been defined. Please run readProteomics and constrainEnzConcs')
end

protFlex         = proteins(frequence>0);
[~,protUsageIdx] = ismember(strcat('usage_prot_', protFlex), model.rxns);
[~,ecProtId]     = ismember(protFlex, model.ec.enzymes);
oldConcs         = model.ec.concs(ecProtId);
%frequence        = frequence(ecProtId);

if flexBreak == false && ~isempty(protFlex)
    % Not all flexibilized proteins require flexibilization in the end
    % Test flux distribution with minimum prot_pool usage
    modelTemp       = setParam(model,'var',params.bioRxn,expGrowth,0.5);
    % minimise pool usage (solveLP maximises the objective, so use -1)
    modelTemp       = setParam(modelTemp,'obj','prot_pool_exchange',-1);
    sol             = solveLP(modelTemp,1);
    % Check which concentrations can remain unchanged
    newConcs        = abs(sol.x(protUsageIdx));
    keepOld         = newConcs<oldConcs;
    % Set UBs
    model.ub(protUsageIdx(keepOld))  = oldConcs(keepOld);
    model.ub(protUsageIdx(~keepOld)) = newConcs(~keepOld);
    % Output structure
    protFlex(keepOld) = [];
    oldConcs(keepOld) = [];
    newConcs(keepOld) = [];
    frequence=frequence(frequence>0);
    frequence(keepOld) = [];
    
    flexEnz.uniprotIDs = protFlex;
    flexEnz.oldConcs   = oldConcs;
    flexEnz.flexConcs  = newConcs;
    flexEnz.frequence  = frequence(frequence>0);
else
    protFlex           = proteins(frequence>0);
    [~, protUsageIdx]  = ismember(strcat('usage_prot_', protFlex), model.rxns);

    flexEnz.uniprotIDs = protFlex;
    flexEnz.oldConcs   = model.ec.concs(ecProtId);
    flexEnz.flexConcs  = model.ub(protUsageIdx);
    flexEnz.frequence  = frequence(frequence>0);
end
if changedProtPool
    flexEnz.uniprotIDs(end+1) = {'prot_pool'};
    flexEnz.oldConcs(end+1)   = oldProtPool;
    flexEnz.flexConcs(end+1)  = newProtPool;
    flexEnz.frequence(end+1)  = 1;
end
% Sort by biggest ratio increase
flexEnz.ratioIncr  = flexEnz.flexConcs./flexEnz.oldConcs;
[~,b] = sort(flexEnz.ratioIncr,'descend');
flexEnz.ratioIncr  = flexEnz.ratioIncr(b);
flexEnz.uniprotIDs = flexEnz.uniprotIDs(b);
flexEnz.oldConcs   = flexEnz.oldConcs(b);
flexEnz.flexConcs  = flexEnz.flexConcs(b);
flexEnz.frequence  = flexEnz.frequence(b);
end

function [model, flexEnz] = flexibilizePerComplex(model, proteins, protConcs, predGrowth, expGrowth, foldChange, iterPerEnzyme, params, verbose)
% Flexibilize with every complex as one unit, see the perComplex option.
if isempty(protConcs)
    error('Protein concentrations have not been defined. Please run readProteomics and constrainEnzConcs')
end
limit = 0.75;
nProt = numel(proteins);
[~, usageIdx] = ismember(strcat('usage_prot_', proteins), model.rxns);
bioRxnIdx = getIndexes(model, params.bioRxn, 'rxns');
families = getComplexFamilies(model);

% Positions in proteins of each unit that is relaxed together: the measured
% subunits that belong to a complex, or a single enzyme
units    = {};
assigned = false(nProt, 1);
for i = 1:nProt
    if assigned(i)
        continue
    end
    [~, e] = ismember(proteins{i}, model.ec.enzymes);
    c = families.owner(e);
    members = i;
    if c > 0
        cplx = families.complexes{c};
        [~, members] = ismember(cplx.enzymes(cplx.owned), proteins);
        members = members(members > 0)';
    end
    assigned(members) = true;
    units{end+1, 1} = members; %#ok<AGROW>
end
unitFrequence = zeros(numel(units), 1);
frequence     = zeros(nProt, 1);
touched       = false(nProt, 1);

% Backstop on the iterations when the number per unit is not limited
if isinf(iterPerEnzyme)
    capPerEnzyme = 1000;
else
    capPerEnzyme = iterPerEnzyme;
end
maxTotalIter = nProt * capPerEnzyme + nProt;
iterCount    = 0;
flexBreak    = false;
changedProtPool = false;

while predGrowth < expGrowth
    iterCount = iterCount + 1;
    if iterCount > maxTotalIter
        printOrange(['WARNING: The global iteration limit has been reached without reaching expGrowth. ' ...
                     'The target may be unreachable even with relaxed enzyme bounds.\n']);
        flexBreak = true;
        break
    end
    scores = getUnitScores(model, usageIdx, units, foldChange, limit);

    if ~any(scores > 0)
        disp('No (more) limiting enzymes have been found. Attempting to increase protein pool exchange...')
        protPoolIdx = getIndexes(model,'prot_pool_exchange','rxns');
        oldProtPool = model.ub(protPoolIdx);
        tempModel = setParam(model,'ub',protPoolIdx,1000);
        tempModel = setParam(tempModel,'ub',bioRxnIdx,expGrowth);
        sol = solveLP(tempModel);
        if (sol.f-predGrowth)>1e-3
            predGrowth = sol.f;
            tempModel = setParam(tempModel,'lb',bioRxnIdx,predGrowth);
            tempModel = setParam(tempModel,'obj',protPoolIdx,-1);
            sol = solveLP(tempModel);
            newProtPool = sol.x(protPoolIdx);
            model.ub(protPoolIdx) = newProtPool;
            fprintf(['Changing the upper bound of protein pool exchange from %s ',...
                     'to %s enabled a\ngrowth rate of %s. It can be helpful to set ', ...
                     'the upper bound of protein\npool exchange to 1000 before ',...
                     'running flexibilizeEnzConcs.\n'],num2str(oldProtPool), ...
                     num2str(newProtPool),num2str(predGrowth));
            changedProtPool = true;
        else
            fprintf(['Protein pool exchange was also not limiting. Inability to reach growth ',...
                    'rate is not related to\nenzyme constraints. Maximum growth rate ',...
                    'is %s.\n'], num2str(predGrowth))
        end
        expGrowth = predGrowth;
        break
    end

    [~, unit] = max(scores);
    members = units{unit};
    frequence(members(1)) = frequence(members(1)) + 1;
    unitFrequence(unit)   = unitFrequence(unit) + 1;
    if unitFrequence(unit) > iterPerEnzyme
        if numel(members) > 1
            what = ['complex of ' strjoin(proteins(members)', ', ')];
        else
            what = ['Protein ' proteins{members(1)}];
        end
        printOrange(['WARNING: Limit has been reached. ' what ' seems to be problematic. Consider changing the kcat.\n']);
        flexBreak = true;
        break
    end

    increase = foldChange * unitFrequence(unit);
    for i = members
        touched(i) = true;
        model.ub(usageIdx(i)) = max(model.ub(usageIdx(i)), model.ec.concs(protConcs(i)) * (1+increase));
    end
    sol = solveLP(model);
    if sol.stat == 1
        predGrowth = sol.f;
    else
        predGrowth = 0;
    end
    if verbose
        disp(['Protein ' proteins{members(1)} ' UB adjusted. Grow: ' num2str(predGrowth)])
    end
end

flexIdx  = find(frequence > 0 | touched);
protFlex = proteins(flexIdx);
oldConcs = model.ec.concs(protConcs(flexIdx));
freq     = frequence(flexIdx);
if isempty(flexIdx)
    flexConcs = zeros(0, 1);
elseif flexBreak
    flexConcs = model.ub(usageIdx(flexIdx));
else
    % Not all flexibilized proteins require flexibilization in the end.
    % Test flux distribution with minimum prot_pool usage
    modelTemp = setParam(model,'var',params.bioRxn,expGrowth,0.5);
    modelTemp = setParam(modelTemp,'obj','prot_pool_exchange',-1);
    sol       = solveLP(modelTemp,1);
    newConcs  = abs(sol.x(usageIdx(flexIdx)));
    % One factor per complex: every flexibilized subunit takes the largest
    % new/old of its complex
    seen = false(numel(flexIdx), 1);
    for k = 1:numel(flexIdx)
        if seen(k)
            continue
        end
        [~, e] = ismember(protFlex{k}, model.ec.enzymes);
        c = families.owner(e);
        group = k;
        if c > 0
            cplx = families.complexes{c};
            [~, group] = ismember(cplx.enzymes(cplx.owned), protFlex);
            group = group(group > 0)';
        end
        seen(group) = true;
        if numel(group) < 2
            continue
        end
        withOld = group(oldConcs(group) > 0);
        if isempty(withOld)
            continue
        end
        factor = max(newConcs(withOld) ./ oldConcs(withOld));
        newConcs(withOld) = oldConcs(withOld) * factor;
    end
    keepOld = newConcs < oldConcs;
    model.ub(usageIdx(flexIdx(keepOld)))  = oldConcs(keepOld);
    model.ub(usageIdx(flexIdx(~keepOld))) = newConcs(~keepOld);
    protFlex(keepOld) = [];
    oldConcs(keepOld) = [];
    freq(keepOld)     = [];
    flexConcs         = newConcs(~keepOld);
end

flexEnz.uniprotIDs = protFlex;
flexEnz.oldConcs   = oldConcs;
flexEnz.flexConcs  = flexConcs;
flexEnz.frequence  = freq;
if changedProtPool
    flexEnz.uniprotIDs(end+1) = {'prot_pool'};
    flexEnz.oldConcs(end+1)   = oldProtPool;
    flexEnz.flexConcs(end+1)  = newProtPool;
    flexEnz.frequence(end+1)  = 1;
end
flexEnz.ratioIncr = flexEnz.flexConcs./flexEnz.oldConcs;
[~,b] = sort(flexEnz.ratioIncr,'descend');
flexEnz.ratioIncr  = flexEnz.ratioIncr(b);
flexEnz.uniprotIDs = flexEnz.uniprotIDs(b);
flexEnz.oldConcs   = flexEnz.oldConcs(b);
flexEnz.flexConcs  = flexEnz.flexConcs(b);
flexEnz.frequence  = flexEnz.frequence(b);
end

function scores = getUnitScores(model, usageIdx, units, foldChange, limit)
% Growth gained per fold of raising the usage upper bounds of each unit
% together. With reduced costs, the sum over the unit's subunits that run
% above limit of their cap of |reduced cost of the usage reaction| times its
% upper bound. Unlike the shadow price of the protein metabolite, the
% reduced cost of the usage reaction leaves out the price of the protein
% pool, which would otherwise be counted once per subunit. When no such
% score is positive (no reduced costs, or the protein pool binds and the
% solver puts the whole price on it), the candidate units are probed one by
% one with all their upper bounds raised by foldChange and the pool lifted.
scores = zeros(numel(units), 1);
sol = solveLP(model);
if sol.stat ~= 1
    return
end
candidates = zeros(0, 1);
active     = cell(numel(units), 1);
for k = 1:numel(units)
    rx = usageIdx(units{k});
    ub = model.ub(rx);
    act = ub > 0;
    act(act) = sol.x(rx(act)) ./ ub(act) > limit;
    if any(act)
        candidates(end+1, 1) = k; %#ok<AGROW>
        active{k} = rx(act);
    end
end
if isempty(candidates)
    return
end
needed = vertcat(active{candidates});
if ~isempty(sol.rCost) && all(isfinite(sol.rCost(needed)))
    for k = candidates'
        scores(k) = sum(abs(sol.rCost(active{k})) .* model.ub(active{k}));
    end
    scores(scores < 1e-12) = 0;
    if any(scores > 0)
        return
    end
    scores = zeros(numel(units), 1);
end
poolIdx = getIndexes(model, 'prot_pool_exchange', 'rxns');
model.ub(poolIdx) = max(model.ub(poolIdx), 1000);
base = solveLP(model);
if base.stat ~= 1
    return
end
for k = candidates'
    tempModel = model;
    rx = usageIdx(units{k});
    tempModel.ub(rx) = tempModel.ub(rx) * foldChange;
    probe = solveLP(tempModel);
    if probe.stat == 1
        gain = probe.f - base.f;
        if gain > 1e-3
            scores(k) = gain / (foldChange - 1);
        end
    end
end
end

