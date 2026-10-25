function [model, rxnsAdded, enzAdded] = addNewRxnsToEC(model, newRxns, varargin)
% addNewRxnsToEC  Add new reactions to an enzyme-constrained model.
%
% Adds new reactions to an enzyme-constrained model. This function is
% useful to simulate metabolic manipulations done to an organism such as
% integration of new genes that add a new reaction/pathway.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% newRxns : struct
%     structure with the new reaction information (see Notes for fields).
%
% Name-Value Arguments
% --------------------
% newEnzymes : struct
%     structure with the new enzymes information (see Notes for fields).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% model : struct
%     ecModel with new reactions.
% rxnsAdded : cell
%     cell array with the reactions added.
% enzAdded : cell
%     cell array with the enzymes added.
%
% Notes
% -----
% The newRxns structure has the following fields:
%
% - rxns : cell array with unique strings that identifies each reaction.
% - rxnNames : cell array with the names of each reaction.
% - equations : cell array with equation strings. Decimal coefficients are expressed as "1.2". Reversibility is indicated by "<=>" or "=>".
% - grRules : cell array with the gene-reaction relationship for each reaction. For example "(A and B) or (C)" means that the reaction could be catalyzed by a complex between A & B or by C on its own. A rule with brackets is split into its isozymes in disjunctive normal form, so "A and (B or C)" is the two isozymes "A and B" and "A and C". All the genes have to be present in model.genes. Add genes with addGenesRaven before calling this function if needed (opt, default '').
%
% The newEnzymes structure has the following fields:
%
% - enzymes : cell array with uniprot id.
% - genes : cell array with the respective gene.
% - mw : cell array with the MW.
%
% Additional usage notes:
%
% 1. Write equations as follows:
%        'xylitol[c] + NAD[c] => D-xylulose[c] + NADH[c] + H+[c]'
%    Note that the reversibility defined in the equation will be used to
%    split to construct irreversible reactions (add _REV) and the rules
%    defined will be used to expand the model (add _EXP_n).
% 2. After adding the new reactions, setKcatForReactions or
%    applyCustomKcats should be run to add the kcat to the reactions, and
%    subsequently applyKcatConstraints.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [model, rxnsAdded, enzAdded] = addNewRxnsToEC(model, newRxns, newEnzymes, modelAdapter);
%     [model, rxnsAdded, enzAdded] = addNewRxnsToEC(model, newRxns, 'newEnzymes', newEnzymes);
%
% See also
% --------
% applyKcatConstraints, applyCustomKcats

p = parseGECKOargs(varargin, { ...
    'newEnzymes',   []; ...
    'modelAdapter', []});
newEnzymes   = p.newEnzymes;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(newEnzymes)
    newEnzymes = [];
end

% Validate if GECKO version 3
if ~isfield(model, 'ec')
    error('Input model is not an ecModel in GECKO 3 format (with ecModel.ec structure)')
end

if model.ec.geckoLight
    error('addNewRxnsToEC only works on non-light ecModels')
else
    % Check if the new proteins are not already present in the model
    if ~isempty(newEnzymes)
        [toRemove,~] = ismember(newEnzymes.enzymes, model.ec.enzymes);

        % Remove from the list those enzymes that already present
        if any(toRemove)
            printOrange(['WARNING: Enzymes ' strjoin(newEnzymes.enzymes(toRemove),', ') ' are already present in the model and will not be added.\n']);
            newEnzymes.enzymes(toRemove) = [];
            newEnzymes.genes(toRemove) = [];
            newEnzymes.mw(toRemove) = [];
        end

        if ~isempty(newEnzymes.enzymes)
            
            % 1. Add new genes
            genesToAdd.genes = newEnzymes.genes;
            genesToAdd.geneShortNames = genesToAdd.genes;
            model = addGenesRaven(model, genesToAdd);

            % 2. Add new Enzyme-usage pseudometabolite

            % Validate compartment to add protein pseudoreactions
            compartmentID = strcmp(model.compNames,params.enzyme_comp);
            if ~any(compartmentID)
                error(['Compartment ' params.enzyme_comp ' (specified in params.enzyme_comp) '...
                    'cannot be found in model.compNames'])
            end
            compartmentID = model.comps(compartmentID);

            proteinMets.mets         = strcat('prot_', newEnzymes.enzymes);
            proteinMets.metNames     = proteinMets.mets;
            proteinMets.compartments = compartmentID;
            if isfield(model,'metMiriams')
                proteinMets.metMiriams = repmat({struct('name',{{'sbo'}},'value',{{'SBO:0000252'}})},numel(proteinMets.mets),1);
            end
            if isfield(model,'metCharges')
                proteinMets.metCharges = zeros(numel(proteinMets.mets),1);
            end
            proteinMets.metNotes = repmat({'Enzyme-usage pseudometabolite'},numel(proteinMets.mets),1);
            model = addMets(model,proteinMets);
        end
        enzAdded = newEnzymes.enzymes(:);  
    end

    % 3. Add new reactions

    % Before add the reactions validate that all the genes in grRules are
    % present in the model.
    genesInRules = strjoin(convertCharArray(newRxns.grRules));
    genesInRules = regexp(genesInRules,'[^\s()]+','match'); % Gene ids and the words 'and' and 'or', without spaces and brackets
    genesInRules = genesInRules(~ismember(lower(genesInRules),{'and','or'}));  % Only the gene ids
    genesInRules = setdiff(unique(genesInRules),model.genes); % Only keep new genes

    % If there is any gene which is not present, as to be pass as
    % newEnzymes input
    if ~isempty(genesInRules)
        errorText = 'The following genes from grRules are not present in the model, add them as newEnzymes input:';
        error('%s', ravenList(errorText,genesInRules,false))
    end

    % Get reversible reactions to split
    toIrrev = find(cellfun(@(x) contains(x,'<=>'),newRxns.equations));
    if ~isempty(toIrrev)
        for i=1:numel(toIrrev)
            idx = toIrrev(i);
            newEquation = split(newRxns.equations(idx),'<=>');
            % Remove directionality
            newRxns.equations(idx) = erase(newRxns.equations(idx), '<');

            % Add a new reaction to the list to be added
            newRxns.rxns{end+1} = [newRxns.rxns{idx} '_REV'];
            newRxns.rxnNames{end+1} = [newRxns.rxnNames{idx} ' (reversible)'];
            newRxns.equations{end+1} = [strtrim(newEquation{2, 1}) ' => ' strtrim(newEquation{1, 1})];
            newRxns.grRules{end+1} = newRxns.grRules{idx};
        end
    end

    % Get reactions to separate isozymes
    isozymes = cellfun(@gprToDnf, newRxns.grRules, 'UniformOutput', false);
    toSplit = find(cellfun(@numel,isozymes) > 1);
    if ~isempty(toSplit)
        splitRxns = {}; splitRxnNames = {}; splitEquations = {}; splitGrRules = {};
        for i=1:numel(toSplit)
            idx = toSplit(i);
            % One reaction per AND-clause of the gene rule in disjunctive normal form
            newRules = cellfun(@(x) strjoin(x,' and '), isozymes{idx}, 'UniformOutput', false);

            for j=1:numel(newRules)
                splitRxns{end+1} = [newRxns.rxns{idx} '_EXP_' num2str(j)];
                splitRxnNames(end+1) = newRxns.rxnNames(idx);
                splitEquations(end+1) = newRxns.equations(idx);
                splitGrRules(end+1) = newRules(j);
            end
        end

        % Remove every original that was split, in one indexed operation, then append
        % the replacements. Removing (or appending) one at a time while iterating a list
        % of indices computed up front lets each removal shift the indices of everything
        % after it, so a second (or later) original silently gets read from the wrong
        % position once more than one reaction needs splitting in the same call.
        newRxns.rxns(toSplit) = [];
        newRxns.rxnNames(toSplit) = [];
        newRxns.equations(toSplit) = [];
        newRxns.grRules(toSplit) = [];

        newRxns.rxns = [newRxns.rxns(:); splitRxns(:)];
        newRxns.rxnNames = [newRxns.rxnNames(:); splitRxnNames(:)];
        newRxns.equations = [newRxns.equations(:); splitEquations(:)];
        newRxns.grRules = [newRxns.grRules(:); splitGrRules(:)];
    end

    newRxns.lb = zeros(length(newRxns.rxns),1);
    newRxns.ub = zeros(length(newRxns.rxns),1) + 1000;

    model = addRxns(model, newRxns, 3, '', true, true);

    % 4. Add protein usage reactions.
    usageRxns.rxns            = strcat('usage_',proteinMets.mets);
    usageRxns.rxnNames        = usageRxns.rxns;
    usageRxns.mets            = cell(numel(usageRxns.rxns),1);
    usageRxns.stoichCoeffs    = cell(numel(usageRxns.rxns),1);
    for i=1:numel(usageRxns.mets)
        usageRxns.mets{i}         = {proteinMets.mets{i}, 'prot_pool'};
        usageRxns.stoichCoeffs{i} = [1,-1];
    end
    usageRxns.lb              = zeros(numel(usageRxns.rxns),1);
    usageRxns.ub              = zeros(numel(usageRxns.rxns),1) + 1000;
    usageRxns.rev             = ones(numel(usageRxns.rxns),1);
    usageRxns.grRules         = newEnzymes.genes;
    model = addRxns(model,usageRxns);

    % 5. Update .ec field with enzyme data
    model.ec.enzymes    = [model.ec.enzymes; newEnzymes.enzymes(:)];
    model.ec.genes      = [model.ec.genes; newEnzymes.genes(:)];
    model.ec.mw         = [model.ec.mw; newEnzymes.mw(:)];
    model.ec.sequence   = [model.ec.sequence; repmat({''},numel(newEnzymes.enzymes),1)]; 
    if isfield(model.ec,'concs')
        model.ec.concs  = [model.ec.concs; NaN(numel(newEnzymes.enzymes),1)];
    end

    % 6. Expand the enzyme rxns matrix
    rxnWithRule = find(~cellfun('isempty',newRxns.grRules));

    model.ec.rxnEnzMat =  [model.ec.rxnEnzMat, zeros(length(model.ec.rxns), length(newEnzymes.enzymes))]; % add n new enzymes
    model.ec.rxnEnzMat =  [model.ec.rxnEnzMat; zeros(length(rxnWithRule), length(model.ec.enzymes))]; % add n new rxns

    % 7. Update .ec field with rxn data

    % Get reactions to add 
    numRxns = length(model.ec.rxns);
    for i = 1:numel(rxnWithRule)
        rxnIdx = rxnWithRule(i);

        model.ec.rxns(end+1)    = newRxns.rxns(rxnIdx);
        model.ec.kcat(end+1)    = 0;
        model.ec.source(end+1)  = {''};
        model.ec.notes(end+1)   = {''};
        model.ec.eccodes(end+1) = {''};

        % Get rules
        genes = regexp(newRxns.grRules{rxnIdx},'[^\s()]+','match'); % Gene ids and the word 'and', without spaces and brackets
        genes = genes(~ismember(lower(genes),{'and','or'}));
        [~,enzIdx] = ismember(genes, model.ec.genes);

        % Update the enzyme rxns matrix
        model.ec.rxnEnzMat(numRxns+i, enzIdx) = 1;
    end
    
    rxnsAdded = newRxns.rxns(:);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that expands a gene rule into disjunctive normal form: a cell array
%with one cell per isozyme, each holding the gene ids that are ANDed in it.
%Gene ids are whole words, so an id that contains 'or' or 'and' is kept
%intact. An empty rule gives an empty cell array.
function dnf = gprToDnf(rule)
tokens = regexp(rule,'\(|\)|[^\s()]+','match');
isOperator = strcmpi(tokens,'and') | strcmpi(tokens,'or');
tokens(isOperator) = lower(tokens(isOperator)); % Operators are case insensitive
if isempty(tokens)
    dnf = {};
    return
end
[dnf,pos] = parseOr(tokens,1,rule);
if pos <= numel(tokens)
    error('Cannot parse the gene rule "%s"',rule)
end
end

function [dnf,pos] = parseOr(tokens,pos,rule)
[dnf,pos] = parseAnd(tokens,pos,rule);
while pos <= numel(tokens) && strcmp(tokens{pos},'or')
    [rhs,pos] = parseAnd(tokens,pos+1,rule);
    dnf = [dnf,rhs];
end
end

function [dnf,pos] = parseAnd(tokens,pos,rule)
[dnf,pos] = parseAtom(tokens,pos,rule);
while pos <= numel(tokens) && strcmp(tokens{pos},'and')
    [rhs,pos] = parseAtom(tokens,pos+1,rule);
    product = cell(1,0);
    for a = 1:numel(dnf)
        for b = 1:numel(rhs)
            product{end+1} = [dnf{a},rhs{b}]; %#ok<AGROW>
        end
    end
    dnf = product;
end
end

function [dnf,pos] = parseAtom(tokens,pos,rule)
if pos > numel(tokens) || any(strcmp(tokens{pos},{')','and','or'}))
    error('Cannot parse the gene rule "%s"',rule)
elseif strcmp(tokens{pos},'(')
    [dnf,pos] = parseOr(tokens,pos+1,rule);
    if pos > numel(tokens) || ~strcmp(tokens{pos},')')
        error('Cannot parse the gene rule "%s"',rule)
    end
    pos = pos + 1;
else
    dnf = {tokens(pos)};
    pos = pos + 1;
end
end
