function pairs = extractEnzymeSubstratePairs(model, varargin)
% extractEnzymeSubstratePairs  The enzyme-substrate pairs of an ecModel's reactions.
%
% Lists, for each model.ec reaction, its enzymes and its substrates as
% DLKcat and the homolog search use them: the reactants of the reaction,
% without the metabolites of DLKcatIgnoreMets.tsv (metal ions, small
% molecules), without protein pseudometabolites, and without currency
% metabolites of DLKcatCurrencyMets.tsv where both members of a pair take
% part (unless no other substrate remains). Both files are read from the
% data folder of the project, or else from the databases folder of GECKO.
% geckopy has the same (extract_enzyme_substrate_pairs).
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% ecRxns : logical
%     which model.ec.rxns to list (default: all).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% pairs : struct
%     ecRxn (index into model.ec.rxns), protein (index into model.ec.genes,
%     model.ec.enzymes and model.ec.sequence) and met (index into
%     model.mets), one entry per pair, reaction by reaction.
%
% Examples
% --------
%     pairs = extractEnzymeSubstratePairs(ecModel);
%     model.metNames(pairs.met(1:5))
%
% See also
% --------
% writeDLKcatInput, kcatsFromHomologs

[geckoPath, ~] = findGECKOroot();
p = parseGECKOargs(varargin, {'ecRxns', []; 'modelAdapter', []});
ecRxnsSel = p.ecRxns;
if isempty(ecRxnsSel)
    ecRxnsSel = true(numel(model.ec.rxns),1);
elseif ~islogical(ecRxnsSel)
    error('ecRxns should be provided as logical vector')
elseif numel(ecRxnsSel)~=numel(model.ec.rxns)
    error('Length of ecRxns is not the same as model.ec.rxns')
end
ecRxnsSel = find(ecRxnsSel);
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if ~model.ec.geckoLight
   origRxns = model.ec.rxns;
else
   origRxns = extractAfter(model.ec.rxns,4);
end
origRxnsToInclude = origRxns(ecRxnsSel);

% Map back to original reactions, to extract substrates
[sanityCheck,origRxnIdxs] = ismember(origRxnsToInclude,model.rxns);
if ~all(sanityCheck)
    error('Not all reactions in model.ec.rxns are found in model.rxns')
end

% Ignore selected metabolites (metal ions, proteins etc.). First check by
% name (case insensitive, without white spaces and special characters),
% then also try to match with metSmiles (if available).
metsNoSpecialChars = lower(regexprep(model.metNames,'[^0-9a-zA-Z]+',''));
if exist(fullfile(params.path,'data','DLKcatIgnoreMets.tsv'),'file')
    fID        = fopen(fullfile(params.path,'data','DLKcatIgnoreMets.tsv'));
else
    fID        = fopen(fullfile(geckoPath,'databases','DLKcatIgnoreMets.tsv'));
end
fileData   = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[ignoreMets, ignoreSmiles] = deal(fileData{[1,2]});
ignoreMets = lower(regexprep(ignoreMets,'[^0-9a-zA-Z]+',''));
ignoreSmiles(cellfun(@isempty,ignoreSmiles)) = [];

ignoreMetsIdx  = logical(ismember(metsNoSpecialChars,ignoreMets));
if isfield(model,'metSmiles')
    ignoreMetsIdx = ignoreMetsIdx | logical(ismember(model.metSmiles,ignoreSmiles));
end
% Also leave out protein-usage pseudometabolites
ignoreMetsIdx = ignoreMetsIdx | startsWith(model.mets,'prot_');
reducedS = model.S;
reducedS(ignoreMetsIdx,:) = 0;

% Ignore currency metabolites if they occur in pairs. First check by
% name (case insensitive, without white spaces and special characters),
% then also try to match with metSmiles (if available).
if exist(fullfile(params.path,'data','DLKcatCurrencyMets.tsv'),'file')
    fID = fopen(fullfile(params.path,'data','DLKcatCurrencyMets.tsv'));
else
    fID = fopen(fullfile(geckoPath,'databases','DLKcatCurrencyMets.tsv'));
end
fileData = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);
[currencyMets(:,1), currencyMets(:,2)] = deal(fileData{[1,2]});
currencyMets = lower(regexprep(currencyMets,'[^0-9a-zA-Z]+',''));

for i=1:size(currencyMets,1)
    subs = strcmp(currencyMets(i,1),metsNoSpecialChars);
    prod = strcmp(currencyMets(i,2),metsNoSpecialChars);
    [~,subsRxns]=find(reducedS(subs,:));
    [~,prodRxns]=find(reducedS(prod,:));
    pairRxns = intersect(subsRxns,prodRxns);
    tempRedS=reducedS;
    tempRedS([find(subs);find(prod)],pairRxns) = 0;
    % Do not remove currency mets if no substrate remains
    rxnsWithRemainingSubstrates = any(tempRedS(:,pairRxns) < 0,1);
    reducedS([find(subs);find(prod)],intersect(pairRxns,pairRxns(rxnsWithRemainingSubstrates))) = 0;
end

% origRxnIdxs already selects exactly the model.rxns columns for the requested
% ecRxns subset (one column per requested reaction, via origRxns(ecRxns) above), so
% clearedRedS needs no further filtering here.
clearedRedS = reducedS(:,origRxnIdxs);

% Enumerate all substrates for each reaction. clearedRedS has one column per entry
% of ecRxnsSel, in that same order, so `reactions` here indexes *that subset*, not
% model.ec.rxns directly --- translate back to a real model.ec.rxns index before
% using it against model.ec.rxns / model.ec.rxnEnzMat, which are always full-length.
[substrates, reactions] = find(clearedRedS<0);
globalEcRxns = ecRxnsSel(reactions);

% Enumerate all proteins for each reaction
[proteins, enzRows] = find(transpose(model.ec.rxnEnzMat(globalEcRxns,:)));
pairs.ecRxn   = reshape(globalEcRxns(enzRows), [], 1);
pairs.protein = reshape(proteins, [], 1);
pairs.met     = reshape(substrates(enzRows), [], 1);
end
