function [model, rxnsMissingGPR, standardMW, standardKcat, rxnsNoKcat] = getStandardKcat(model, varargin)
% getStandardKcat  Assign standard kcat and MW to reactions without enzymes.
%
% Calculates a standard kcat and standard molecular weight (MW) that can be
% used to apply enzyme constraints to reactions without any associated
% enzymes. Such reactions have either an empty model.grRules field, or they
% have no match in model.ec.rxns (which can be the case if the genes in the
% model.grRules field could not be mapped to enzymes). This is done by
% adding those reactions to model.ec, assigning a "standard" pseudoenzyme
% with the standard MW (median of all proteins in the organism) and standard
% kcat (median from all kcat, or subsystem specific kcat).
%
% A reaction is assigned a subSystem specific kcat value if the model has a
% subSystems field and the reaction is annotated with a subSystem. Only the
% first subSystem will be considered if multiple are annotated to the same
% reaction.
%
% Exchange (RAVEN's getExchangeRxns: no substrate or no product), transport
% and pseudoreactions (a reaction name that contains 'pseudoreaction' or
% 'SLIME rxn', in any letter case) are filtered out, plus any
% reaction identifiers specified in /data/pseudoRxns.tsv in the model
% adapter folder.
%
% In addition, reactions that are annotated with an enzyme (and therefore
% already in model.ec), but not assigned any reaction-specific kcat value
% (their model.ec.kcat entry is either 0 or NaN), can be assigned standard
% kcat values by a similar approach. However, those reactions will not be
% linked to the "standard" pseudoenzyme, but will use the enzyme that they
% had already been associated with.
%
% Any pre-existing standard kcat assignments (identified by 'standard'
% entries in model.ec.source) are removed when applying this function.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% threshold : double
%     a threshold to determine whether to use a kcat value based on the
%     median kcat of the reactions in the same subSystem or based on the
%     median value of all the kcat in the model. The second option is used
%     when the number of reactions in a determined subSystem with a real
%     (nonzero, non-NaN) kcat value is < threshold (default 10).
% fillZeroKcat : logical
%     whether zero kcat values should be replaced with standard kcat values
%     (default true).
%
% Returns
% -------
% model : struct
%     ecModel where model.ec is expanded with a standard protein with
%     standard kcat and standard MW, assigned to reactions without gene
%     associations.
% rxnsMissingGPR : cell
%     a list of updated rxns identifiers with a standard value.
% standardMW : double
%     the standard MW value calculated.
% standardKcat : double
%     the standard kcat value calculated.
% rxnsNoKcat : cell
%     a list of rxns identifiers whose zero kcat has been replaced.
%
% Notes
% -----
% While model.ec.kcat is populated, applyKcatConstraints would still need to
% be run to apply the new constraints to the S-matrix.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [model, rxnsMissingGPR, standardMW, standardKcat, rxnsNoKcat] = getStandardKcat(model);
%     [model, rxnsMissingGPR, standardMW, standardKcat, rxnsNoKcat] = getStandardKcat(model, 'threshold', 20);
%
% See also
% --------
% removeStandardKcat, applyKcatConstraints

p = parseGECKOargs(varargin, { ...
    'modelAdapter', []; ...
    'threshold',    []; ...
    'fillZeroKcat', []});
modelAdapter = p.modelAdapter;
threshold    = p.threshold;
fillZeroKcat = p.fillZeroKcat;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(threshold)
    threshold = 10;
end

if isempty(fillZeroKcat)
    fillZeroKcat = true;
end

databases = loadDatabases('uniprot', modelAdapter);

% An standard MW is defined for all the rxns which does not have a GPR
% rule defined. This is based in all the proteins reported for the specific
% organism in uniprot
standardMW = median(databases.uniprot.MW, 'omitnan');

% An standard Kcat is defined for all the rxns which does not have a GPR
% rule defined. In this case, the kcat value for a particular reaction is
% defined as the median of the kcat values of the reactions (with a real,
% nonzero, non-NaN kcat) involved in the same subsystem in which the given
% reaction is involved. Nevertheless, if a subSystem have a number of such
% reactions lower than a treshold, the kcat value will be the median of the
% kcat in all the reactions of the model.

% Remove from the list those with kcat zero
rxnsKcatZero = model.ec.kcat > 0;

% Get the kcat value based on all the kcats in the model
standardKcat = median(model.ec.kcat(rxnsKcatZero), 'omitnan');

% If the model have subSystems assigned calculate kcat based on subSystem
if isfield(model,'subSystems') && ~all(cellfun(@isempty, model.subSystems))
    standard = false;
    if model.ec.geckoLight
        modRxns = extractAfter(model.ec.rxns,4);
    else
        modRxns = model.ec.rxns;
    end
    % Map ec-rxns to model.rxns
    [~,rxnIdx]  = ismember(modRxns,model.rxns);
    % Choose first subSystem ('' when a reaction has none)
    enzSubSystems = cellfun(@firstSubSystem, model.subSystems(rxnIdx), 'UniformOutput', false);
    if ~all(cellfun(@isempty, enzSubSystems))

    % Make list of unique subsystems, and which rxns are linked to them.
    % Names that differ only in letter case are one subSystem.
    [~, firstOccurrence, rxnToSub] = unique(lower(enzSubSystems));
    enzSubSystem_names = enzSubSystems(firstOccurrence);
    % Calculate the median kcat per subSystem, counting only reactions
    % with a real (nonzero, non-NaN) kcat toward both the median itself
    % and the threshold that decides whether a subSystem-specific value
    % is used at all.
    numSubSystem = numel(enzSubSystem_names);
    kcatsPerSubSystem = zeros(numSubSystem,1);
    kcatSubSystem = zeros(numSubSystem,1);
    for i = 1:numSubSystem
        subKcat = model.ec.kcat(rxnToSub == i);
        subKcat = subKcat(subKcat ~= 0 & ~isnan(subKcat));
        kcatsPerSubSystem(i) = numel(subKcat);
        if kcatsPerSubSystem(i) >= threshold
            kcatSubSystem(i) = median(subKcat);
        else
            kcatSubSystem(i) = standardKcat;
        end
    end
    else
        standard = true;
        printOrange('WARNING: No subSystem-specific kcat values can be calculated');
    end
else
    standard = true;
    printOrange('WARNING: No subSystem-specific kcat values can be calculated');
end

% Find reactions without GPR
rxnsMissingGPR = find(cellfun(@isempty, model.grRules));

% Find reactions with GPR but without model.ec entry (for instance due to
% no protein matching)
rxnsMissingEnzyme = find(~cellfun(@isempty, model.grRules));
if model.ec.geckoLight
    ecRxnsList = unique(extractAfter(model.ec.rxns,4));
else
    ecRxnsList = model.ec.rxns;
end
rxnsMissingEnzyme = find(and(~ismember(model.rxns(rxnsMissingEnzyme),ecRxnsList), ~contains(model.rxns(rxnsMissingEnzyme),'usage_prot_')));
rxnsMissingGPR = [rxnsMissingGPR;rxnsMissingEnzyme];

% Get custom list of reaction IDs to ignore, if existing. First column
% contains reaction IDs, second column contains reaction names for
% reference only.
if exist(fullfile(params.path,'data','pseudoRxns.tsv'),'file')
    fID        = fopen(fullfile(params.path,'data','pseudoRxns.tsv'));
    fileData   = textscan(fID,'%s %s','delimiter','\t');
    fclose(fID);
    customRxns = fileData{1};
    customRxns = find(ismember(model.rxns,customRxns));
else
    customRxns = [];
end
% Get and remove exchange, transport, spontaneous and pseudo reactions
[~, exchangeRxns]  = getExchangeRxns(model);
transportRxns = getTransportRxns(model);
[spontaneousRxns, ~] = modelAdapter.getSpontaneousReactions(model);
pseudoRxns = contains(lower(model.rxnNames),'pseudoreaction');
slimeRxns = contains(lower(model.rxnNames),'slime rxn');
rxnsToIgnore = [customRxns; exchangeRxns; find(transportRxns); ...
                find(spontaneousRxns); find(pseudoRxns); find(slimeRxns)];
rxnsMissingGPR(ismember(rxnsMissingGPR, rxnsToIgnore)) = [];

% Only add if getStandardKcat was not run earlier
if ~any(strcmp(model.mets,'prot_standard')) && ~any(strcmp(model.ec.enzymes,'standard'))
    % Add a new gene to be consistent with ec field named standard
    proteinStdGenes.genes = 'standard';
    if isfield(model,'geneShortNames')
        proteinStdGenes.geneShortNames = 'std';
    end
    model = addGenesRaven(model, proteinStdGenes);

    if ~model.ec.geckoLight
        % Add a new metabolite named prot_standard
        proteinStdMets.mets         = 'prot_standard';
        proteinStdMets.metNames     = proteinStdMets.mets;
        % Validate compartment
        compartmentID = strcmp(model.compNames,params.enzyme_comp);
        if ~any(compartmentID)
            error(['Compartment ' params.enzyme_comp ' (specified in params.enzyme_comp) '...
                'cannot be found in model.compNames'])
        end
        proteinStdMets.compartments = model.comps(compartmentID);

        if isfield(model,'metNotes')
            proteinStdMets.metNotes = 'Standard enzyme-usage pseudometabolite';
        end
        model = addMets(model, proteinStdMets);

        % Add a protein usage reaction if not a light version
        proteinStdUsageRxn.rxns         = {'usage_prot_standard'};
        proteinStdUsageRxn.rxnNames     = proteinStdUsageRxn.rxns;
        proteinStdUsageRxn.mets         = {proteinStdMets.mets, 'prot_pool'};
        proteinStdUsageRxn.stoichCoeffs = [1, -1];
        proteinStdUsageRxn.lb           = 0;
        proteinStdUsageRxn.ub           = 1000;
        proteinStdUsageRxn.grRules      = proteinStdGenes.genes;

        model = addRxns(model, proteinStdUsageRxn);
    end
    % Update .ec structure in model
    model.ec.genes(end+1)      = {'standard'};
    model.ec.enzymes(end+1)    = {'standard'};
    model.ec.mw(end+1)         = standardMW;
    model.ec.sequence(end+1)   = {''};
    % Additional info
    if isfield(model.ec,'concs')
        model.ec.concs(end+1)  = nan();
    end

    % Expand the enzyme rxns matrix
    model.ec.rxnEnzMat =  [model.ec.rxnEnzMat, zeros(length(model.ec.rxns), 1)]; % 1 new enzyme
    model.ec.rxnEnzMat =  [model.ec.rxnEnzMat; zeros(length(rxnsMissingGPR), length(model.ec.enzymes))]; % new rxns
end
stdMetIdx = find(strcmpi(model.ec.enzymes, 'standard'));

% Remove previous standard kcat assignment
oldStandardEnz = find(strcmp(model.ec.source,'standard'));
if ~isempty(oldStandardEnz)
    oldStandardProt = logical(model.ec.rxnEnzMat(oldStandardEnz,stdMetIdx));
    % If annotated with real enzyme => set kcat to zero
    model.ec.kcat(oldStandardEnz(~oldStandardProt))        = 0;
    model.ec.source(oldStandardEnz(~oldStandardProt))      = {''};
    % If annotated with standard protein => remove entry
    model.ec.rxns(oldStandardEnz(oldStandardProt))        = [];
    model.ec.kcat(oldStandardEnz(oldStandardProt))        = [];
    model.ec.source(oldStandardEnz(oldStandardProt))      = [];
    model.ec.notes(oldStandardEnz(oldStandardProt))       = [];
    model.ec.eccodes(oldStandardEnz(oldStandardProt))     = [];
    model.ec.rxnEnzMat(oldStandardEnz(oldStandardProt),:) = [];
end

numRxns = length(model.ec.rxns);
for i = 1:numel(rxnsMissingGPR)
    rxnIdx = rxnsMissingGPR(i);

    % Update .ec structure in model
    if ~model.ec.geckoLight
        model.ec.rxns(end+1)     = model.rxns(rxnIdx);
        % Add prefix in case is light version
    else
        model.ec.rxns{end+1}     = ['001_' model.rxns{rxnIdx}];
    end

    if ~standard
        % kcatSubSystemIdx is a one-hot vector (the reaction's first subSystem
        % matched against every subSystem that has a computed kcat); any(), not
        % all(), is the right test -- all() only ever succeeds when the model
        % has exactly one subSystem in total, silently falling back to the
        % model-wide standardKcat for every reaction otherwise.
        rxnSubSystem = firstSubSystem(model.subSystems{rxnIdx});
        kcatSubSystemIdx = strcmpi(enzSubSystem_names, rxnSubSystem);
        if ~isempty(rxnSubSystem) && any(kcatSubSystemIdx)
            model.ec.kcat(end+1) = kcatSubSystem(kcatSubSystemIdx);
        else
            model.ec.kcat(end+1) = standardKcat;
        end
    else
        model.ec.kcat(end+1) = standardKcat;
    end

    model.ec.source(end+1)   = {'standard'};
    model.ec.notes(end+1)    = {''};
    model.ec.eccodes(end+1)  = {''};

    % Update the enzyme rxns matrix
    model.ec.rxnEnzMat(numRxns+i, stdMetIdx) = 1;
end
% Get the rxns identifiers of the updated rxns
rxnsMissingGPR = model.rxns(rxnsMissingGPR);

if fillZeroKcat
    zeroKcat = model.ec.kcat == 0 | isnan(model.ec.kcat);
    model.ec.kcat(zeroKcat)     = standardKcat;
    model.ec.source(zeroKcat)   = {'standard'};
    rxnsNoKcat = model.ec.rxns(zeroKcat);
else
    rxnsNoKcat = [];
end
end

function name = firstSubSystem(subSystems)
% First subSystem of a reaction as a character vector; '' when the reaction
% has none.
if ischar(subSystems)
    name = subSystems;
elseif iscell(subSystems) && ~isempty(subSystems) && ischar(subSystems{1})
    name = subSystems{1};
else
    name = '';
end
end
