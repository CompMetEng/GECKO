function fseof = ecFSEOF(model,prodTargetRxn,csRxn,varargin)
% ecFSEOF  Run Flux-Scanning with Enforced Objective Function for a target.
%
% Runs Flux-Scanning with Enforced Objective Function (FSEOF) for a
% specified production target, to identify reaction and gene targets for
% metabolic engineering. Reactions are classified as targets with the same
% rule as RAVEN's FSEOF: the correlation between the absolute flux of a
% reaction and the enforced production flux.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% prodTargetRxn : char
%     rxn ID for the production target reaction, an exchange reaction is
%     recommended.
% csRxn : char
%     rxn ID for the main carbon source uptake reaction.
%
% Name-Value Arguments
% --------------------
% nSteps : double
%     number of steps for suboptimal objective in FSEOF (default 16).
% outputFile : logical
%     boolean option to save results in a file (default false).
% filePath : char
%     file path for results output. It will store two files: at the genes
%     level, ecFSEOF_genes.tsv; and at the reactions level,
%     ecFSEOF_rxns.tsv (default: the 'output' sub-folder taken from
%     modelAdapter, e.g. output/ecFSEOF_rxns.tsv).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% corrThreshold : double
%     minimum absolute Pearson correlation between the absolute flux of a
%     reaction and the enforced production flux for the reaction to be a
%     target (default 0.9).
%
% Returns
% -------
% fseof : struct
%     a structure with all results.
%
% Notes
% -----
% The fseof structure has the following fields:
%
% - alpha : target production used for enforced objective limits (from minimum to maximum production).
% - v_matrix : fluxes for each target reaction predicted and each alpha.
% - rxnTargets : a list with all reactions whose absolute flux is correlated with the enforced production, with an absolute correlation of at least corrThreshold. Contains ID, name, slope, gene rule, and equation. The slope is the absolute slope of the regression of the absolute flux on the enforced production flux. The list is ordered by the product of slope and correlation, highest first.
% - transportTargets : the same for transport reactions.
% - geneTargets : a list with all genes of the targets. Contains gene, shortName, slope, action, and essentiality.
%
% Only reactions with a gene rule are considered, and not those that are
% assigned to the standard enzyme. The enforced production reaction is
% included when it has a gene rule, as its flux follows the enforced
% production by definition. The reactions that draw enzyme usage
% (usage_prot_*) take part in the classification and are not reported.
% A reaction is over expressed (OE) when its flux increases with production,
% knocked out (KO) when it decreases and is at most 1e-8 at the last step,
% and knocked down (KD) when it decreases otherwise.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     fseof = ecFSEOF(model,prodTargetRxn,csRxn,nSteps,outputFile,filePath,modelAdapter,corrThreshold);
%     fseof = ecFSEOF(model,prodTargetRxn,csRxn,'nSteps',32);

p = parseGECKOargs(varargin, { ...
    'nSteps',       16; ...
    'outputFile',   false; ...
    'filePath',     []; ...
    'modelAdapter', []; ...
    'corrThreshold',0.9});
nSteps       = p.nSteps;
outputFile   = p.outputFile;
filePath     = p.filePath;
modelAdapter  = p.modelAdapter;
corrThreshold = p.corrThreshold;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(filePath)
    filePath = fullfile(params.path,'output');
end

% Get relevant rxn indexes
prodTargetRxnIdx = getIndexes(model, prodTargetRxn,'rxns');
csRxnIdx         = getIndexes(model, csRxn,'rxns');
bioRxnIdx        = getIndexes(model, params.bioRxn,'rxns');

% Standardize grRules and rxnGeneMat in model
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules        = grRules;
model.rxnGeneMat     = rxnGeneMat;

% Check carbon source uptake rate and LB defined
model = setParam(model, 'obj', params.bioRxn, 1);
sol   = solveLP(model);
if model.lb(csRxnIdx) < sol.x(csRxnIdx)
    printOrange(['WARNING: Carbon source lower bound was set to ' num2str(model.lb(csRxnIdx)) ...
        ', but the uptake rate after model optimization is ' num2str(sol.x(csRxnIdx)) '.\n']);
end

% run FSEOF analysis

% Find out the initial target production.
iniTarget = sol.x(prodTargetRxnIdx);

% Find out the maximum theoretical yield of target reaction.
model     = setParam(model,'obj',prodTargetRxn,1);
sol       = solveLP(model,1);
% Set to 90%% based on https://doi.org/10.1128/AEM.00115-10
maxTarget = sol.x(prodTargetRxnIdx) * 0.9;

% Define alpha vector for suboptimal enforced objective values between
% minimal production and 90%% of the maximum theoretical yield, initialize
% fluxes matriz
alpha    = iniTarget:((maxTarget-iniTarget)/(nSteps-1)):maxTarget;
v_matrix = zeros(length(model.rxns),length(alpha));

% Enforce objective flux iteratively
PB = progressReport(nSteps, 'Flux Scanning with Enforced Objective Function');
for i = 1:nSteps
    PB.count;
    % Enforce the objetive flux of product formation
    model = setParam(model,'eq',prodTargetRxnIdx,alpha(i));

    % Restore minimum biomass lb to zero and set it as objective
    model.lb(bioRxnIdx) = 0;
    model = setParam(model,'obj',params.bioRxn,1);
    sol   = solveLP(model,1);

    % Store flux distribution
    v_matrix(:,i) = sol.x;

end
PB.done;

% Take out rxns with no grRule and standard gene
withGR         = ~cellfun(@isempty,model.grRules);
stdIdx         = contains(model.grRules,'standard');
withGR(stdIdx) = 0;
target_rxns    = model.rxns(withGR);
v_matrix       = v_matrix(withGR,:);
rxnGeneM       = model.rxnGeneMat(withGR,:);

% Classify the reactions as in RAVEN's FSEOF. Slope and correlation are taken
% on the absolute flux against the enforced production flux, so a reversible
% reaction running backwards counts by its magnitude. A reaction is a target
% when some absolute flux exceeds 1e-10, the standard deviation of its
% absolute flux is above 0 and the absolute Pearson correlation reaches
% corrThreshold. It is suggested for over expression (OE) when the slope is
% positive, for knock out (KO) when the slope is negative and the absolute
% flux at the last step is at most 1e-8, and for knock down (KD) for any
% other negative slope.
absV        = abs(v_matrix);
slope_rxns  = zeros(size(target_rxns));
corr_rxns   = zeros(size(target_rxns));
target_type = cell(size(target_rxns));
targets     = false(size(target_rxns));
if numel(alpha) > 1 && std(alpha) > 0
    for i = 1:length(target_rxns)
        track = absV(i,:);
        if any(track > 1e-10) && std(track) > 0
            fit = polyfit(alpha,track,1);
            r   = corrcoef(alpha,track);
            if abs(r(1,2)) >= corrThreshold && fit(1) ~= 0
                targets(i)     = true;
                slope_rxns(i)  = abs(fit(1));
                corr_rxns(i)   = abs(r(1,2));
                if fit(1) > 0
                    target_type(i) = {'OE'};
                elseif track(end) <= 1e-8
                    target_type(i) = {'KO'};
                else
                    target_type(i) = {'KD'};
                end
            end
        end
    end
end

% Only keep those reactions that are classified as a target
target_rxns = target_rxns(targets);
v_matrix    = v_matrix(targets,:);
rxnGeneM    = rxnGeneM(targets,:);
slope_rxns  = slope_rxns(targets);
corr_rxns   = corr_rxns(targets);
target_type = target_type(targets);

% Order from highest to lowest score, the product of slope and correlation
[~,order]   = sort(slope_rxns .* corr_rxns,'descend');
target_rxns = target_rxns(order);
v_matrix    = v_matrix(order,:);
rxnGeneM    = rxnGeneM(order,:);
slope_rxns  = slope_rxns(order);
target_type = target_type(order);

% Create gene list of those connected to the remaining rxns
genes             = model.genes(sum(rxnGeneM,1) > 0);
slope_genes       = zeros(size(genes));
rxnGeneM          = rxnGeneM(:,sum(rxnGeneM,1) > 0);
target_type_genes = cell(size(genes));
essentiality      = cell(size(genes));

% Validate for gene essentiality
PB = progressReport(numel(genes), 'Checking for gene essentiality');
for i = 1:numel(genes)
    PB.count;

    % Block protein usage to KO al the reactions associated to it.
    usage_rxn_idx = strcmpi(model.ec.genes, genes{i});
    usage_rxn = strcat('usage_prot_', model.ec.enzymes(usage_rxn_idx));
    tempModel = setParam(model, 'eq', usage_rxn, 0);
    solKO     = solveLP(tempModel);
    % Check if no feasible solution was found
    if solKO.stat == -1 || solKO.x(bioRxnIdx) < 1e-8
        essentiality(i) = {'essential'};
    end

    % Since a gene can be involved in multiple reactions, multiple
    % engineering manipulations can be suggested for the same gene
    % e.g. (OE and KD). So, reconcilie them, and report only one.
    
    % Get reaction index, in targets set, for each gene
    rxns_for_gene = find(rxnGeneM(:,i) > 0);
    actions = unique(target_type(rxns_for_gene));
    if numel(actions) > 1
        % If any OE is suggested give the highest priority over KD or KO
        if any(startsWith(actions, 'OE'))
            target_type_genes(i) = {'OE'};
            % Otherwise change to KO if not essential
        elseif ~any(startsWith(actions, 'OE')) && ~isequal(essentiality(i),{'essential'})
            target_type_genes(i) = {'KO'};
        else
            target_type_genes(i) = {'KD'};
        end
    else
        target_type_genes(i) = actions;
    end
    % Extract all the slope (from rxns across alphas) conected to
    % each remaining gene
    slope_set      = slope_rxns(rxns_for_gene);
    slope_genes(i) = mean(slope_set);

end
PB.done;

% Order genes from highest to lowest slope:
[~,order]         = sort(slope_genes,'descend');
genes             = genes(order,:);
slope_genes       = slope_genes(order,:);
target_type_genes = target_type_genes(order,:);
essentiality      = essentiality(order,:);

% Create output:

% Report values used to enforce objective target production
fseof.alpha = alpha;

% Exclude enzyme usage reactions
toKeep = ~startsWith(target_rxns,'usage_prot_');
rxnIdx = getIndexes(model,target_rxns(toKeep),'rxns');

% Report the fluxes at each alpha
fseof.v_matrix = array2table(v_matrix(toKeep,:), ...
    'VariableNames', string(alpha), 'RowNames', model.rxns(rxnIdx));

% Report reaction targets
fseof.rxnTargets(:,1) = model.rxns(rxnIdx);
fseof.rxnTargets(:,2) = model.rxnNames(rxnIdx);
fseof.rxnTargets(:,3) = num2cell(slope_rxns(toKeep));
fseof.rxnTargets(:,4) = model.grRules(rxnIdx);
fseof.rxnTargets(:,5) = constructEquations(model,rxnIdx);

% Identify transport reactions; defined as reactions involving (at least) 
% one metabolite name in more than one compartment. 
transportRxns = false(numel(rxnIdx),1);
for i = 1:numel(rxnIdx)
    % Get the involved metabolites in each reaction
    mets = model.metNames(model.S(:,rxnIdx(i))~=0);
    % Remove enzyme metabolite
    mets(startsWith(mets,'prot_')) = [];
    % Validate if it is a transport reaciton
    transportRxns(i) = numel(mets) ~= numel(unique(mets));
end

% Create separate table for transport reactions
fseof.transportTargets = fseof.rxnTargets(transportRxns,:);
fseof.rxnTargets(transportRxns,:) = [];

% Report gene targets
geneIdx              = cellfun(@(x) find(strcmpi(model.genes,x)),genes);
fseof.geneTargets(:,1) = genes;
fseof.geneTargets(:,2) = model.geneShortNames(geneIdx);
fseof.geneTargets(:,3) = num2cell(slope_genes);
fseof.geneTargets(:,4) = target_type_genes;
fseof.geneTargets(:,5) = essentiality;

% Save results in a file if defined
if outputFile
    % Write file with gene targets
    writetable(cell2table(fseof.geneTargets, ...
        'VariableNames', {'gene_IDs' 'gene_names' 'slope' 'action' 'essentiality'}), ...
        fullfile(filePath, 'ecFSEOF_genes.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);

    % Write file with rxn targets
    writetable(cell2table(fseof.rxnTargets, ...
        'VariableNames', {'rxn_IDs' 'rxnNames' 'slope' 'grRules' 'rxn_formula'}), ...
        fullfile(filePath, 'ecFSEOF_rxns.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);

    % Write file with transport targets
    writetable(cell2table(fseof.transportTargets, ...
        'VariableNames', {'rxn_IDs' 'rxnNames' 'slope' 'grRules' 'rxn_formula'}), ...
        fullfile(filePath, 'ecFSEOF_transporter.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'QuoteStrings', false);

    disp(['ecFSEOF results stored at: ' newline filePath]);
end

end
