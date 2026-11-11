% This script builds and analyses a light ecModel with GECKO 4, starting
% from Human-GEM, the genome-scale metabolic model (GEM) of Homo sapiens,
% and makes a context-specific ecModel of the HT-29 cell line. The sections
% follow the GECKO guide
% (https://gecko-docs.readthedocs.io/en/latest/guide/): project and model
% adapter, ecModel structure, kcat sources and assignment, protein pool,
% kcat tuning, sigma fit, saving, proteomics and analysis. The published
% GECKO 3 protocol (https://doi.org/10.1038/s41596-023-00931-7) describes
% GECKO 3.
%
% The ecModel made here demonstrates the functions of GECKO. It has not
% been curated or evaluated for use outside this tutorial.
%
% The comments explain what differs for a light ecModel.
% tutorials/full_ecModel/protocol.m explains each step in more detail.

%% Installation
% Install RAVEN and GECKO, see tutorials/full_ecModel/protocol.m. The
% light ecModel of Human-GEM is large, and the script uses Gurobi as
% solver.
checkRaven; % On RAVEN 2, run checkInstallation.
GECKOInstaller.install
setRavenSolver('gurobi')

%% Project and model adapter
% The project of this tutorial is GECKO/tutorials/light_ecModel.
% model_adapter.toml sets, among others, conv_gem (Human-GEM 1.15.0 in
% models/human-GEM.yml, reduced in size as described in the file),
% c_source (MAR09034, glucose exchange), bio_rxn (MAR13082, biomass), the
% protein pool parameters of the NCI-60 cell lines, and [uniprot] and
% [kegg] (UniProt proteome UP000005640, KEGG organism hsa), which have no
% default and are set for each organism.
adapterLocation = fullfile(findGECKOroot,'tutorials','light_ecModel');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

%% ecModel structure
model = loadConventionalGEM();

% The gene identifiers of Human-GEM are Ensembl identifiers without a
% version suffix, and no UniProt field has this format.
% data/uniprotConversion.tsv maps them to UniProt identifiers: it holds the
% first and the fourth column of genes.tsv of Human-GEM 1.15.0
% (https://github.com/SysBioChalmers/Human-GEM/blob/194ebe5431c83e25f78df61caacad2fa485b5cb4/model/genes.tsv).
% makeEcModel converts the genes with this table when the file exists.
%
% With 'geckoLight' true, makeEcModel makes a light ecModel. Reactions with
% isozymes are not copied: ec.rxns has one row per isozyme, with a number
% prefix (001_, 002_, ...), and each reaction draws from the protein pool
% (prot_pool) with the cost of its isozyme with the smallest MW/kcat. A
% light ecModel has no pseudometabolite or usage reaction per enzyme, so
% it cannot hold proteomics constraints or report the usage of each
% enzyme.
%
% noUniprot lists the genes that neither data/uniprot.tsv nor
% data/kegg.tsv matches; it is empty for this model.
[ecModel, noUniprot] = makeEcModel(model,'geckoLight',true);
fprintf('%d genes without a UniProt or KEGG match.\n', numel(noUniprot))

% Subunit stoichiometry of enzyme complexes, from data/ComplexPortal.json.
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

%% Kcat sources and assignment
% The project uses the default preset brenda_ec_tree, as the full_ecModel
% tutorial: the BRENDA median of measurements of the organism on the
% substrate of the reaction, the EC tree, then the predictions of the
% OpenKineticsPredictor (OKP) method KinForm-H. The project has no
% data/kcats/custom.tsv, so assignKcats reports the custom stage as left
% out.

% EC numbers from Human-GEM, and from UniProt and KEGG for the reactions
% without one.
ecModel = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel = getECfromDatabase(ecModel,'ecRxns',noEC);

% SMILES of the substrates, for the OKP predictions.
[ecModel, noSmiles] = findMetSmiles(ecModel);

% The kcat source files are keyed by the reaction of the GEM and the
% enzyme set. They do not depend on the number prefixes of a light
% ecModel, so the same kind of file serves full and light ecModels. The
% tutorial includes brenda.tsv, ec_tree.tsv and okp_kinform_h.tsv, and
% updateKcatSources reports them as current. In a new project it gathers
% them; tutorials/full_ecModel/protocol.m describes the BRENDA download
% and the OKP API key.
kcatReport = updateKcatSources(ecModel);

% assignKcats selects and applies the kcats. A light ecModel skips the
% step that gives a reaction without a kcat the mean kcat of its isozymes.
% Reactions without a kcat get the standard kcat (getStandardKcat). In a
% light ecModel, a reaction without a gene association draws from the
% protein pool with the cost standardMW/standardKcat; it has no
% pseudo-enzyme prot_standard.
[ecModel, kcatDecisions] = assignKcats(ecModel);
groupcounts(kcatDecisions,'stage')

%% Protein pool
% Ptot, f and sigma of model_adapter.toml. sigma is set to 0.1.
ecModel = setProtPoolSize(ecModel);
sol = solveLP(ecModel,1);
fprintf('Growth rate: %g /hour.\n', sol.f)

%% Kcat tuning and sigma fit
% sensitivityTuning and sigmaFitter work on light ecModels: in a light
% ecModel, sensitivityTuning raises the kcat of the reaction that draws
% most from the protein pool. They are not run here. The ecModel of the
% generic Human-GEM does not represent a cell type, so neither the
% exchange reactions that the medium (for cultured cells) or the tissue
% opens nor the maximum growth rate that the model should reach is
% defined. Tune the kcats and fit sigma on an ecModel of a specific cell
% type with its measured growth rate, as described in
% tutorials/full_ecModel/protocol.m.

%% Saving
% saveEcModel without a file name writes models/ecModel.yml.
saveEcModel(ecModel);

%% Proteomics
% Proteomics constraints need a usage reaction per enzyme, which a light
% ecModel does not have; constrainEnzConcs raises an error on a light
% ecModel. Build a full ecModel for proteomics, see
% tutorials/full_ecModel/protocol.m.

%% Analysis: objectives
% Maximise growth, then fix growth at 99 % of the maximum and minimise the
% protein pool usage. This gives the enzyme allocation with the smallest
% total enzyme mass at that growth rate.
ecModel = loadEcModel();
ecModel = setParam(ecModel,'obj',params.bioRxn,1);
sol = solveLP(ecModel);
fprintf('Maximum growth rate: %g /hour.\n', sol.f)
ecModel = setParam(ecModel,'lb',params.bioRxn,0.99*sol.f);
ecModel = setParam(ecModel,'obj','prot_pool_exchange',-1);
sol = solveLP(ecModel);
fprintf('Minimum protein pool usage: %g mg/gDCW.\n', -sol.f)

% A light ecModel has no usage reaction per enzyme, so getEnzymeUsage
% raises an error. mapRxnsToConv maps the fluxes to the reactions of the
% GEM; enzUsageFlux holds only the flux of prot_pool_exchange.
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel,model,sol.x);

% tutorials/full_ecModel/code/plotlightVSfull.m compares a light and a
% full ecModel of yeast-GEM.

%% Analysis: context-specific ecModel
% A context-specific GEM describes the metabolism of one tissue, cell line
% or strain. models/HT29-GEM.yml is a GEM of the HT-29 cell line, made from
% the same Human-GEM with ftINIT of RAVEN (code/HT29_ftINIT.m).
HT29 = readYAMLmodel(fullfile(params.path,'models','HT29-GEM.yml'));

% getSubsetEcModel trims the generic ecModel to the genes and reactions of
% the context-specific GEM. The context-specific ecModel keeps the kcats
% and the protein pool of the generic ecModel; setProtPoolSize sets the
% protein pool of the context when its total protein content is known.
ecModel = loadEcModel();
ecHT29 = getSubsetEcModel(ecModel,HT29);

% The maximum growth rate of the three models.
sol = solveLP(HT29);
fprintf('Growth rate of HT29-GEM: %.3f /hour.\n', sol.f)
sol = solveLP(ecModel);
fprintf('Growth rate of the ecModel of Human-GEM: %.3f /hour.\n', sol.f)
sol = solveLP(ecHT29);
fprintf('Growth rate of the ecModel of HT29-GEM: %.3f /hour.\n', sol.f)
