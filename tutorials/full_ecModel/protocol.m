% This script builds and analyses a full ecModel with GECKO 4, starting
% from yeast-GEM, the genome-scale metabolic model (GEM) of Saccharomyces
% cerevisiae. The sections follow the GECKO guide
% (https://gecko-docs.readthedocs.io/en/latest/guide/): project and model
% adapter, ecModel structure, kcat sources and assignment, protein pool,
% kcat tuning, sigma fit, saving, proteomics and analysis. The published
% GECKO 3 protocol (https://doi.org/10.1038/s41596-023-00931-7) describes
% GECKO 3.
%
% The ecModel made here demonstrates the functions of GECKO. It has not
% been curated or evaluated for use outside this tutorial. The ecModels of
% S. cerevisiae are distributed by the yeast-GEM repository since release
% 9.1.0 (https://github.com/SysBioChalmers/yeast-GEM).
%
% tutorials/light_ecModel/protocol.m builds a light ecModel of Human-GEM
% and a context-specific ecModel.

%% Installation
% Install RAVEN (https://raven-docs.readthedocs.io/en/latest/installation/raven.html)
% and GECKO (https://gecko-docs.readthedocs.io/en/latest/installation/),
% for example as MATLAB Add-Ons. The script runs with GLPK, the solver
% that RAVEN includes; setRavenSolver('gurobi') selects Gurobi when it is
% installed.
checkRaven; % On RAVEN 2, run checkInstallation.

% GECKOInstaller adds the GECKO folders to the MATLAB path. It is in the
% GECKO root folder.
GECKOInstaller.install

%% Project and model adapter
% A GECKO project is a folder with the subfolders code, data, models and
% output, and a model_adapter.toml with the parameters of the organism and
% the model. createGECKOproject makes a new project with a template
% model_adapter.toml. The project of this tutorial is
% GECKO/tutorials/full_ecModel, so the next line is a comment.
%createGECKOproject('name', 'ecYeastGEM', 'path', pwd);

% model_adapter.toml of this project sets, among others:
% - conv_gem: the starting GEM, yeast-GEM 8.6.2 in models/yeast-GEM.yml;
% - org_name, c_source (r_1714, glucose exchange), bio_rxn (r_4041,
%   biomass) and enzyme_comp (cytoplasm);
% - sigma, p_tot, f and gr_exp for the protein pool and the target growth
%   rate (0.41 /hour);
% - [uniprot] and [kegg]: the UniProt proteome UP000002311 and the KEGG
%   organism sce. These identifiers have no default and are set for each
%   organism;
% - [complex]: the taxonomy identifier for Complex Portal;
% - [kcat_assignment] and [kcat_sources.okp]: how kcats are selected (see
%   the section on kcat sources);
% - [evotune]: the settings of the advanced kcat tuning;
% - [matlab]: anaerobicModel_GECKO in the code folder, which the advanced
%   kcat tuning calls for anaerobic conditions.
% The same file configures geckopy, the Python version of GECKO.
% src/model_adapter/adapterParameters.tsv lists every key and its default.

% setDefault reads model_adapter.toml and makes it the default model
% adapter: every GECKO function that needs a model adapter uses it, until
% another adapter is set or MATLAB is closed. The code folder of the
% project is added to the MATLAB path.
adapterLocation = fullfile(findGECKOroot,'tutorials','full_ecModel');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();
% params holds the values that model_adapter.toml had when it was read. A
% change to params in the workspace has no effect on the default adapter.
% After a change to model_adapter.toml, run ModelAdapterManager.setDefault
% again.

%% ecModel structure
% loadConventionalGEM loads the GEM that conv_gem names. A GEM in another
% location can be loaded with importModel (SBML) or readYAMLmodel (YAML).
model = loadConventionalGEM();

% makeEcModel makes a full ecModel: each isozyme gets its own copy of the
% reaction (suffix _EXP_<n>), and each enzyme gets a pseudometabolite
% prot_<UniProt> and a usage reaction usage_prot_<UniProt> that draws from
% the protein pool. The genes are matched to data/uniprot.tsv and
% data/kegg.tsv, which makeEcModel downloads when they do not exist.
% noUniprot lists the genes without a match; it is empty for yeast-GEM.
[ecModel, noUniprot] = makeEcModel(model,'geckoLight',false);

% applyComplexData sets the subunit stoichiometry of the enzyme complexes
% from data/ComplexPortal.json (Complex Portal,
% https://www.ebi.ac.uk/complexportal/). getComplexData downloads that
% file for the taxonomy identifier in [complex]; the tutorial has the
% file, so the next line is a comment. foundComplex lists the complexes
% that match the gene rule of a reaction, and proposedComplex the
% complexes that match in part.
%complexInfo = getComplexData();
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

%% Kcat sources and assignment
% assignKcats selects one kcat per reaction from the kcat source files in
% data/kcats, in the order that [kcat_assignment] of model_adapter.toml
% defines. This tutorial uses the default preset brenda_ec_tree:
% 1. custom kcats from data/kcats/custom.tsv;
% 2. the BRENDA median of measurements of the organism on the substrate of
%    the reaction (brenda.tsv);
% 3. the EC tree: the median kcat of the EC number, shrunk towards its EC
%    class when the EC number has few measurements (ec_tree.tsv);
% 4. for the remaining reactions, the predictions of the
%    OpenKineticsPredictor (OKP) method KinForm-H ([kcat_sources.okp]
%    methods), in okp_kinform_h.tsv.
% kcatAssignmentConfig shows the stages of a preset, and the guide page
% "Choose the kcat sources" compares the presets.

% BRENDA and the EC tree read the EC numbers of the reactions. The EC
% numbers of yeast-GEM have not been curated, so all EC numbers are taken
% from UniProt and KEGG. For a GEM with curated EC numbers, getECfromGEM
% fills them from the model first, and getECfromDatabase fills only the
% reactions without one:
%   ecModel = getECfromGEM(ecModel);
%   noEC = cellfun(@isempty, ecModel.ec.eccodes);
%   ecModel = getECfromDatabase(ecModel,'ecRxns',noEC);
ecModel = getECfromDatabase(ecModel);

% The OKP methods read the SMILES of the substrates. findMetSmiles reads
% them from data/smilesDB.tsv and queries PubChem for the metabolite names
% that are not in that file. noSmiles lists the metabolites without a
% SMILES.
[ecModel, noSmiles] = findMetSmiles(ecModel);

% Custom kcats. data/kcats/custom.tsv holds the curated kcats of this
% project, which every preset uses before all other sources. Most rows
% were curated during the development of the yeast-GEM ecModels with
% GECKO 1 and 2, and kept in data/customKcats.tsv by GECKO 3.
% convertCustomKcats(ecModel) wrote that file as data/kcats/custom.tsv,
% keyed by the reaction of the GEM and the enzyme set, so the file also
% applies to light ecModels.
%
% The row for r_0079 (5'-phosphoribosylformyl glycinamidine synthetase,
% enzyme P38972) was added by hand. Without it, the kcat of r_0079 came
% from a BRENDA kcat of Escherichia coli with NH4+ as the nitrogen donor
% (EC 6.3.5.3, https://www.brenda-enzymes.org/enzyme.php?ecno=6.3.5.3).
% The reaction in the model uses glutamine. The paper that reports the
% value (Schendel et al. 1989, https://doi.org/10.1021/bi00432a017) states
% that "NH3 can replace glutamine as a nitrogen donor with a Km = 1 M and a
% turnover of 3 min-1 (2% glutamine turnover)", and it reports a specific
% activity of 2.15 umol/min/mg protein, which BRENDA lists as well. A
% specific activity times the molecular weight (MW) of the enzyme gives a
% kcat:
%   sa = 2.15;                                                % umol/min/mg = mmol/min/g
%   mw = ecModel.ec.mw(strcmp(ecModel.ec.enzymes,'P38972'));  % 148905 g/mol
%   kcat = sa / 1000 / 60 * mw                                % 5.3357625 /s
% custom.tsv holds this value in the row
%   r_0079  P38972  5.3357625  custom
% A kcat set with setKcatForReactions is replaced the next time
% assignKcats runs, so a curated kcat is kept in custom.tsv.

% updateKcatSources lists the source files that the preset reads and
% gathers the files that are missing or outdated. custom.tsv is never
% gathered or changed. The tutorial includes brenda.tsv, ec_tree.tsv and
% okp_kinform_h.tsv, so updateKcatSources reports them as current and
% gathers nothing.
%
% In a new project, updateKcatSources makes brenda.tsv and ec_tree.tsv
% from the BRENDA files, which it downloads on first use. For
% okp_kinform_h.tsv it submits a job to OpenKineticsPredictor, a web
% service that needs an API key: set the environment variable OKP_API_KEY,
% or write the key in data/okpApiKey.txt (generate a key at
% https://predictor.openkinetics.org/, and do not commit the file to a
% public repository). The first call reports the source as pending; call
% updateKcatSources again until no source is pending.
kcatReport = updateKcatSources(ecModel);

% assignKcats selects the kcats, writes ec.kcat and ec.source, and applies
% the kcats to the stoichiometry (applyKcatConstraints). After the stages,
% a reaction without a kcat gets the mean kcat of its isozymes, and then
% the standard kcat (getStandardKcat). Reactions without a gene
% association get the pseudo-enzyme prot_standard with the standard kcat,
% except exchange, transport and pseudo-reactions and the reactions in
% data/pseudoRxns.tsv. The decision per reaction (stage, source, kcat and the
% number of values) is written to output/<model id>_kcatDecisions.tsv.
[ecModel, kcatDecisions, assignReport] = assignKcats(ecModel);
groupcounts(kcatDecisions,'stage')

% assignReport lists, per source file, the rows that were not applied:
% invalid (no valid kcat), unmatched (no reaction of the model) and
% partial (part of an enzyme set).
assignReport.sources.custom

% collectKcats determines the same kcats without changing the model, to
% compare presets before one is applied.

%% Protein pool
% setProtPoolSize sets the upper bound of prot_pool_exchange to
% Ptot * f * sigma * 1000 mg/gDCW. Ptot (total protein content), f (mass
% fraction of the protein that the enzymes of the model account for) and
% sigma (average saturation of the enzymes) are read from
% model_adapter.toml unless they are given. calculateFfactor calculates f
% from the protein abundances of data/paxDB.tsv (PAXdb,
% https://pax-db.org/).
f = calculateFfactor(ecModel);
ecModel = setProtPoolSize(ecModel,'f',f);

%% Kcat tuning
% With the uptake of glucose open and growth as the objective, the growth
% rate is compared with gr_exp (0.41 /hour), the growth rate of
% S. cerevisiae when the nutrient uptake does not limit growth.
ecModel = setParam(ecModel,'lb',params.c_source,-1000);
ecModel = setParam(ecModel,'obj',params.bioRxn,1);
sol = solveLP(ecModel,1);
fprintf('Growth rate before kcat tuning: %f /hour.\n', sol.f)

% The ecModel before tuning is kept for a comparison in the analysis.
% Files whose name contains _stage are not tracked by git.
saveEcModel(ecModel,'ecYeastGEM_stage_untuned.yml');

% Standard kcat tuning. sensitivityTuning selects the enzyme with the
% highest usage, multiplies the kcat of the reaction that uses most of it
% by 10, and repeats this until the model reaches gr_exp. When the model
% already reaches gr_exp, it changes no kcat. tunedKcats lists each changed
% reaction with its kcat before and after tuning, and the source of the
% kcat before tuning. ec.notes of a tuned reaction records the kcat and
% the source before tuning.
[ecModel, tunedKcats] = sensitivityTuning(ecModel);
struct2table(tunedKcats)

% A tuning raises the kcats that limit growth; it does not check whether
% the value before tuning was wrong. For each tuned kcat, look up where the
% value before tuning came from and search the literature for a measured
% value. loadKcatSource returns the rows of a source file in data/kcats for
% the reactions of the model. The detail column of brenda.tsv holds the EC
% wildcard level and the origin of the match, and that of ec_tree.tsv the
% number of measurements and the depth in the EC tree. A measured value
% that supports a tuned kcat is added to custom.tsv, as for r_0079 above.
for i = 1:numel(tunedKcats.rxns)
    fprintf('%s (%s): %g -> %g /s, source before tuning: %s\n', tunedKcats.rxns{i}, ...
        tunedKcats.rxnNames{i}, tunedKcats.oldKcat(i), tunedKcats.newKcat(i), tunedKcats.source{i});
    sourceName = strtok(tunedKcats.source{i});
    if isfile(fullfile(params.path,'data','kcats',[sourceName '.tsv']))
        kcatList = loadKcatSource(ecModel,sourceName);
        idx = strcmp(kcatList.rxns, tunedKcats.rxns{i});
        disp(table(kcatList.kcats(idx), kcatList.eccodes(idx), kcatList.detail(idx), ...
            'VariableNames', {'kcat','eccodes','detail'}))
    end
end

% Advanced kcat tuning. cmaesKcatTuning fits kcats to measured maximum
% growth rates (data/evotuneMaxGrowth.tsv) and exchange fluxes
% (data/evotuneFluxData.tsv, with data/evotuneZeroExch.tsv), measured on
% several carbon sources and in aerobic and anaerobic conditions. It scores
% each set of kcats by the carbon-weighted root-mean-square error (RMSE)
% against these data and searches with CMA-ES, an evolutionary optimiser.
% screenKcatLeverage first ranks the kcats by their effect on the fit, and
% the search varies the top kcats of the screen. The trust in each kcat
% source ([evotune.source_groups] and [evotune.sigma0_log_source] of
% model_adapter.toml) limits how far each kcat moves, and
% prior_penalty_weight penalises large changes. The settings depend on the
% model, on the certainty of its kcats and on the data, and the search runs
% many FBA problems, so it is shown here and not run:
%   screen = screenKcatLeverage(ecModel, 'nProc', 8);
%   [ecModel, result] = cmaesKcatTuning(ecModel, 'screen', screen, 'nProc', 8, 'seed', 0);
% The standard and the advanced tuning are alternatives: cmaesKcatTuning
% starts from the kcats before sensitivityTuning. The guide page "Tune
% kcats against experimental data" describes the data files and the
% settings, and the yeast-GEM repository distributes ecModels tuned this
% way. The rest of this script uses the result of sensitivityTuning.

%% Sigma fit
% sensitivityTuning multiplies each kcat by 10, so the growth rate after
% tuning can be above gr_exp. sigmaFitter sets the protein pool for sigma
% = 0.01, 0.02, ..., 1 and selects the sigma at which the growth rate is
% closest to gr_exp. f is given, because sigmaFitter otherwise reads f from
% model_adapter.toml. The returned ecModel has the protein pool of the
% fitted sigma, and the bounds of the biomass reaction are 0 and 1000.
[ecModel, sigma] = sigmaFitter(ecModel,'f',f,'makePlot',false);
sol = solveLP(ecModel,1);
fprintf('Fitted sigma: %.2f; growth rate: %f /hour.\n', sigma, sol.f)

%% Saving
% saveEcModel writes the ecModel to the models folder of the project, in
% YAML format, which GECKO and geckopy both read. GECKO includes this
% ecModel as tutorials/full_ecModel/models/ecYeastGEM.yml. A file name
% ending in .xml writes SBML, without the ec fields.
saveEcModel(ecModel,'ecYeastGEM.yml');

%% Proteomics
% Absolute proteomics data bound the usage reaction of each measured
% enzyme. data/proteomics.tsv has one condition with three replicates, in
% mg/gDCW. loadProtData filters the proteins and averages the replicates;
% fillEnzConcs writes the concentrations to ec.concs, and
% constrainEnzConcs sets them as upper bounds of the usage reactions.
protData = loadProtData(3);
ecModel = fillEnzConcs(ecModel,protData);
ecModel = constrainEnzConcs(ecModel);

% All usage reactions draw from the protein pool, also those with a
% measured bound. The protein pool of the condition uses the Ptot of the
% condition from data/fluxData.tsv and the f of the proteomics data.
fluxData = loadFluxData();
f = calculateFfactor(ecModel,'protData',protData);
ecModel = setProtPoolSize(ecModel,'Ptot',fluxData.Ptot(1),'f',f,'sigma',sigma);

% constrainFluxData sets the growth rate and the exchange fluxes of the
% first condition of data/fluxData.tsv. 'max' sets the measured growth rate
% as the upper bound; 'loose' bounds each exchange reaction between 0 and
% the measured value.
ecModel = constrainFluxData(ecModel,'fluxData',fluxData,'condition',1,'maxMinGrowth','max','looseStrictFlux','loose');
sol = solveLP(ecModel);
fprintf('Growth rate with proteomics: %f /hour; measured: %f /hour.\n', sol.f, fluxData.grRate(1))

% When the model does not reach the measured growth rate, the enzyme
% concentrations, the kcats, the protein pool or the exchange fluxes limit
% growth. flexibilizeEnzConcs raises the bound of the enzyme that limits
% growth most by foldChange, and repeats this until the model reaches the
% measured growth rate. ec.concs is not changed. relaxProteomicsGreedy and
% relaxEnzConcs are other methods (guide page "Relax proteomics
% constraints").
[ecModel, flexEnz] = flexibilizeEnzConcs(ecModel,'expGrowth',fluxData.grRate(1),'foldChange',10);
% flexEnz lists the relaxed enzymes, ordered by the fold change of their
% bound, from high to low.
struct2table(flexEnz)
sol = solveLP(ecModel);
fprintf('Growth rate after flexibilizeEnzConcs: %f /hour.\n', sol.f)

% The conventional GEM with the same exchange bounds shows whether the
% network itself reaches the measured growth rate.
modelFlux = constrainFluxData(model,'fluxData',fluxData,'condition',1,'maxMinGrowth','max','looseStrictFlux','loose');
sol = solveLP(modelFlux);
fprintf('Growth rate of the conventional GEM: %f /hour.\n', sol.f)

saveEcModel(ecModel,'ecYeastGEM_stage_proteomics.yml');

%% Analysis: Crabtree effect
% The analyses start from the saved ecModel without proteomics.
ecModel = loadEcModel('ecYeastGEM.yml');

% plotCrabtree, in the code folder of the project, minimises glucose
% uptake and then the protein pool usage at growth rates from 0 to 0.4
% /hour. The left graph shows the exchange fluxes of the simulations
% (lines) and of chemostat experiments (circles, data/vanHoek1998.tsv,
% https://doi.org/10.1128/AEM.64.11.4226-4233.1998), with uptake fluxes as
% positive numbers. The right graph shows the fraction of the protein pool
% that the enzymes use. When the protein pool is fully used, the model
% shifts from respiration to fermentation: oxygen uptake falls and ethanol
% excretion rises (the Crabtree effect).
[fluxes, gRate] = plotCrabtree(ecModel);
saveas(gcf,fullfile(params.path,'output','crabtree.pdf'))

% With an unlimited protein pool, the ecModel behaves as a conventional
% GEM: the protein pool never limits growth, and the model does not shift
% to fermentation.
ecModel_infProt = setProtPoolSize(ecModel,'Ptot',Inf);
plotCrabtree(ecModel_infProt);
saveas(gcf,fullfile(params.path,'output','crabtree_infProt.pdf'))

% The ecModel before kcat tuning and sigma fit, for comparison.
ecModel_untuned = loadEcModel('ecYeastGEM_stage_untuned.yml');
plotCrabtree(ecModel_untuned);
saveas(gcf,fullfile(params.path,'output','crabtree_untuned.pdf'))

%% Analysis: objectives and enzyme usage
% Maximise growth, then fix growth at 99 % of the maximum and minimise the
% protein pool usage. setParam with 'obj' and -1 makes solveLP minimise
% the flux of that reaction.
ecModel = setParam(ecModel,'obj',params.bioRxn,1);
sol = solveLP(ecModel);
fprintf('Maximum growth rate: %f /hour.\n', sol.f)
ecModel = setParam(ecModel,'lb',params.bioRxn,0.99*sol.f);
ecModel = setParam(ecModel,'obj','prot_pool_exchange',-1);
sol = solveLP(ecModel);
fprintf('Minimum protein pool usage: %.2f mg/gDCW.\n', -sol.f)

% Enzyme usage at a growth rate of 0.25 /hour with minimal glucose uptake.
% getEnzymeUsage gives the usage of each enzyme, and reportEnzymeUsage
% ranks the enzymes.
ecModel = setParam(ecModel,'obj',params.c_source,1);
ecModel = setParam(ecModel,'lb',params.bioRxn,0.25);
sol = solveLP(ecModel,1);
usageData = getEnzymeUsage(ecModel,sol.x);
usageReport = reportEnzymeUsage(ecModel,usageData);
usageReport.topAbsUsage

% mapRxnsToConv maps the fluxes of the ecModel to the reactions of the
% conventional GEM: the fluxes of the isozyme copies are added, and the
% flux of a _REV reaction is subtracted from the flux of its reaction.
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel,model,sol.x);
fprintf('%d mapped fluxes for %d reactions of the GEM.\n', numel(mappedFlux), numel(model.rxns))

%% Analysis: ecFVA
% Flux variability analysis of the conventional GEM, the ecModel and the
% ecModel with proteomics, with the exchange fluxes of the first condition
% of data/fluxData.tsv. The growth rate is set to the maximum that the
% ecModel with proteomics reaches.
model = loadConventionalGEM();
ecModel = loadEcModel('ecYeastGEM.yml');
ecModelProt = loadEcModel('ecYeastGEM_stage_proteomics.yml');
ecModelProt = setParam(ecModelProt,'obj',params.bioRxn,1);
sol = solveLP(ecModelProt);
fluxData.grRate(1) = floor(sol.f*1e4)/1e4;
fprintf('Growth rate for ecFVA: %f /hour.\n', fluxData.grRate(1))

model = constrainFluxData(model,'fluxData',fluxData,'condition',1,'maxMinGrowth','max','looseStrictFlux','loose');
ecModel = constrainFluxData(ecModel,'fluxData',fluxData,'condition',1,'maxMinGrowth','max','looseStrictFlux','loose');
ecModelProt = constrainFluxData(ecModelProt,'fluxData',fluxData,'condition',1,'maxMinGrowth','max','looseStrictFlux','loose');

% ecFVA returns the flux ranges mapped to the reactions of the GEM.
minFlux = zeros(numel(model.rxns),3);
maxFlux = minFlux;
[minFlux(:,1), maxFlux(:,1)] = ecFVA(model,model);
[minFlux(:,2), maxFlux(:,2)] = ecFVA(ecModel,model);
[minFlux(:,3), maxFlux(:,3)] = ecFVA(ecModelProt,model);

output = [model.rxns, model.rxnNames, num2cell([minFlux(:,1), maxFlux(:,1), ...
    minFlux(:,2), maxFlux(:,2), minFlux(:,3), maxFlux(:,3)])]';
fID = fopen(fullfile(params.path,'output','ecFVA.tsv'),'w');
fprintf(fID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','rxnIDs','rxnNames','minFlux', ...
    'maxFlux','ec-minFlux','ec-maxFlux','ecP-minFlux','ecP-maxFlux');
fprintf(fID,'%s\t%s\t%g\t%g\t%g\t%g\t%g\t%g\n',output{:});
fclose(fID);

% plotEcFVA plots the cumulative distribution of the flux ranges of the
% three models.
plotEcFVA(minFlux,maxFlux);
saveas(gca,fullfile(params.path,'output','ecFVA.pdf'))

%% Analysis: light and full ecModel
% plotlightVSfull, in the code folder, builds a light and a full ecModel
% of yeast-GEM with the same BRENDA kcats and no tuning, maximises growth
% in both, and plots the fluxes mapped to the reactions of the GEM. It
% also prints the time of each step for both structures.
[fluxLight, fluxFull] = plotlightVSfull();

% Reactions whose flux differs by more than 0.1 % between the two
% structures:
fluxRatio = fluxFull ./ fluxLight;
changedFlux = abs(fluxRatio-1) > 0.001;
model.rxnNames(changedFlux)

% tutorials/light_ecModel/protocol.m makes a context-specific ecModel.
