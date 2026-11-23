# Changelog

## 4.0.0b3 (2026-11-23)

Third beta toward 4.0. The main changes are the removal of GECKO 3 functions, a stricter `model_adapter.toml`, scoring of the evotune functions with the `experiments` functions of RAVEN, and functions that compare the kcats of all sources. Changes below are relative to 4.0.0b2.

GECKO 4.0.0b3 requires RAVEN 3.0.0b3, which has the `experiments` functions that GECKO calls, and `GECKOInstaller` checks for them.

A `model_adapter.toml` written for 4.0.0b2 needs manual changes. Remove the keys `okp.method`, `[kcat_sources.okp_<x>] method`, `match_okp`, `[evotune.sigma0_log_source]` and `kcat_aggregate_brenda`, which are errors now. Set `gr_exp`, `[uniprot] id` and `[kegg] id`, which have no default. Run `showAdapterParams` to list the values the adapter resolves to.

### Breaking changes

* `updateProtPool`, `mergeDLKcatAndFuzzyKcats`, `applyCustomKcats` and the adapter template `adapterTemplate.m` moved to `src/deprecated` and warn. `truncateValues` is removed ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `bayesianSensitivityTuning`, its helpers `abc_max`, `getrSample` and `updateprior`, `getcarbonnum` and `findMaxValue` raise `GECKO:removedFunction` and name the replacement ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `runDLKcat` and `writeDLKcatInput` warn and then raise an error that names `kcatsFromDLKcat` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `collectKcats` and `assignKcats` raise `collectKcats:legacyCustomKcats` on a project with `data/customKcats.tsv` only. `convertCustomKcats` converts the file ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* A GECKO 3 adapter class warns `GECKO:deprecatedAdapterClass`. `model_adapter.toml` configures GECKO and geckopy, and `convertAdapterToToml` converts a class ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `[uniprot] id` and `[kegg] id` are empty by default (was: the identifiers of *Saccharomyces cerevisiae*). `loadDatabases` raises `loadDatabases:noUniprotID` when `data/uniprot.tsv` is missing, and an empty KEGG identifier skips KEGG with a warning ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `loadDatabases` and `getECfromDatabase` always split a gene cell of `uniprot.tsv` with several names into one row per name, and the option `splitGeneCells` is removed ([#36](https://github.com/CompMetEng/GECKO/pull/36)).
* `[kcat_sources.okp] methods` is the only setting that names OpenKineticsPredictor methods. `okp.method`, `[kcat_sources.okp_<x>] method`, `match_okp`, `[evotune.sigma0_log_source]` and `kcat_aggregate_brenda` are removed ([#34](https://github.com/CompMetEng/GECKO/pull/34)).
* The default OpenKineticsPredictor method is KinForm-H alone (was: CataPro) ([#34](https://github.com/CompMetEng/GECKO/pull/34)).
* Evotune source groups are `[evotune.source_groups.<name>]` tables with `sources` and `sigma0_log` ([#34](https://github.com/CompMetEng/GECKO/pull/34)).
* A BRENDA call without an aggregation reads `[kcat_sources.brenda] aggregation`, the median for the default preset ([#34](https://github.com/CompMetEng/GECKO/pull/34)).
* `gr_exp` has no default, and unknown or malformed keys in `model_adapter.toml` are errors ([#34](https://github.com/CompMetEng/GECKO/pull/34)).
* A `[kcat_sources.<name>]` table overrides only the keys it sets of the preset (was: it replaced the preset) ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `getKcatAcrossIsozymes` is now `fillKcatsFromIsozymes`. The old name is a stub that warns and forwards ([#34](https://github.com/CompMetEng/GECKO/pull/34)).
* `setProtPoolSize` without `prot_pool_exchange` raises `setProtPoolSize:noPoolExchange` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `sigmaFitter` with a growth rate of 0 or below raises `sigmaFitter:growthRate` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `getConcControlCoeffs` raises an error after a first solve that is not optimal ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `ecFSEOF` without a range to scan raises `ecFSEOF:noRange` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `cmaesKcatTuning` and `screenKcatLeverage` without evotune data raise `evotune:noData` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `flexibilizeEnzConcs` on an infeasible input raises `flexibilizeEnzConcs:infeasible` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `getComplexData` without a taxonomic id warns and returns an empty result ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `parseGECKOargs` is removed and every function uses `parseRAVENargs` of RAVEN, so the error identifiers change from `parseGECKOargs:*` to `parseRAVENargs:*` ([#38](https://github.com/CompMetEng/GECKO/pull/38)).
* `loadFluxData` and `fillCarbonNum` are removed, and `addCarbonNum` and `getcarbonnum` raise `GECKO:removedFunction` ([#38](https://github.com/CompMetEng/GECKO/pull/38)).
* The experimental data are `fluxData.tsv`, `maxGrowth.tsv` and `zeroFlux.tsv`, with the column `weight` ([#38](https://github.com/CompMetEng/GECKO/pull/38)).
* The evotune functions take `experimentData` (was: `evoData`), and `evotuneScore` returns the carbon weights as its fourth output ([#38](https://github.com/CompMetEng/GECKO/pull/38)).
* The scoring maximises the biomass reaction of the adapter, whatever the objective of the model is ([#38](https://github.com/CompMetEng/GECKO/pull/38)).
* `ec.source` has one format. `setKcatForReactions` writes the origin `manual` (was: `setKcatForReactions`), and `sensitivityTuning` and `cmaesKcatTuning` keep the origin of the kcats they change and add the tuning method ([#32](https://github.com/CompMetEng/GECKO/pull/32)).

### Experiments

* `constrainFluxData`, `loadEvotuneData` and `evotuneScore` call the experiments functions of RAVEN ([#38](https://github.com/CompMetEng/GECKO/pull/38)).

### Kcat sources

* New: `kcatEvidence` shows, for selected reactions, the kcat of every source next to the kcat of the ecModel ([#37](https://github.com/CompMetEng/GECKO/pull/37)).
* New: `kcatsForReactions` gets the kcats of every source for a few reactions ([#37](https://github.com/CompMetEng/GECKO/pull/37)).
* `kcatsFromOKP` has the option `jobTag` ([#37](https://github.com/CompMetEng/GECKO/pull/37)).
* New: `parseKcatSource` and `formatKcatSource` read and write the provenance of `ec.source` ([#32](https://github.com/CompMetEng/GECKO/pull/32)).
* `kcatsFromDLKcat` takes `via` (`'okp'`, `'file'`, `'docker'`), so the presets `gecko_median` and `gecko_classic` no longer need Docker ([#33](https://github.com/CompMetEng/GECKO/pull/33)).

### Proteomics and enzyme usage

* `getEnzymeUsage` reports the demand ratio of the usage over the measured concentrations ([#29](https://github.com/CompMetEng/GECKO/pull/29)).
* New: `getIsozymeCoverage` and `zeroUnmeasuredIsozymes` ([#30](https://github.com/CompMetEng/GECKO/pull/30)).
* `loadProtData` has the option `zerosAsMeasured`, and an abundance of 0 is read as not measured by default ([#30](https://github.com/CompMetEng/GECKO/pull/30)).
* New: `relaxEnzConcs` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `constrainEnzConcs` has the option `restrictTo` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `getConcControlCoeffs` has the option `method` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).

### Model adapter

* `showAdapterParams` lists every parameter with its value and origin ([#34](https://github.com/CompMetEng/GECKO/pull/34)).

### Fixes

* `downloadUniProt` keeps unreviewed entries by default, as `uniprot.reviewed = false` ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* `loadAdapterParams` keeps the default of an empty adapter value ([#33](https://github.com/CompMetEng/GECKO/pull/33)).

### Documentation and repository

* The tutorials `full_ecModel` and `light_ecModel` are rewritten for GECKO 4 ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* Each change for GECKO 3 users is also on the [GECKO 3 to GECKO 4](https://gecko-docs.readthedocs.io/en/latest/gecko3-to-gecko4.html) page ([#33](https://github.com/CompMetEng/GECKO/pull/33)).
* GECKO is licensed under the GNU General Public License v3, as RAVEN.
* Help texts, comments and the README describe the code literally ([#31](https://github.com/CompMetEng/GECKO/pull/31)).

## 4.0.0b2 (2026-11-02)

Second beta toward 4.0. The main changes are kcat source files with the preset `brenda_ec_tree`, `model_adapter.toml` as the adapter, and agreement with geckopy. Changes below are relative to 4.0.0b1.

### Breaking changes

* The default kcat preset is `brenda_ec_tree`, which estimates kcats from the EC tree ([#10](https://github.com/CompMetEng/GECKO/pull/10)).
* `calculateFfactor` averages each protein over the samples in which it is measured ([#18](https://github.com/CompMetEng/GECKO/pull/18)).
* `loadProtData` raises an error when the sum of `replPerCond` differs from the number of replicate columns ([#18](https://github.com/CompMetEng/GECKO/pull/18)).
* `getEnzymeBottlenecks` ranks the enzymes by the re-solved change of the objective ([#20](https://github.com/CompMetEng/GECKO/pull/20)).
* GECKO requires RAVEN 3.0.0 and ignores the pre-release suffix of the installed RAVEN ([#24](https://github.com/CompMetEng/GECKO/pull/24)).

### Kcat sources

* Gathered kcats are kept as kcat source files, one per source, and assigned as `model_adapter.toml` sets ([#3](https://github.com/CompMetEng/GECKO/pull/3), [#9](https://github.com/CompMetEng/GECKO/pull/9)).
* The adapter is read from `model_adapter.toml`, the file geckopy reads ([#1](https://github.com/CompMetEng/GECKO/pull/1)).
* BRENDA data come from toolbox-assets, and GECKO 3 projects can be upgraded ([#11](https://github.com/CompMetEng/GECKO/pull/11)).
* OpenKineticsPredictor results of CatRange that are returned as class ranges are read ([#2](https://github.com/CompMetEng/GECKO/pull/2)).
* `matchKcat` divides a BRENDA kcat by the stoichiometric coefficient of the substrate it was measured on (was: the smallest coefficient in the reaction) ([SysBioChalmers/GECKO#484](https://github.com/SysBioChalmers/GECKO/pull/484)).
* `loadBRENDAdata` and `fuzzyKcatMatching` take the name-value argument `kcatAggregation` (`'max'` or `'median'`), and `params.kcatAggregation` of the model adapter sets it for a project. The default is `'max'` ([SysBioChalmers/GECKO#485](https://github.com/SysBioChalmers/GECKO/pull/485)).

### Ec models and complexes

* New: `levelComplexConcs` and `getComplexFamilies` ([#19](https://github.com/CompMetEng/GECKO/pull/19)).
* `flexibilizeEnzConcs` has the option `perComplex` ([#19](https://github.com/CompMetEng/GECKO/pull/19)).
* Saved ecModels record the ecModel format and the GECKO release ([#12](https://github.com/CompMetEng/GECKO/pull/12)).

### Fixes

* `getPfbaEnzymes` with `fractionOfOptimum` below 1 works for a negative optimum ([#25](https://github.com/CompMetEng/GECKO/pull/25)).
* `sigmaFitter`, `getStandardKcat` and `getSubsetEcModel` handle a plateau, subsystem case, empty subsystems and light ecModels ([#15](https://github.com/CompMetEng/GECKO/pull/15)).
* `getECstring` splits EC numbers on semicolons as well as spaces ([#14](https://github.com/CompMetEng/GECKO/pull/14)).
* GECKO gives the same result as geckopy in the scenarios where they differed ([#13](https://github.com/CompMetEng/GECKO/pull/13), [#21](https://github.com/CompMetEng/GECKO/pull/21)).

### Documentation and repository

* `version.txt` and the README version exist on `main` only ([#26](https://github.com/CompMetEng/GECKO/pull/26)).
* Links point to CompMetEng ([#8](https://github.com/CompMetEng/GECKO/pull/8)).
* References to the retired Bayesian kcat tuning are removed ([#16](https://github.com/CompMetEng/GECKO/pull/16)).
* The HTML documentation is regenerated ([#17](https://github.com/CompMetEng/GECKO/pull/17)).
* The README has the new logo ([#22](https://github.com/CompMetEng/GECKO/pull/22)).

## 4.0.0b1 (2026-09-12)

First beta of GECKO 4, which works with RAVEN 3. The main changes are CMA-ES kcat tuning (evotune) in place of Bayesian tuning, forward-direction protein usage reactions, name-value arguments, verb-first function names, and functions ported from geckopy.

### Breaking changes

* The protein usage reactions `usage_prot_<id>` and the exchange `prot_pool_exchange` have the forward direction, as in geckopy (was: the reverse direction). `usage_prot_<id>` is `prot_pool -> prot_<id>` with bounds `[0, 1000]`, and `setProtPoolSize` sets the upper bound of `prot_pool_exchange` ([SysBioChalmers/GECKO#419](https://github.com/SysBioChalmers/GECKO/pull/419)).
* The Bayesian ABC-SMC kcat tuning (`bayesianSensitivityTuning`, `abc_max` and their helpers) and `params.bayesian` are removed. `params.evotune` and the functions in `src/kcat_tuning/evotune` replace them ([SysBioChalmers/GECKO#477](https://github.com/SysBioChalmers/GECKO/pull/477)).
* `src/kcat_sensitivity_analysis` is renamed to `src/kcat_tuning` ([SysBioChalmers/GECKO#477](https://github.com/SysBioChalmers/GECKO/pull/477)).
* Eight functions have verb-first names: `getEnzymeUsage` (was: `enzymeUsage`), `getPfbaEnzymes` (was: `pfbaEnzymes`), `assignKcatValues` (was: `selectKcatValue`), `createGECKOproject` (was: `startGECKOproject`), `buildGECKOdoc` (was: `updateGECKOdoc`), `applyECcodes` (was: `copyECtoGEM`), `fillCarbonNum` (was: `addCarbonNum`) and `updatePrior` (was: `updateprior`). The old names are wrappers that warn `GECKO:deprecatedName` ([SysBioChalmers/GECKO#468](https://github.com/SysBioChalmers/GECKO/pull/468)).
* `src/geckomat` is flattened into `src`. `get_enzyme_data` is split into `download_databases` and `get_enzyme_data`, and `utilities` into `io`, `flux_analysis` and further folders ([SysBioChalmers/GECKO#473](https://github.com/SysBioChalmers/GECKO/pull/473)).
* BRENDA data are `kcat.tsv`, `mw.tsv` and `sa.tsv`, the BRENDA 2026.1 snapshots of geckopy (was: `max_KCAT.txt`, `max_MW.txt` and `max_SA.txt`), and `loadBRENDAdata` reads them ([SysBioChalmers/GECKO#455](https://github.com/SysBioChalmers/GECKO/pull/455)).
* The legacy Python BRENDA scraper in `src/geckopy` is removed ([SysBioChalmers/GECKO#436](https://github.com/SysBioChalmers/GECKO/pull/436)).
* `getStandardKcat` takes the median of the kcats of a subsystem, counting only reactions with a real kcat (was: the mean of all mapped rows, including zero and NaN kcats) ([SysBioChalmers/GECKO#475](https://github.com/SysBioChalmers/GECKO/pull/475)).

### Kcat sources and tuning

* New: `submitOpenKineticsPredictor` and `fetchOpenKineticsPredictor` submit a job to the OpenKineticsPredictor REST API and fetch the result ([SysBioChalmers/GECKO#441](https://github.com/SysBioChalmers/GECKO/pull/441)).
* New: `mergeKcats` merges any number of kcat lists and resolves each reaction by a `sourcePriority` ([SysBioChalmers/GECKO#457](https://github.com/SysBioChalmers/GECKO/pull/457)).
* New: `reviewKcatAssignment` and `ecStatsFromBrenda` flag kcat assignments to curate before a tuning run ([SysBioChalmers/GECKO#477](https://github.com/SysBioChalmers/GECKO/pull/477)).
* New: the evotune functions `loadEvotuneData`, `evotuneScore`, `screenKcatLeverage`, `selectTunableMask` and `cmaesKcatTuning` tune kcats with CMA-ES ([SysBioChalmers/GECKO#477](https://github.com/SysBioChalmers/GECKO/pull/477)).

### Enzyme data and analysis

* New: `downloadKEGG` and `downloadUniProt` are standalone functions, and `loadDatabases` calls them ([SysBioChalmers/GECKO#456](https://github.com/SysBioChalmers/GECKO/pull/456)).
* `makeEcModel` falls back to KEGG for the protein data of a gene that UniProt does not match ([SysBioChalmers/GECKO#458](https://github.com/SysBioChalmers/GECKO/pull/458)).
* New: `getEnzymeBottlenecks` ranks enzymes by the shadow price of their mass balance ([SysBioChalmers/GECKO#438](https://github.com/SysBioChalmers/GECKO/pull/438)).
* New: `pfbaEnzymes` (now `getPfbaEnzymes`) finds the flux distribution with the lowest enzyme usage at the optimum ([SysBioChalmers/GECKO#439](https://github.com/SysBioChalmers/GECKO/pull/439)).
* New: `relaxProteomicsGreedy` relaxes proteomics constraints in the order of the shadow price until a target growth rate is reached ([SysBioChalmers/GECKO#440](https://github.com/SysBioChalmers/GECKO/pull/440)).
* `calculateMW` matches `calculate_mw` of geckopy ([SysBioChalmers/GECKO#459](https://github.com/SysBioChalmers/GECKO/pull/459)).

### Interface

* Optional arguments of the functions in `src` can be given as name-value pairs, and positional calls still work ([SysBioChalmers/GECKO#416](https://github.com/SysBioChalmers/GECKO/pull/416)).
* The help texts of the functions are NumPy-style docstrings ([SysBioChalmers/GECKO#415](https://github.com/SysBioChalmers/GECKO/pull/415)).
* `ModelAdapterManager` warns once per session that GECKO 4 is not fully backward compatible with GECKO 3, with a link to the GECKO 3 to GECKO 4 page ([SysBioChalmers/GECKO#476](https://github.com/SysBioChalmers/GECKO/pull/476), [SysBioChalmers/GECKO#482](https://github.com/SysBioChalmers/GECKO/pull/482)).
* `version.txt` is `4.0.0b1`, and `GECKOInstaller` compares versions correctly, also when the installed version is a beta ([SysBioChalmers/GECKO#482](https://github.com/SysBioChalmers/GECKO/pull/482)).

### Fixes

* `downloadKEGG` trims each entry after splitting on `///`, and retries a download at most five times with exponential backoff ([SysBioChalmers/GECKO#423](https://github.com/SysBioChalmers/GECKO/pull/423)).
* `GECKOInstaller` identifies RAVEN 2.11.0 as the newer version ([SysBioChalmers/GECKO#413](https://github.com/SysBioChalmers/GECKO/pull/413)).
* The tutorial `full_ecModel` runs end to end on GECKO 4 and RAVEN 3 ([SysBioChalmers/GECKO#481](https://github.com/SysBioChalmers/GECKO/pull/481), [SysBioChalmers/GECKO#483](https://github.com/SysBioChalmers/GECKO/pull/483)).
