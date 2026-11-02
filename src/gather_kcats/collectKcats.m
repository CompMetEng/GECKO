function [decisions, report, ecModel] = collectKcats(model, varargin)
% collectKcats  Decide one kcat per reaction from the kcat source files.
%
% Combines the kcat source files in data/kcats (see kcatsFromBrenda,
% kcatsFromDLKcat, kcatsFromOKP and saveKcatSource) as [kcat_assignment] in
% model_adapter.toml sets, without changing the model given. assignKcats
% applies the decisions.
%
% The order is a preset or a list of stages:
%
%     [kcat_assignment]
%     preset = "brenda_ec_tree"    # the default; or ec_tree_homologs, ec_tree_okp, okp,
%                                  # gecko_median, gecko_classic
%     isozymes = true              # fill from isozymes after the stages
%     standard = true              # standard kcat for the rest
%
%     # instead of a preset:
%     [[kcat_assignment.stage]]
%     action = "select"            # or "combine"
%     sources = ["custom", "database_top", "ec_tree"]
%     criteria = "max"             # select: among several values of one source
%     top_origin_limit = 1
%
%     [kcat_sources.okp]
%     methods = ["CataPro", "UniKP"]   # the okp source group
%
% A stage decides the reactions that are still open:
% - select: the first listed source with a value (see mergeKcats). Sources
%   are the names of source files (custom, ec_tree, dlkcat, okp_catapro,
%   ...) and the BRENDA quality tiers database_exact, database_top and
%   database_bottom, which read brenda.tsv.
% - combine: all listed sources with a value, each corrected for its bias
%   and weighted by 1 / its error variance, both for the source's evidence
%   for that reaction (the EC tree: measurements of the EC and its EC
%   level; an OKP method: the enzyme's similarity to the method's training
%   sequences). The parameters come with the BRENDA files
%   (kcatCombineParams.tsv). model.ec.source reads 'combined (sd=<log10
%   sd>, n=<sources>)'; sd is the uncertainty for a new measurement.
% In both, okp stands for the OKP methods of [kcat_sources.okp] methods
% (default: all methods with combination parameters).
%
% The presets:
% - brenda_ec_tree (the default): custom kcats, the BRENDA median of
%   measurements of the organism on the reaction's substrate (origin 1,
%   without EC wildcards), the EC tree (kcatsFromECTree); then the OKP
%   methods combined for the reactions left (without an EC number).
% - ec_tree_homologs: custom kcats, the organism's own BRENDA value; then
%   the EC tree and homologs (kcatsFromHomologs) combined; then the OKP
%   methods combined.
% - ec_tree_okp: custom kcats, the organism's own BRENDA value; then the EC
%   tree, homologs and the OKP methods combined; then the OKP methods
%   combined.
% - okp: custom kcats; then the OKP methods combined.
% - gecko_median: custom kcats, the BRENDA median of matches without EC
%   wildcards from the organism or other organisms, with or without the
%   substrate (origin 1 to 4), DLKcat, then the other BRENDA matches (EC
%   wildcards, and specific activities, origin 5 and 6).
% - gecko_classic: as gecko_median with the BRENDA maximum and the
%   specific-activity matches above DLKcat (origin up to 6), as the GECKO 3
%   tutorial.
% Against measured kcats hidden from the search, the EC tree was closer
% than the BRENDA median at origins 2 and 4 and with EC wildcards and as
% close at origins 1 and 3; the EC tree combined with the OKP methods was
% slightly closer than the EC tree, and the OKP methods combined as close.
% In ecModels, values combined with predictions set high-turnover enzymes
% such as ATP synthase too low, which is why the default uses predictions
% only where there is no value from measurements. The BRENDA aggregation is
% [kcat_sources.brenda] aggregation (see updateKcatSources).
%
% After the stages:
% - custom kcats of data/customKcats.tsv are applied (applyCustomKcats),
%   when a stage lists custom but there is no data/kcats/custom.tsv;
% - with isozymes, reactions without a kcat take one from their isozymes
%   (getKcatAcrossIsozymes);
% - with standard, the remaining reactions get the standard kcat
%   (getStandardKcat).
%
% A source file a stage names that does not exist is reported and left out.
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
% writeReport : logical
%     write the decisions to output/<model id>_kcatDecisions.tsv (default
%     true).
%
% Returns
% -------
% decisions : table
%     one row per model.ec reaction: ecRxn, stage (the stage number, or
%     custom, isozymes, standard, or empty when no kcat was found), source
%     (as in model.ec.source), kcat, candidates (the number of values the
%     deciding source had for the reaction, or of sources combined),
%     spreadLog10 (their range, in log10 units) and sdLog10 (the
%     uncertainty of a combined value).
% report : struct
%     - preset: the preset used, or '' with stages;
%     - missing: source files the stages name that do not exist;
%     - sources: per source file loaded, its unmatched, partial and invalid
%       rows (see loadKcatSource).
% ecModel : struct
%     the model with the decided kcats, before applyKcatConstraints; this
%     is what assignKcats returns, after applying the constraints.
%
% Examples
% --------
%     decisions = collectKcats(ecModel);
%     decisions(strcmp(decisions.stage, ''), :)   % reactions without a kcat
%
% See also
% --------
% assignKcats, updateKcatSources, loadKcatSource, mergeKcats, assignKcatValues

checkOptions(varargin, {'modelAdapter','writeReport'}, 'collectKcats');
p = parseGECKOargs(varargin, {'modelAdapter', []; 'writeReport', true});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
config = kcatAssignmentConfig(params);
kcatDir = fullfile(params.path, 'data', 'kcats');

ecModel = model;
nEc = numel(ecModel.ec.rxns);
ecModel.ec.kcat(:) = 0;
ecModel.ec.source(:) = {''};
stage = repmat({''}, nEc, 1);
candidates = zeros(nEc, 1);
spread = nan(nEc, 1);
sd = nan(nEc, 1);
combineParams = [];

report.preset = config.preset;
report.missing = {};
report.sources = struct();
loaded = struct();
legacyCustom = false;

    function list = loadSource(name)
        % The kcat list of a source file, or [] when it does not exist.
        list = [];
        file = fullfile(kcatDir, [name '.tsv']);
        if isfile(file)
            [list, rep] = loadKcatSource(ecModel, file);
            report.sources.(name) = rep;
        elseif strcmp(name, 'custom') && isfile(fullfile(params.path, 'data', 'customKcats.tsv'))
            legacyCustom = true;
        else
            report.missing{end+1} = name;
            printOrange(sprintf(['No kcat source file %s; it is left out. ' ...
                'updateKcatSources gathers it.\n'], file));
        end
    end

for s = 1:numel(config.stages)
    st = config.stages(s);
    if strcmp(st.action, 'combine')
        names = {};
        lists = {};
        for t = 1:numel(st.files)
            name = st.files{t};
            if ~isfield(loaded, name)
                loaded.(name) = loadSource(name);
            end
            if ~isempty(loaded.(name))
                names{end+1} = name; %#ok<AGROW>
                lists{end+1} = loaded.(name); %#ok<AGROW>
            end
        end
        if isempty(names)
            continue
        end
        if isempty(combineParams)
            combineParams = loadCombineParams(modelAdapter);
        end
        openRxns = ecModel.ec.rxns(ecModel.ec.kcat == 0);
        c = combineKcats(names, lists, openRxns, combineParams);
        if isempty(c.rxns)
            continue
        end
        [~, idx] = ismember(c.rxns, ecModel.ec.rxns);
        ecModel.ec.kcat(idx) = c.kcats;
        ecModel.ec.source(idx) = cellstr(compose('combined (sd=%.2f, n=%d)', c.sd, c.n));
        stage(idx) = {sprintf('%d', s)};
        candidates(idx) = c.n;
        spread(idx) = c.spread;
        sd(idx) = c.sd;
        continue
    end
    lists = {};
    usedFiles = {};
    priority = {};
    for t = 1:numel(st.files)
        token = st.files{t};
        name = sourceFile(token);
        if ~isfield(loaded, name)
            loaded.(name) = loadSource(name);
        end
        list = loaded.(name);
        if isempty(list)
            continue
        end
        if ~ismember(name, usedFiles)
            lists{end+1} = list; %#ok<AGROW>
            usedFiles{end+1} = name; %#ok<AGROW>
        end
        if startsWith(token, 'database_')
            priority{end+1} = token; %#ok<AGROW>
        else
            % The labels of the file's rows, at the file's place in the order.
            labels = unique(cellfun(@normalizeSource, list.kcatSource, 'UniformOutput', false), 'stable');
            priority = [priority, labels(~ismember(labels, {'brenda','sabio_rk'}))']; %#ok<AGROW>
        end
    end
    if isempty(lists) || isempty(priority)
        continue
    end
    merged = mergeKcats(lists, priority, 'topOriginLimit', st.topOriginLimit, ...
        'bottomOriginLimit', st.bottomOriginLimit, 'wildcardLimit', st.wildcardLimit);
    if isempty(merged.rxns)
        continue
    end
    open = ecModel.ec.kcat == 0;
    [ecModel, idx] = assignKcatValues(ecModel, merged, 'criteria', st.criteria, 'overwrite', 'false');
    idx = idx(open(idx));
    stage(idx) = {sprintf('%d', s)};
    [~, where] = ismember(merged.rxns, ecModel.ec.rxns);
    valid = merged.kcats > 0;
    for i = idx(:)'
        values = merged.kcats(where == i & valid);
        candidates(i) = numel(values);
        spread(i) = log10(max(values)) - log10(min(values));
    end
end

if legacyCustom
    before = ecModel.ec.kcat;
    ecModel = applyCustomKcats(ecModel, 'modelAdapter', modelAdapter);
    changed = ecModel.ec.kcat ~= before;
    stage(changed) = {'custom'};
    candidates(changed) = 1;
    spread(changed) = 0;
end
if config.isozymes && ~ecModel.ec.geckoLight
    before = ecModel.ec.kcat == 0;
    ecModel = getKcatAcrossIsozymes(ecModel);
    stage(before & ecModel.ec.kcat > 0) = {'isozymes'};
end
if config.standard
    ecModel = getStandardKcat(ecModel, 'modelAdapter', modelAdapter);
    n = numel(ecModel.ec.rxns);
    stage(end+1:n) = {''};
    candidates(end+1:n) = 0;
    spread(end+1:n) = NaN;
    sd(end+1:n) = NaN;
    stage(strcmp(ecModel.ec.source, 'standard') & ecModel.ec.kcat > 0) = {'standard'};
end

decisions = table(ecModel.ec.rxns(:), stage(:), ecModel.ec.source(:), ecModel.ec.kcat(:), ...
    candidates(:), spread(:), sd(:), 'VariableNames', {'ecRxn','stage','source','kcat','candidates', ...
    'spreadLog10','sdLog10'});
if p.writeReport
    outDir = fullfile(params.path, 'output');
    if ~isfolder(outDir)
        mkdir(outDir);
    end
    id = 'ecModel';
    if isfield(model, 'id') && ~isempty(model.id)
        id = char(model.id);
    end
    writetable(decisions, fullfile(outDir, [id '_kcatDecisions.tsv']), 'FileType', 'text', ...
        'Delimiter', '\t');
end
nOpen = sum(strcmp(stage, ''));
fprintf('Kcats decided for %d of %d reactions.\n', numel(stage) - nOpen, numel(stage));
end

function s = normalizeSource(label)
% As mergeKcats folds source labels: lower case, other characters as '_',
% without a leading okp_ (the tag of OpenKineticsPredictor rows).
s = lower(strtrim(char(label)));
s = regexprep(s, '[^0-9a-z]+', '_');
s = regexprep(s, '^_+|_+$', '');
s = regexprep(s, '^okp_', '');
end

function name = sourceFile(token)
% The source file a stage token reads: the BRENDA tiers read brenda.tsv.
if startsWith(token, 'database_')
    name = 'brenda';
else
    name = token;
end
end

function checkOptions(args, names, functionName)
% Options are name-value pairs with known names, so a misspelled option is
% an error instead of being taken as a positional argument.
id = [functionName ':arguments'];
given = args(1:2:end);
if mod(numel(args), 2) ~= 0 || ~all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), given))
    error(id, 'Options of %s are name-value pairs. Options: %s.', functionName, strjoin(names, ', '))
end
unknown = setdiff(cellfun(@char, given, 'UniformOutput', false), names);
if ~isempty(unknown)
    error(id, 'Unknown option of %s: %s. Options: %s.', functionName, strjoin(unknown, ', '), ...
        strjoin(names, ', '))
end
end

function params = loadCombineParams(modelAdapter)
% kcatCombineParams.tsv of the BRENDA folder, or of GECKO's BRENDA release.
file = fullfile(modelAdapter.getBrendaDBFolder(), 'kcatCombineParams.tsv');
if ~isfile(file)
    file = fullfile(fetchBrendaAsset({'kcatCombineParams.tsv'}), 'kcatCombineParams.tsv');
end
lines = strsplit(fileread(file), {'\r\n', '\n'});
params.measurementSd = 0.48;
params.minVariance = 0.01;
params.scaleKeys = {};
params.scales = [];
k = 1;
while k <= numel(lines) && startsWith(lines{k}, '#')
    tok = regexp(lines{k}, '^#\s*([^:]+):\s*(.*)$', 'tokens', 'once');
    if ~isempty(tok)
        key = strtrim(tok{1});
        value = str2double(strtrim(tok{2}));
        if strcmp(key, 'measurement_sd')
            params.measurementSd = value;
        elseif strcmp(key, 'min_variance')
            params.minVariance = value;
        elseif startsWith(key, 'interval_scale ')
            params.scaleKeys{end+1} = key(numel('interval_scale ')+1:end);
            params.scales(end+1) = value;
        end
    end
    k = k + 1;
end
header = strsplit(lines{k}, '\t');
rows = lines(k+1:end);
rows = rows(~cellfun(@isempty, strtrim(rows)));
fields = cellfun(@(r) strsplit(r, '\t', 'CollapseDelimiters', false), rows, 'UniformOutput', false);
fields = vertcat(fields{:});
col = @(name) fields(:, strcmp(header, name));
params.source = col('source');
params.bin = col('bin');
params.bias = str2double(col('bias'));
params.variance = str2double(col('variance'));
end

function c = combineKcats(names, lists, openRxns, params)
% One kcat per reaction from several sources: each source's value corrected
% for its bias and weighted by 1 / its error variance, both for the
% source's evidence bin (see evidenceBin). Several rows of one source for a
% reaction give their median, with the bin of best evidence; a value
% without a known bin gets the source's largest variance. geckopy's
% combine_kcats does the same.
unknown = setdiff(names, params.source);
if ~isempty(unknown)
    error('collectKcats:combine', 'No combination parameters for %s; known sources: %s.', ...
        strjoin(unknown, ', '), strjoin(unique(params.source), ', '))
end
R = {}; S = {}; X = []; W = []; V = [];
for k = 1:numel(names)
    name = names{k};
    list = lists{k};
    keep = ismember(list.rxns, openRxns) & list.kcats > 0 & isfinite(list.kcats);
    rxns = list.rxns(keep);
    values = log10(list.kcats(keep));
    details = list.detail(keep);
    mine = strcmp(params.source, name);
    srcBins = params.bin(mine);
    srcBias = params.bias(mine);
    srcVar = params.variance(mine);
    [~, worst] = max(srcVar);
    [u, ~, g] = unique(rxns, 'stable');
    for r = 1:numel(u)
        rows = find(g == r);
        bins = cellfun(@(d) evidenceBin(name, d), details(rows), 'UniformOutput', false);
        [found, where] = ismember(bins, srcBins);
        if any(found)
            [~, best] = min(srcVar(where(found)));
            pick = where(found);
            pick = pick(best);
        else
            pick = worst;
        end
        R{end+1, 1} = u{r}; %#ok<AGROW>
        S{end+1, 1} = name; %#ok<AGROW>
        V(end+1, 1) = median(values(rows)); %#ok<AGROW>
        X(end+1, 1) = median(values(rows)) - srcBias(pick); %#ok<AGROW>
        W(end+1, 1) = 1 / max(srcVar(pick), params.minVariance); %#ok<AGROW>
    end
end
[u, ~, g] = unique(R, 'stable');
c.rxns = u;
n = numel(u);
c.kcats = zeros(n, 1); c.sd = zeros(n, 1); c.n = zeros(n, 1); c.spread = zeros(n, 1);
for r = 1:n
    rows = g == r;
    w = W(rows);
    c.kcats(r) = 10 ^ (sum(w .* X(rows)) / sum(w));
    srcs = S(rows);
    hasOkp = any(startsWith(srcs, 'okp_'));
    hasTree = any(strcmp(srcs, 'ec_tree'));
    if hasTree && hasOkp
        key = 'ec_tree+okp';
    elseif hasOkp
        key = 'okp';
    elseif hasTree
        key = 'ec_tree+homologs';
    else
        key = '';
    end
    scale = 1;
    hit = strcmp(params.scaleKeys, key);
    if any(hit)
        scale = params.scales(find(hit, 1));
    end
    c.sd(r) = sqrt(1 / sum(w) + params.measurementSd ^ 2) * scale;
    c.n(r) = sum(rows);
    c.spread(r) = max(V(rows)) - min(V(rows));
end
end

function b = evidenceBin(source, detail)
% The evidence bin of a source's value, from its detail column.
b = '';
if strcmp(source, 'ec_tree')
    n = regexp(detail, 'n=(\d+)', 'tokens', 'once');
    depth = regexp(detail, 'depth=(\d)', 'tokens', 'once');
    if ~isempty(depth) && strcmp(depth{1}, '4') && ~isempty(n)
        v = str2double(n{1});
        edges = [2 5 20 100];
        labels = {'1-2', '3-5', '6-20', '21-100', '>100'};
        b = ['ec4:' labels{find([v <= edges, true], 1)}];
    elseif ~isempty(depth)
        b = ['ec' depth{1}];
    else
        b = 'ec0';
    end
elseif strcmp(source, 'homologs')
    t = regexp(detail, 'identity=([\d.]+)', 'tokens', 'once');
    if ~isempty(t) && str2double(t{1}) >= 0.95
        b = 'exact'; % the same protein: an assumed bin, see kcatCombineParams.tsv
    elseif ~isempty(t)
        labels = {'40', '50', '60', '70', '80', '90'};
        b = labels{find([str2double(t{1}) <= [0.5 0.6 0.7 0.8 0.9], true], 1)};
    end
elseif startsWith(source, 'okp_')
    t = regexp(detail, 'similarity=([\d.]+)', 'tokens', 'once');
    if ~isempty(t)
        labels = {'<40', '40-60', '60-80', '>80'};
        b = labels{find([str2double(t{1}) <= [40 60 80], true], 1)};
    end
end
end
