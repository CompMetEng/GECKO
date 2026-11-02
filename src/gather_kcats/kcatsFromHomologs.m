function kcatList = kcatsFromHomologs(model, varargin)
% kcatsFromHomologs  Estimate kcats from measured homologs and keep them as a source file.
%
% For each reaction, enzyme and substrate (the substrates DLKcat would use,
% see extractEnzymeSubstratePairs), the enzyme's sequence is searched with
% DIAMOND against the sequences of the BRENDA measurements (at least 40%
% identity, 80% coverage of both sequences). Each homolog measured on the
% same substrate (same InChIKey, without the protonation flag) contributes
% the median of its measurements; the estimate is the mean of the k most
% identical such homologs, weighted by identity. Homologs at 95% identity or
% more are measurements of (nearly) the same protein: when there are any,
% the estimate is the median of those alone. For an enzyme complex the
% estimate is the median over its subunits. Model metabolites are matched to
% InChIKeys through substrates.tsv of the BRENDA files: their ChEBI
% annotations, then MetaNetX annotations (model.metMiriams), then their
% names. Reactions without a homolog measured on one of their substrates get
% no value. The result is written to data/kcats/homologs.tsv (see
% saveKcatSource). DIAMOND is RAVEN's (getDiamond; downloaded when
% missing). geckopy has the same (kcats_from_homologs).
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure), with enzyme
%     sequences (model.ec.sequence).
%
% Name-Value Arguments
% --------------------
% k : double
%     the number of most identical homologs averaged (default 10).
% minIdentity : double
%     DIAMOND identity threshold, as a fraction (default 0.4).
% coverage : double
%     DIAMOND coverage threshold of both sequences, as a fraction (default 0.8).
% threads : double
%     DIAMOND threads (default: DIAMOND's).
% save : logical
%     write the source file (default true).
% saveAs : char
%     name of the source file, without .tsv (default 'homologs').
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     the kcat values, with fields source, rxns, substrates, kcats and
%     detail (the identity of the nearest homolog, identity, and the number
%     of homologs used, n).
%
% Examples
% --------
%     kcatsFromHomologs(ecModel);
%
% See also
% --------
% kcatsFromECTree, extractEnzymeSubstratePairs, saveKcatSource, collectKcats

checkOptions(varargin, {'k','minIdentity','coverage','threads','save','saveAs','modelAdapter'}, 'kcatsFromHomologs');
p = parseGECKOargs(varargin, { ...
    'k',            10; ...
    'minIdentity',  0.4; ...
    'coverage',     0.8; ...
    'threads',      []; ...
    'save',         true; ...
    'saveAs',       'homologs'; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

measurements = loadBrendaMeasurements('folder', fileparts(brendaFile(modelAdapter, 'measurements.tsv')));
substrates = readSubstrates(brendaFile(modelAdapter, 'substrates.tsv'));
sequencesFile = brendaFile(modelAdapter, 'sequences.fasta');

pairs = extractEnzymeSubstratePairs(model, 'modelAdapter', modelAdapter);
iks = metaboliteInchikeys(model, substrates);
pairIk = iks(pairs.met);
seqs = model.ec.sequence(pairs.protein);
keep = ~cellfun(@isempty, pairIk) & ~cellfun(@isempty, seqs);
rxnIdx = pairs.ecRxn(keep); metIdx = pairs.met(keep); pairIk = pairIk(keep); seqs = seqs(keep);

kcatList.source = 'homologs';
kcatList.rxns = cell(0, 1);
kcatList.substrates = cell(0, 1);
kcatList.kcats = zeros(0, 1);
kcatList.detail = cell(0, 1);
if ~isempty(seqs)
    [useqs, ~, seqNo] = unique(seqs, 'stable');
    queryIds = cellstr(compose('q%d', (0:numel(useqs)-1)'));
    hits = diamondHits(queryIds, useqs, sequencesFile, p.minIdentity, p.coverage, p.threads);
    % Median log10 kcat per measured sequence and substrate.
    hasSeq = ~cellfun(@isempty, measurements.seq_id) & ~cellfun(@isempty, measurements.substrate_inchikey);
    key = strcat(measurements.seq_id(hasSeq), '|', cellfun(@(s) s(1:min(25, end)), ...
        measurements.substrate_inchikey(hasSeq), 'UniformOutput', false));
    [ukey, ~, g] = unique(key);
    med = containers.Map(ukey, num2cell(accumarray(g, measurements.log10kcat(hasSeq), [], @median)));
    % Estimate per enzyme sequence and substrate.
    [up, ~, pairNo] = unique(strcat(queryIds(seqNo), '|', pairIk), 'stable');
    est = nan(numel(up), 1); best = nan(numel(up), 1); nHom = zeros(numel(up), 1);
    for u = 1:numel(up)
        parts = strsplit(up{u}, '|');
        rows = find(strcmp(hits.query, parts{1}));
        found = zeros(0, 2);
        for r = rows(:)'
            mk = [hits.target{r} '|' parts{2}];
            if isKey(med, mk)
                found(end+1, :) = [hits.identity(r), med(mk)]; %#ok<AGROW>
            end
        end
        if isempty(found)
            continue
        end
        exact = found(found(:, 1) >= 0.95, 2);
        best(u) = found(1, 1);
        if ~isempty(exact)
            est(u) = median(exact);
            nHom(u) = numel(exact);
        else
            found = found(1:min(p.k, end), :);
            est(u) = sum(found(:, 1) .* found(:, 2)) / sum(found(:, 1));
            nHom(u) = size(found, 1);
        end
    end
    % One row per reaction and substrate: the median over its enzymes.
    [ur, ~, rs] = unique(strcat(model.ec.rxns(rxnIdx), '|', model.metNames(metIdx)), 'stable');
    for r = 1:numel(ur)
        members = unique(pairNo(rs == r));
        members = members(~isnan(est(members)));
        if isempty(members)
            continue
        end
        [~, b] = max(best(members));
        first = find(rs == r, 1);
        kcatList.rxns{end+1, 1} = model.ec.rxns{rxnIdx(first)};
        kcatList.substrates{end+1, 1} = model.metNames{metIdx(first)};
        kcatList.kcats(end+1, 1) = 10 ^ median(est(members));
        kcatList.detail{end+1, 1} = sprintf('identity=%.3f; n=%d', best(members(b)), nHom(members(b)));
    end
end

if p.save
    parameters = struct('coverage', p.coverage, 'k', p.k, 'min_identity', p.minIdentity);
    file = saveKcatSource(model, kcatList, p.saveAs, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromHomologs', 'parameters', parameters, 'data', measurements.release);
    fprintf('Homolog kcats for %d reactions written to %s\n', numel(unique(kcatList.rxns)), file);
end
end

function file = brendaFile(modelAdapter, name)
% A file of the adapter's BRENDA folder, or of GECKO's BRENDA release.
file = fullfile(modelAdapter.getBrendaDBFolder(), name);
if ~isfile(file)
    file = fullfile(fetchBrendaAsset({name}), name);
end
end

function s = readSubstrates(file)
fid = fopen(file, 'r', 'n', 'UTF-8');
c = textscan(fid, '%s %s %s', 'Delimiter', '\t', 'CommentStyle', '#', 'HeaderLines', 0);
fclose(fid);
header = strcmp(c{1}, 'inchikey');
s.inchikey = c{1}(~header); s.kind = c{2}(~header); s.value = c{3}(~header);
end

function iks = metaboliteInchikeys(model, substrates)
% InChIKey (first 25 characters) per model metabolite, '' when unknown:
% from its ChEBI annotations, then its MetaNetX annotations, then its name.
maps = struct();
for kind = {'chebi', 'metanetx', 'name'}
    k = strcmp(substrates.kind, kind{1});
    [values, first] = unique(substrates.value(k), 'stable');
    ik = substrates.inchikey(k);
    maps.(kind{1}) = containers.Map('KeyType', 'char', 'ValueType', 'char');
    if ~isempty(values)
        maps.(kind{1}) = containers.Map(values, ik(first));
    end
end
iks = repmat({''}, numel(model.mets), 1);
for i = 1:numel(model.mets)
    for pair = {{'chebi', 'chebi'}, {'metanetx', 'metanetx.chemical'}}
        kind = pair{1}{1};
        for v = annotations(model, i, pair{1}{2})
            if isKey(maps.(kind), v{1})
                iks{i} = maps.(kind)(v{1});
                break
            end
        end
        if ~isempty(iks{i})
            break
        end
    end
    if isempty(iks{i}) && isKey(maps.name, compactName(model.metNames{i}))
        iks{i} = maps.name(compactName(model.metNames{i}));
    end
end
end

function values = annotations(model, i, key)
% The identifiers of one annotation of a metabolite, without their prefix.
values = {};
if ~isfield(model, 'metMiriams') || numel(model.metMiriams) < i || ~isstruct(model.metMiriams{i})
    return
end
m = model.metMiriams{i};
hit = strcmpi(m.name, key);
values = regexprep(strtrim(m.value(hit)), '^.*:', '');
values = values(:)';
end

function s = compactName(name)
% A name as the compact keys of substrates.tsv: without markup, Greek letters
% spelled out, accents dropped, lower case, only letters, digits and +.
s = regexprep(char(name), '<[^>]+>', '');
s = strrep(strrep(strrep(strrep(strrep(s, '&amp;', '&'), '&lt;', '<'), '&gt;', '>'), '&quot;', '"'), '&#39;', '''');
greek = {char(945), 'alpha'; char(946), 'beta'; char(947), 'gamma'; char(948), 'delta'; char(949), 'epsilon'; ...
    char(969), 'omega'; char(954), 'kappa'; char(955), 'lambda'; char(956), 'mu'; char(968), 'psi'; ...
    char(958), 'xi'; char(950), 'zeta'};
for g = 1:size(greek, 1)
    s = strrep(strrep(s, greek{g, 1}, greek{g, 2}), upper(greek{g, 1}), greek{g, 2});
end
s = char(javaMethod('normalize', 'java.text.Normalizer', java.lang.String(s), ...
    javaMethod('valueOf', 'java.text.Normalizer$Form', 'NFKD')));
s = s(double(s) < 128);
s = regexprep(lower(s), '[^a-z0-9+]', '');
end

function hits = diamondHits(ids, sequences, targetFasta, minIdentity, coverage, threads)
% DIAMOND hits of the query sequences against a FASTA file: query, target,
% identity (a fraction), best identity first per query, then by target.
ravenPath = findRAVENroot();
if ismac
    binEnd = '.mac';
elseif ispc
    binEnd = '.exe';
else
    binEnd = '';
end
exe = fullfile(ravenPath, 'src', 'external', 'diamond', ['diamond' binEnd]);
if ~isfile(exe)
    downloadRavenBinaries({'diamond'});
end
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
fid = fopen(fullfile(tmp, 'queries.fasta'), 'w');
entries = [cellstr(ids(:))'; cellstr(sequences(:))'];
fprintf(fid, '>%s\n%s\n', entries{:});
fclose(fid);
run = @(args) system(['"' exe '" ' args]);
[status, out] = run(sprintf('makedb --in "%s" --db "%s" --quiet', targetFasta, fullfile(tmp, 'targets')));
if status ~= 0
    error('kcatsFromHomologs:diamond', 'DIAMOND makedb failed: %s', out)
end
args = sprintf(['blastp --query "%s" --db "%s" --out "%s" --ultra-sensitive --id %g --query-cover %g ' ...
    '--subject-cover %g --evalue 0.001 --max-target-seqs 0 --outfmt 6 qseqid sseqid pident --quiet'], ...
    fullfile(tmp, 'queries.fasta'), fullfile(tmp, 'targets'), fullfile(tmp, 'hits.tsv'), ...
    100 * minIdentity, 100 * coverage, 100 * coverage);
if ~isempty(threads)
    args = sprintf('%s --threads %d', args, threads);
end
[status, out] = run(args);
if status ~= 0
    error('kcatsFromHomologs:diamond', 'DIAMOND blastp failed: %s', out)
end
fid = fopen(fullfile(tmp, 'hits.tsv'), 'r');
c = textscan(fid, '%s %s %f', 'Delimiter', '\t');
fclose(fid);
q = c{1}; t = c{2}; ident = c{3} / 100;
% Best identity per query and target.
[pairKeys, ~, g] = unique(strcat(q, '|', t));
ident = accumarray(g, ident, [], @max);
parts = split(pairKeys, '|');
if isempty(pairKeys)
    parts = cell(0, 2);
end
parts = reshape(parts, [], 2);
T = table(parts(:, 1), parts(:, 2), ident, 'VariableNames', {'query', 'target', 'identity'});
T = sortrows(T, {'query', 'identity', 'target'}, {'ascend', 'descend', 'ascend'});
hits = struct('query', {T.query}, 'target', {T.target}, 'identity', T.identity);
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
