function [kcatList, report] = kcatsFromWILDkCAT(model, filename, varargin)
% kcatsFromWILDkCAT  Read the kcats of a WILDkCAT run and keep them as a source file.
%
% WILDkCAT (https://github.com/sysbiolux/WILDkCAT) retrieves kcat values
% from BRENDA and SABIO-RK for each reaction, EC number, enzyme and
% direction of a conventional GEM, and fills gaps with CataPro
% predictions. It keeps one value per row, the measurement with the lowest
% penalty score (0 is an exact match of organism, enzyme, substrate,
% temperature and pH; higher scores are less close matches). Run WILDkCAT
% on the conventional GEM the ecModel was made from, then give its output
% table (kcat_full.tsv, or kcat_retrieved.tsv without predictions) to this
% function.
%
% A WILDkCAT row is matched to the model.ec rows of the same reaction (its
% _REV for direction 'reverse') whose genes (model.ec.genes of the enzymes
% in model.ec.rxnEnzMat) are the genes of the row. When several WILDkCAT
% rows match the same model.ec row, the one with the lowest penalty score
% is kept, and of those the highest kcat, as WILDkCAT itself chooses.
% Rows without a kcat, and rows that match no model.ec row, are counted in
% the report. The result is written to data/kcats/wildkcat.tsv (see
% saveKcatSource), which a stage of [kcat_assignment] in
% model_adapter.toml can name as a source. geckopy has the same
% (kcats_from_wildkcat).
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure).
% filename : char
%     the WILDkCAT output table (kcat_full.tsv or kcat_retrieved.tsv).
%
% Name-Value Arguments
% --------------------
% maxPenalty : double
%     keep only values with this penalty score or lower (default Inf: all
%     values). A predicted value has the penalty score of the measurement
%     it replaced.
% predictions : logical
%     keep the CataPro predictions, the rows whose db is 'catapro'
%     (default true).
% save : logical
%     write the source file (default true).
% saveAs : char
%     name of the source file, without .tsv (default 'wildkcat').
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     the kcat values, with fields source, rxns, substrates, kcats,
%     eccodes and detail (db, penalty score, and the organism, enzyme,
%     temperature, pH and sequence identity of the measurement, as far as
%     WILDkCAT gives them).
% report : struct
%     numbers of WILDkCAT rows: rows, noValue (no kcat), excluded (by
%     maxPenalty or predictions), matched, and unmatched, the
%     '<rxn>|<direction>|<genes>' of the rows that match no model.ec row.
%
% Examples
% --------
%     kcatsFromWILDkCAT(ecModel, 'wildkcat_output/kcat_full.tsv');
%     kcatsFromWILDkCAT(ecModel, 'kcat_full.tsv', 'maxPenalty', 4, 'predictions', false);
%
% See also
% --------
% saveKcatSource, loadKcatSource, getEcRxnKeys, collectKcats

checkOptions(varargin, {'maxPenalty','predictions','save','saveAs','modelAdapter'}, 'kcatsFromWILDkCAT');
p = parseGECKOargs(varargin, { ...
    'maxPenalty',   Inf; ...
    'predictions',  true; ...
    'save',         true; ...
    'saveAs',       'wildkcat'; ...
    'modelAdapter', []});
filename = char(filename);
if ~isfile(filename)
    error('kcatsFromWILDkCAT:file', 'WILDkCAT output %s does not exist.', filename)
end

t = readWildkcatTable(filename);
required = {'rxn','direction','genes','kcat','db','penalty_score'};
missing = setdiff(required, fieldnames(t));
if ~isempty(missing)
    error('kcatsFromWILDkCAT:columns', '%s lacks the WILDkCAT column(s) %s.', ...
        filename, strjoin(missing, ', '))
end
n = numel(t.rxn);
kcat = str2double(t.kcat);
penalty = str2double(t.penalty_score);
db = lower(strtrim(t.db));

report.rows = n;
hasValue = isfinite(kcat) & kcat > 0;
report.noValue = sum(~hasValue);
use = hasValue & ~(penalty > p.maxPenalty);
if ~p.predictions
    use = use & ~strcmp(db, 'catapro');
end
report.excluded = sum(hasValue & ~use);

% Keys of the model.ec rows: conventional reaction and sorted gene set.
baseRxns = getEcRxnKeys(model);
geneSets = repmat({''}, numel(model.ec.rxns), 1);
[r, c] = find(model.ec.rxnEnzMat);
if ~isempty(r)
    perRxn = accumarray(r(:), c(:), [numel(model.ec.rxns) 1], @(x) {x});
    for i = find(~cellfun(@isempty, perRxn))'
        geneSets{i} = strjoin(unique(model.ec.genes(perRxn{i}))', ';');
    end
end
ecKeys = strcat(baseRxns, '|', geneSets);

% Keys of the WILDkCAT rows.
rxn = strtrim(t.rxn);
isRev = strcmpi(strtrim(t.direction), 'reverse');
rxn(isRev) = strcat(rxn(isRev), '_REV');
genes = cellfun(@(g) strjoin(unique(strtrim(strsplit(g, ';'))), ';'), t.genes, 'UniformOutput', false);
genes = regexprep(genes, '^;|;$', '');
rowKeys = strcat(rxn, '|', genes);

matched = false(n, 1);
best = zeros(numel(model.ec.rxns), 1);  % chosen WILDkCAT row per model.ec row
for i = find(use)'
    hits = find(strcmp(ecKeys, rowKeys{i}) & ~cellfun(@isempty, geneSets));
    if isempty(hits)
        continue
    end
    matched(i) = true;
    for h = hits(:)'
        b = best(h);
        if b == 0 || penalty(i) < penalty(b) || (penalty(i) == penalty(b) && kcat(i) > kcat(b))
            best(h) = i;
        end
    end
end
report.matched = sum(matched);
unmatched = use & ~matched;
report.unmatched = unique(strcat(t.rxn(unmatched), '|', t.direction(unmatched), '|', t.genes(unmatched)));

ecRows = find(best);
src = best(ecRows);
kcatList.source = 'wildkcat';
kcatList.rxns = model.ec.rxns(ecRows);
kcatList.substrates = column(t, 'kcat_substrate', src);
kcatList.kcats = kcat(src);
kcatList.eccodes = column(t, 'ec_code', src);
kcatList.detail = cell(numel(src), 1);
fields = {'db','db'; 'penalty','penalty_score'; 'organism','kcat_organism'; ...
    'enzyme','kcat_enzyme'; 'temperature','kcat_temperature'; 'pH','kcat_ph'; ...
    'identity','kcat_id_percent'};
for k = 1:numel(src)
    parts = {};
    for f = 1:size(fields, 1)
        v = column(t, fields{f, 2}, src(k));
        v = strtrim(v{1});
        if ~isempty(v) && ~strcmpi(v, 'nan')
            parts{end+1} = sprintf('%s=%s', fields{f, 1}, strrep(v, ';', ',')); %#ok<AGROW>
        end
    end
    kcatList.detail{k} = strjoin(parts, '; ');
end

fprintf('WILDkCAT: %d rows, %d without a kcat, %d excluded, %d matched to %d model.ec reactions, %d unmatched.\n', ...
    report.rows, report.noValue, report.excluded, report.matched, numel(ecRows), sum(unmatched));
if p.save
    modelAdapter = p.modelAdapter;
    if isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end
    [~, base, ext] = fileparts(filename);
    parameters = struct('file', [base ext], 'max_penalty', p.maxPenalty, 'predictions', logical(p.predictions));
    if isinf(p.maxPenalty)
        parameters = rmfield(parameters, 'max_penalty');
    end
    file = saveKcatSource(model, kcatList, p.saveAs, 'modelAdapter', modelAdapter, ...
        'generator', 'kcatsFromWILDkCAT', 'parameters', parameters);
    fprintf('WILDkCAT kcats for %d reactions written to %s\n', numel(unique(kcatList.rxns)), file);
end
end

function t = readWildkcatTable(filename)
% The tab-separated table as a struct of n x 1 cell columns, all text.
lines = splitlines(string(fileread(filename)));
lines = lines(strlength(strtrim(lines)) > 0);
if isempty(lines)
    error('kcatsFromWILDkCAT:file', '%s is empty.', filename)
end
header = strtrim(split(lines(1), sprintf('\t')));
cells = strings(numel(lines) - 1, numel(header));
for i = 2:numel(lines)
    parts = split(lines(i), sprintf('\t'));
    m = min(numel(parts), numel(header));
    cells(i - 1, 1:m) = parts(1:m);
end
t = struct();
for j = 1:numel(header)
    t.(matlab.lang.makeValidName(char(header(j)))) = cellstr(cells(:, j));
end
end

function v = column(t, name, idx)
% Entries idx of a column, or empty text when WILDkCAT did not write it.
if isfield(t, name)
    v = t.(name)(idx);
else
    v = repmat({''}, numel(idx), 1);
end
v = v(:);
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
