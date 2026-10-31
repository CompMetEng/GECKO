function [kcatList, report] = loadKcatSource(model, name, varargin)
% loadKcatSource  Read a kcat source file and map its rows to an ecModel.
%
% Reads data/kcats/<name>.tsv (see saveKcatSource) and maps every row to
% the model.ec reactions it applies to, by conventional reaction and
% enzyme set (see getEcRxnKeys):
%
% 1. A row with rxn and enzymes applies to the model.ec rows with that
%    conventional reaction and that enzyme set.
% 2. A row with only enzymes applies to every model.ec row with that
%    enzyme set; a row with only rxn to every model.ec row of that
%    conventional reaction (curated rows; saveKcatSource always writes
%    both).
% 3. A row whose kcat is missing, not a number, infinite, or not above 0
%    is reported as invalid and not applied.
% 4. A row that matches nothing is reported, never dropped silently. If a
%    model.ec row has at least half of the row's enzymes, that partial
%    match is reported as well, but not applied.
%
% The enzymes of a row may be listed in any order.
%
% The file is read as UTF-8 (a leading byte-order mark is ignored). The
% '#' lines before the column header are the header; '#' lines and blank
% lines further down are skipped, as are rows without a filled field (a row
% of empty quoted fields included). Fields wrapped in double quotes, as a
% spreadsheet may save them, are unwrapped; a field that opens a double quote
% and does not close it is an error. kcat and log10_sd must be plain
% numbers (12, 0.5, 1e-05); other text, such as '5,5' with a decimal comma,
% is reported as invalid. A row with more filled fields than there are
% columns is an error.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
% name : char
%     name of the source (data/kcats/<name>.tsv in the project folder), or
%     the path to a kcat source file.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% kcatList : struct
%     kcat values per model.ec reaction, as used by mergeKcats and
%     assignKcatValues: source (the source name), rxns, kcats,
%     substrates, eccodes, wildcardLvl and origin (from the detail
%     column), kcatSource (the per-row source), log10sd, detail (the detail
%     column), and line (the line of the file each entry comes from).
% report : struct
%     - invalid: table of the file rows with a missing or invalid kcat;
%     - unmatched: table of the file rows that apply to no model.ec row;
%     - partial: table of rows that only partly match an enzyme set, with
%       the model.ec reaction and its enzymes;
%     - header: the '#' header lines of the file, as field: value (for
%       instance header.source, header.parameters);
%     - records: the rows of the file as read, one field per column (kcat
%       and log10_sd as double, other columns as cell of char) and line,
%       the line number of each row.
%
% Examples
% --------
%     kcatList = loadKcatSource(ecModel, 'brenda');
%     [kcatList, report] = loadKcatSource(ecModel, 'custom');
%
% See also
% --------
% saveKcatSource, kcatsFromBrenda, getEcRxnKeys

checkOptions(varargin, {'modelAdapter'}, 'loadKcatSource');
p = parseGECKOargs(varargin, {'modelAdapter', []});

if isfile(name)
    filename = name;
else
    modelAdapter = p.modelAdapter;
    if isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end
    filename = fullfile(modelAdapter.params.path, 'data', 'kcats', [name '.tsv']);
    if ~isfile(filename)
        error('loadKcatSource:missing', ['No kcat source file %s. Make it with the ' ...
            'kcatsFrom<Source> function of this source (kcatsFromBrenda, kcatsFromDLKcat, ' ...
            'kcatsFromOKP), or with saveKcatSource.'], filename)
    end
end
[records, meta] = readSourceFile(filename);
sourceName = '';
if isfield(meta, 'source')
    sourceName = meta.source;
end

[baseRxns, enzymeSets] = getEcRxnKeys(model);
nEc = numel(baseRxns);
nRec = numel(records.kcat);

% Lookup tables: key -> model.ec rows.
byPair = groupIndex(strcat(baseRxns, '|', enzymeSets));
byBase = groupIndex(baseRxns);
byEnz  = groupIndex(enzymeSets);
% accession -> model.ec rows containing it, for partial matches.
accRows = cell(nEc, 1);
for i = 1:nEc
    if ~isempty(enzymeSets{i})
        accRows{i} = strsplit(enzymeSets{i}, ' + ');
    end
end
accIdx = repelem((1:nEc)', cellfun(@numel, accRows));
byAcc = groupIndex([accRows{:}]', accIdx);

invalid = ~(records.kcat > 0) | ~isfinite(records.kcat);
matchRec = cell(nRec, 1);
matchEc  = cell(nRec, 1);
unmatched = false(nRec, 1);
partialRec = zeros(0, 1); partialEc = zeros(0, 1);
for r = find(~invalid)'
    rxn = records.rxn{r};
    enz = normalizeEnzymes(records.enzymes{r});
    if ~isempty(rxn) && ~isempty(enz)
        hit = lookup(byPair, [rxn '|' enz]);
    elseif ~isempty(enz)
        hit = lookup(byEnz, enz);
    elseif ~isempty(rxn)
        hit = lookup(byBase, rxn);
    else
        hit = [];
    end
    if isempty(hit)
        unmatched(r) = true;
        if ~isempty(enz)
            wanted = strsplit(enz, ' + ');
            candidates = unique(cell2mat(cellfun(@(a) lookup(byAcc, a), wanted(:), 'UniformOutput', false)));
            if ~isempty(rxn)
                candidates = intersect(candidates, lookup(byBase, rxn));
            end
            for c = candidates(:)'
                if mean(ismember(wanted, accRows{c})) >= 0.5
                    partialRec(end+1, 1) = r; %#ok<AGROW>
                    partialEc(end+1, 1) = c; %#ok<AGROW>
                end
            end
        end
    else
        matchRec{r} = repmat(r, numel(hit), 1);
        matchEc{r} = hit(:);
    end
end
rec = vertcat(matchRec{:}, zeros(0, 1));
ec  = vertcat(matchEc{:}, zeros(0, 1));

kcatList.source      = sourceName;
kcatList.rxns        = reshape(model.ec.rxns(ec), [], 1);
kcatList.kcats       = records.kcat(rec);
kcatList.substrates  = records.substrate(rec);
kcatList.eccodes     = records.ec_code(rec);
kcatList.kcatSource  = records.source(rec);
kcatList.log10sd     = records.log10_sd(rec);
kcatList.detail      = records.detail(rec);
[kcatList.wildcardLvl, kcatList.origin] = parseDetail(records.detail(rec));
kcatList.line        = records.line(rec);

report.invalid = table(records.line(invalid), records.rxn(invalid), records.enzymes(invalid), ...
    records.kcat(invalid), 'VariableNames', {'line','rxn','enzymes','kcat'});
report.unmatched = table(records.line(unmatched), records.rxn(unmatched), records.enzymes(unmatched), ...
    records.kcat(unmatched), 'VariableNames', {'line','rxn','enzymes','kcat'});
report.partial = table(records.line(partialRec), records.rxn(partialRec), records.enzymes(partialRec), ...
    reshape(model.ec.rxns(partialEc), [], 1), enzymeSets(partialEc), ...
    'VariableNames', {'line','rxn','enzymes','ecRxn','ecEnzymes'});
if any(invalid)
    printOrange(sprintf('%d of %d rows of %s have no valid kcat; see report.invalid.\n', ...
        sum(invalid), nRec, filename));
end
if any(unmatched)
    printOrange(sprintf('%d of %d rows of %s apply to no reaction of this model; see report.unmatched.\n', ...
        sum(unmatched), nRec, filename));
end
report.header = meta;
report.records = records;
end

function map = groupIndex(keys, rows)
% containers.Map from each key to the rows (default 1:n) that have it.
if nargin < 2
    rows = (1:numel(keys))';
end
map = containers.Map('KeyType', 'char', 'ValueType', 'any');
if isempty(keys)
    return
end
[u, ~, g] = unique(keys(:));
groups = accumarray(g, rows(:), [], @(x) {sort(x)});
for i = 1:numel(u)
    if ~isempty(u{i})
        map(u{i}) = groups{i};
    end
end
end

function rows = lookup(map, key)
if isKey(map, key)
    rows = map(key);
else
    rows = zeros(0, 1);
end
end

function enz = normalizeEnzymes(enz)
% Sorted and joined as getEcRxnKeys does, so hand-written rows match too.
if isempty(enz)
    return
end
parts = strtrim(strsplit(enz, '+'));
parts = unique(parts(~cellfun(@isempty, parts)));
enz = strjoin(sort(parts(:)'), ' + ');
end

function [wc, origin] = parseDetail(detail)
wc = nan(numel(detail), 1);
origin = nan(numel(detail), 1);
for i = 1:numel(detail)
    t = regexp(detail{i}, 'wc=(\d+)', 'tokens', 'once');
    if ~isempty(t)
        wc(i) = str2double(t{1});
    end
    t = regexp(detail{i}, 'origin=(\d+)', 'tokens', 'once');
    if ~isempty(t)
        origin(i) = str2double(t{1});
    end
end
end

function [records, meta] = readSourceFile(filename)
% The rows and '#' header lines of a kcat source file, as they are.
text = readUtf8(filename);
lines = regexp(text, '\r?\n', 'split');

meta = struct();
k = 1;
while k <= numel(lines) && (isempty(strtrim(lines{k})) || startsWith(unquoteField(lines{k}), '#'))
    [key, value] = parseHeaderLine(lines{k});
    if ~isempty(key)
        meta.(key) = value;
    end
    k = k + 1;
end
if k > numel(lines)
    error('loadKcatSource:noHeader', '%s has no column header.', filename)
end
columns = cellfun(@unquoteField, strsplit(lines{k}, sprintf('\t'), ...
    'CollapseDelimiters', false), 'UniformOutput', false);
if ~all(ismember({'kcat','source'}, columns)) || ~any(ismember({'rxn','enzymes'}, columns))
    error('loadKcatSource:columns', ['%s needs the columns kcat and source, and rxn ' ...
        'and/or enzymes.'], filename)
end
named = columns(~cellfun(@isempty, columns));
if numel(unique(named)) < numel(named) || any(strcmp(columns, 'line'))
    error('loadKcatSource:columns', '%s has a column name twice, or a column named line.', filename)
end
nCol = numel(columns);

trimmed = strtrim(lines);
skip = cellfun(@isempty, trimmed) | startsWith(trimmed, '#');
commentLike = find(startsWith(trimmed, '"#'));
for q = commentLike(:)'
    skip(q) = startsWith(unquoteField(trimmed{q}), '#');
end
skip(1:k) = true;
kept = find(~skip);
L = string(lines(kept));
L = L(:);
F = strings(numel(kept), nCol);
nTabs = count(L, sprintf('\t'));
regular = nTabs == nCol - 1;
if any(regular)
    % Rows with exactly one field per column, split in one go.
    F(regular, :) = reshape(split(L(regular), sprintf('\t'), 2), [], nCol);
end
for r = find(~regular)'
    fields = strsplit(char(L(r)), sprintf('\t'), 'CollapseDelimiters', false);
    if numel(fields) > nCol
        if any(~cellfun(@isempty, strtrim(fields(nCol+1:end))))
            error('loadKcatSource:fields', '%s, line %d has more fields than the %d columns.', ...
                filename, kept(r), nCol)
        end
        fields = fields(1:nCol);
    end
    fields(end+1:nCol) = {''};
    F(r, :) = string(fields);
end
F = strtrim(F);
quoted = find(startsWith(F, '"'));
unclosed = quoted(~(strlength(F(quoted)) >= 2 & endsWith(F(quoted), '"')));
if ~isempty(unclosed)
    badRow = min(mod(unclosed - 1, size(F, 1)) + 1);
    error('loadKcatSource:quote', '%s, line %d has a field with an opening double quote and no closing one.', ...
        filename, kept(badRow))
end
for q = quoted(:)'
    F(q) = unquoteField(char(F(q)));
end
filled = any(F ~= "", 2); % a row of tabs only is skipped
values = cellstr(F(filled, :));
if isempty(values)
    values = cell(0, nCol);
end
lineNo = kept(filled);
lineNo = lineNo(:);
n = numel(lineNo);

records = struct();
for c = 1:nCol
    if isempty(columns{c})
        continue
    end
    name = matlab.lang.makeValidName(columns{c});
    v = values(:, c);
    if any(strcmp(columns{c}, {'kcat','log10_sd'}))
        v = parseKcatNumber(v);
    end
    records.(name) = v;
end
for c = {'rxn','enzymes','substrate','source','detail','ec_code','note'}
    if ~isfield(records, c{1})
        records.(c{1}) = repmat({''}, n, 1);
    end
end
if ~isfield(records, 'log10_sd')
    records.log10_sd = nan(n, 1);
end
records.line = lineNo;
end

function text = readUtf8(filename)
fid = fopen(filename, 'r', 'n', 'UTF-8');
if fid < 0
    error('loadKcatSource:read', 'Cannot read %s', filename)
end
cleanup = onCleanup(@() fclose(fid));
text = fread(fid, [1 Inf], '*char');
text = regexprep(text, ['^' char(65279)], '');
end

function x = parseKcatNumber(texts)
% Plain numbers only: '5,5' or '12 /s' read as NaN.
x = nan(numel(texts), 1);
ok = ~cellfun(@isempty, regexp(texts(:), '^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$|^[-+]?Inf$', 'once'));
x(ok) = str2double(texts(ok));
end

function [key, value] = parseHeaderLine(line)
key = ''; value = '';
tok = regexp(unquoteField(line), '^\s*#\s*([A-Za-z][A-Za-z0-9_]{0,62})\s*:\s?(.*)$', 'tokens', 'once');
if ~isempty(tok)
    key = lower(tok{1});
    value = strtrim(tok{2});
end
end

function t = unquoteField(t)
t = strtrim(t);
if numel(t) >= 2 && t(1) == '"' && t(end) == '"'
    t = strtrim(strrep(t(2:end-1), '""', '"'));
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
