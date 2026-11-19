function [evidence, summary, rows] = kcatEvidence(model, varargin)
% kcatEvidence  Show the kcats that every source gives for selected reactions.
%
% For each selected model.ec reaction, reads the kcat of every source file
% in data/kcats (see loadKcatSource) and puts it next to the kcat that the
% model has now and, for a tuned kcat, the kcat before tuning. The
% reactions are selected with rxns, enzymes and tunedKcats; a reaction that
% any of them selects is included. Without a selection, all model.ec
% reactions are included; a selection that is given and empty (for example
% no relaxed enzymes) selects none.
%
% The EC-tree estimate can be computed for the selected reactions (compute).
% The source files are only read, never changed.
%
% Name-Value Arguments
% --------------------
% rxns : cell, char, double or logical
%     reactions to include: names of model.ec.rxns or conventional reaction
%     names (see getEcRxnKeys; a reverse direction is named with _REV), or
%     indices or a logical vector over model.ec.rxns.
% enzymes : cell
%     UniProt accessions; every reaction catalysed by one of them is
%     included. flexEnz.uniprotIDs, the result of flexibilizeEnzConcs or
%     relaxEnzConcs, gives the enzymes whose measured concentration was
%     relaxed.
% tunedKcats : struct
%     the tuned kcats returned by sensitivityTuning; the model is the one
%     it returned.
% sources : cell
%     names of the source files, without .tsv (default: all files in
%     data/kcats, including the okp_<method>__query.tsv files that
%     kcatsForReactions wrote; {} for none).
% compute : cell
%     sources to compute for the selected reactions. 'ec_tree' is the one
%     that is supported (default none). The computed estimate replaces a
%     source file of the same name, so it follows measurements and the k
%     of the model adapter (default 5) and is not outdated.
% measurements : struct
%     BRENDA measurement table for the EC-tree estimate (default:
%     loadBrendaMeasurements).
% saveAs : char
%     path of a tab-separated file to write summary to (default: none).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% evidence : table
%     one row per kcat of a source, with the variables ecRxn (model.ec.rxns),
%     source, kcat, substrate, detail and assigned (true when the kcat is the
%     one that was assigned to the reaction: the kcat before tuning for a
%     tuned reaction, and the current kcat otherwise).
% summary : table
%     one row per selected reaction, with the variables ecRxn, enzymes,
%     eccodes, kcat (current), origin (the origin in model.ec.source), tuned
%     (logical), prior (the kcat before the first tuning, NaN when untuned),
%     fold (kcat divided by prior; NaN without a positive prior) and one
%     variable src_<source> per source, with the median of the kcats of that
%     source for the reaction (NaN without one). A conventional reaction name
%     in rxns or tunedKcats selects every isozyme row of the reaction, as
%     tunedKcats names the reactions of a light ecModel by their conventional
%     name.
% rows : double
%     the selected indices of model.ec.rxns.
%
% Examples
% --------
%     [model, tunedKcats] = sensitivityTuning(model);
%     [~, summary] = kcatEvidence(model, 'tunedKcats', tunedKcats);
%     [~, summary] = kcatEvidence(model, 'rxns', {'MAR03905'}, 'compute', {'ec_tree'});
%     [~, summary] = kcatEvidence(model, 'enzymes', flexEnz.uniprotIDs);
%
% See also
% --------
% loadKcatSource, kcatsFromECTree, kcatsForReactions, sensitivityTuning

checkOptions(varargin, {'rxns','enzymes','tunedKcats','sources','compute','measurements','saveAs','modelAdapter'}, 'kcatEvidence');
p = parseGECKOargs(varargin, { ...
    'rxns',         []; ...
    'enzymes',      []; ...
    'tunedKcats',   []; ...
    'sources',      []; ...
    'compute',      {}; ...
    'measurements', []; ...
    'saveAs',       ''; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
compute = {};
if ~isempty(p.compute)
    compute = cellstr(p.compute);
end
unknown = setdiff(compute, {'ec_tree'});
if ~isempty(unknown)
    error('kcatEvidence:compute', 'Cannot compute %s; supported: ec_tree.', strjoin(unknown, ', '))
end

given = cellfun(@char, varargin(1:2:end), 'UniformOutput', false);
rows = selectRows(model, p.rxns, p.enzymes, p.tunedKcats, ismember({'rxns','enzymes','tunedKcats'}, given));
nRows = numel(rows);

kcatDir = fullfile(modelAdapter.params.path, 'data', 'kcats');
sources = p.sources;
if isnumeric(sources) && isempty(sources)
    files = dir(fullfile(kcatDir, '*.tsv'));
    sources = sort(cellfun(@(f) f(1:end-4), {files.name}, 'UniformOutput', false));
end
sources = unique(cellstr(sources), 'stable');
if ismember('ec_tree', compute)
    sources = sources(~strcmp(sources, 'ec_tree'));
end

ecRxn = cell(0, 1); source = cell(0, 1); kcat = zeros(0, 1);
substrate = cell(0, 1); detail = cell(0, 1);
for s = 1:numel(sources)
    kcatList = loadKcatSource(model, sources{s}, 'modelAdapter', modelAdapter);
    [keep, ~] = ismember(kcatList.rxns, model.ec.rxns(rows));
    n = nnz(keep);
    ecRxn = [ecRxn; kcatList.rxns(keep)]; %#ok<AGROW>
    source = [source; repmat(sources(s), n, 1)]; %#ok<AGROW>
    kcat = [kcat; kcatList.kcats(keep)]; %#ok<AGROW>
    substrate = [substrate; fieldOrEmpty(kcatList, 'substrates', keep)]; %#ok<AGROW>
    detail = [detail; fieldOrEmpty(kcatList, 'detail', keep)]; %#ok<AGROW>
end
if ismember('ec_tree', compute) && nRows > 0
    tree = kcatsFromECTree(model, 'ecRxns', rows, 'save', false, 'measurements', p.measurements, ...
        'modelAdapter', modelAdapter);
    n = numel(tree.rxns);
    ecRxn = [ecRxn; tree.rxns(:)];
    source = [source; repmat({'ec_tree'}, n, 1)];
    kcat = [kcat; tree.kcats(:)];
    substrate = [substrate; repmat({''}, n, 1)];
    detail = [detail; tree.detail(:)];
    sources = [sources(:); {'ec_tree'}];
end

info = parseKcatSource(model.ec.source(rows));
if nRows == 0
    info = struct('origin', {}, 'prior', {}, 'tuned', {});
end
current = model.ec.kcat(rows);
current = current(:);
prior = reshape([info.prior], [], 1);
tuned = ~isnan(prior) | reshape(~cellfun(@isempty, {info.tuned}), [], 1);
reference = current;
reference(~isnan(prior)) = prior(~isnan(prior));

[~, at] = ismember(ecRxn, model.ec.rxns(rows));
assigned = at > 0 & abs(log(kcat ./ reference(max(at, 1)))) < 1e-5;
evidence = table(ecRxn, source, kcat, substrate, detail, assigned);

fold = nan(nRows, 1);
positive = prior > 0;
fold(positive) = current(positive) ./ prior(positive);
enzymes = cell(nRows, 1);
for i = 1:nRows
    enzymes{i} = strjoin(model.ec.enzymes(model.ec.rxnEnzMat(rows(i), :) ~= 0)', ';');
end
eccodes = model.ec.eccodes(rows);
eccodes = cellfun(@joinCodes, eccodes(:), 'UniformOutput', false);
summary = table(reshape(model.ec.rxns(rows), [], 1), enzymes, eccodes, current, ...
    reshape({info.origin}, [], 1), tuned, prior, fold, ...
    'VariableNames', {'ecRxn', 'enzymes', 'eccodes', 'kcat', 'origin', 'tuned', 'prior', 'fold'});
for s = 1:numel(sources)
    values = nan(nRows, 1);
    isSource = strcmp(source, sources{s});
    for i = 1:nRows
        hit = isSource & strcmp(ecRxn, model.ec.rxns{rows(i)});
        if any(hit)
            values(i) = median(kcat(hit));
        end
    end
    summary.(matlab.lang.makeValidName(['src_' sources{s}])) = values;
end

if ~isempty(p.saveAs)
    writetable(summary, char(p.saveAs), 'FileType', 'text', 'Delimiter', '\t');
end
end

function rows = selectRows(model, rxns, enzymes, tunedKcats, given)
% Indices of the model.ec rows that rxns, enzymes or tunedKcats select; an
% argument that is given and empty selects none.
nEc = numel(model.ec.rxns);
if ~any(given)
    rows = (1:nEc)';
    return
end
sel = false(nEc, 1);
baseRxns = getEcRxnKeys(model);
if given(1) && ~isempty(rxns)
    if islogical(rxns)
        if numel(rxns) ~= nEc
            error('kcatEvidence:rxns', 'A logical rxns needs %d entries, one per row of model.ec.rxns.', nEc)
        end
        sel(rxns) = true;
    elseif isnumeric(rxns)
        if any(rxns < 1 | rxns > nEc | rxns ~= round(rxns))
            error('kcatEvidence:rxns', 'rxns has an index outside model.ec.rxns (%d rows).', nEc)
        end
        sel(rxns) = true;
    else
        sel = sel | namedRows(model, baseRxns, cellstr(rxns), 'rxns');
    end
end
if given(3) && ~isempty(tunedKcats.rxns)
    sel = sel | namedRows(model, baseRxns, tunedKcats.rxns, 'tunedKcats.rxns');
end
if given(2) && ~isempty(enzymes)
    enzymes = cellstr(enzymes);
    missing = setdiff(enzymes, model.ec.enzymes);
    if ~isempty(missing)
        error('kcatEvidence:enzymes', 'Not in model.ec.enzymes: %s.', strjoin(missing, ', '))
    end
    sel = sel | any(model.ec.rxnEnzMat(:, ismember(model.ec.enzymes, enzymes)) ~= 0, 2);
end
rows = find(sel);
end

function sel = namedRows(model, baseRxns, names, label)
% Rows named by their model.ec.rxns entry or by their conventional reaction.
names = names(:);
sel = ismember(model.ec.rxns(:), names) | ismember(baseRxns(:), names);
missing = names(~ismember(names, model.ec.rxns) & ~ismember(names, baseRxns));
if ~isempty(missing)
    error('kcatEvidence:rxns', 'Not reactions of model.ec (%s): %s.', label, strjoin(missing', ', '))
end
end

function out = fieldOrEmpty(s, field, keep)
if isfield(s, field)
    out = reshape(s.(field)(keep), [], 1);
else
    out = repmat({''}, nnz(keep), 1);
end
end

function s = joinCodes(c)
if iscell(c)
    s = strjoin(c(:)', ';');
elseif ischar(c) || isstring(c)
    s = char(c);
else
    s = '';
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
