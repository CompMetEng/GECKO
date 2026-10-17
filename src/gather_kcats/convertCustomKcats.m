function [file, report] = convertCustomKcats(model, varargin)
% convertCustomKcats  Convert data/customKcats.tsv into data/kcats/custom.tsv.
%
% Writes the custom kcats of a GECKO 3 project as the curated kcat source
% file data/kcats/custom.tsv, which collectKcats reads first. The kcats are
% those applyCustomKcats gives the ecModel, including which of two entries
% for the same reaction wins: one row per reaction, keyed by conventional
% reaction and enzyme set (see getEcRxnKeys), with the entry's notes in the
% note column. The stoicho column is left out: complex stoichiometry comes
% from applyComplexData. Unlike applyCustomKcats, custom.tsv also applies to
% light ecModels. geckopy has the same (convert_custom_kcats).
%
% Parameters
% ----------
% model : struct
%     the full ecModel the custom kcats are for (with EC data).
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% customKcats : char
%     the file to convert (default: data/customKcats.tsv in the project).
% outFile : char
%     the file to write (default: data/kcats/custom.tsv in the project). An
%     existing file is not overwritten.
%
% Returns
% -------
% file : char
%     the file written.
% report : struct
%     - rows: number of rows written;
%     - partial: reactions whose enzymes only partly match an entry, which
%       applyCustomKcats does not change either (its notMatch).
%
% Examples
% --------
%     convertCustomKcats(ecModel);   % then collectKcats/assignKcats read custom.tsv
%
% See also
% --------
% collectKcats, applyCustomKcats, getEcRxnKeys, upgradeGECKOproject

checkOptions(varargin, {'modelAdapter','customKcats','outFile'}, 'convertCustomKcats');
p = parseGECKOargs(varargin, {'modelAdapter', []; 'customKcats', []; 'outFile', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
if model.ec.geckoLight
    error('convertCustomKcats:light', ['applyCustomKcats does not work on light ecModels; give ' ...
        'the full ecModel. The custom.tsv written applies to light ecModels too.'])
end
inFile = p.customKcats;
if isempty(inFile)
    inFile = fullfile(modelAdapter.params.path, 'data', 'customKcats.tsv');
end
file = p.outFile;
if isempty(file)
    file = fullfile(modelAdapter.params.path, 'data', 'kcats', 'custom.tsv');
end
if ~isfile(inFile)
    error('convertCustomKcats:missing', 'No custom kcats file %s.', inFile)
end
if isfile(file)
    error('convertCustomKcats:exists', '%s exists and is not overwritten.', file)
end

base = model;
base.ec.notes(:) = {''};
[applied, rxnUpdated, notMatch] = applyCustomKcats(base, 'customKcats', inFile, 'modelAdapter', modelAdapter);
idx = find(ismember(applied.ec.rxns, rxnUpdated));
[baseRxns, enzymeSets] = getEcRxnKeys(applied);
rows = strings(numel(idx), 9);
for k = 1:numel(idx)
    i = idx(k);
    note = regexprep(applied.ec.notes{i}, '[\t\r\n]', ' ');
    rows(k, :) = [string(baseRxns{i}), string(enzymeSets{i}), "", formatNumber(applied.ec.kcat(i)), ...
        "", "custom", "", "", string(note)];
end
% Isozymic copies of one conventional reaction share a key and a value.
if ~isempty(rows)
    [~, first] = unique(join(rows(:, 1:2), char(1), 2), 'stable');
    rows = rows(sort(first), :);
    [~, order] = sort(join(rows(:, [1 2 3 6 4 5 7 8 9]), char(1), 2));
    rows = rows(order, :);
end

folder = fileparts(file);
if ~isfolder(folder)
    mkdir(folder);
end
gecko = strtrim(fileread(fullfile(findGECKOroot(), 'version.txt')));
[~, name, ext] = fileparts(inFile);
fid = fopen(file, 'w', 'n', 'UTF-8');
fprintf(fid, '# source: custom\n');
fprintf(fid, '# generator: convertCustomKcats (GECKO %s)\n', gecko);
fprintf(fid, '# parameters: {}\n');
fprintf(fid, '# data: %s\n', [name ext]);
fprintf(fid, '# kind: curated\n');
fprintf(fid, '# created: %s\n', char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
fprintf(fid, 'rxn\tenzymes\tsubstrate\tkcat\tlog10_sd\tsource\tdetail\tec_code\tnote\n');
if ~isempty(rows)
    fprintf(fid, '%s\n', join(rows, sprintf('\t'), 2));
end
fclose(fid);
report.rows = size(rows, 1);
report.partial = notMatch;
fprintf('%d custom kcats written to %s.\n', report.rows, file);
end

function t = formatNumber(x)
% The shortest text that reads back as x, as saveKcatSource writes numbers.
if x == round(x) && abs(x) < 1e15
    t = string(sprintf('%d', x));
    return
end
for digits = 1:17
    t = string(sprintf(sprintf('%%.%dg', digits), x));
    if str2double(t) == x
        return
    end
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
