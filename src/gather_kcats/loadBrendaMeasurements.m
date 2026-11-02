function measurements = loadBrendaMeasurements(varargin)
% loadBrendaMeasurements  Load the BRENDA measurement table.
%
% Reads measurements.tsv.xz, or measurements.tsv when there is no .xz file,
% from the BRENDA database folder of the model adapter (next to kcat.tsv,
% sa.tsv and mw.tsv). A .xz file is decompressed with the xz command, which
% has to be on the system path. The table has one row per wild-type
% BRENDA kcat with an EC number, after the outlier screen, with the
% substrate, organism, UniProt accession and sequence id of each
% measurement. geckopy reads the same file (load_brenda_measurements).
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% folder : char
%     folder with measurements.tsv.xz or measurements.tsv, instead of the
%     adapter's BRENDA folder.
%
% Returns
% -------
% measurements : struct
%     one field per column of the file (ec, substrate, substrate_inchikey,
%     organism, taxid, uniprot, seq_id, temperature, ph and pmid as cell
%     arrays of char; kcat as double), plus log10kcat and release (the
%     release line of the file, e.g. 'BRENDA release 2026.1').
%
% Examples
% --------
%     measurements = loadBrendaMeasurements();
%
% See also
% --------
% kcatsFromECTree, ecTreeEstimate, loadBRENDAdata

checkOptions(varargin, {'modelAdapter','folder'}, 'loadBrendaMeasurements');
p = parseGECKOargs(varargin, {'modelAdapter', []; 'folder', []});
folder = p.folder;
if isempty(folder)
    modelAdapter = p.modelAdapter;
    if isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end
    folder = modelAdapter.getBrendaDBFolder();
end
xzFile = fullfile(folder, 'measurements.tsv.xz');
file = fullfile(folder, 'measurements.tsv');
tempFile = '';
if isfile(xzFile)
    tempFile = [tempname() '.tsv'];
    cleanup = onCleanup(@() deleteIfExists(tempFile)); %#ok<NASGU>
    [status, output] = system(sprintf('xz -dc "%s" > "%s"', xzFile, tempFile));
    if status ~= 0
        error('loadBrendaMeasurements:xz', ['Could not decompress %s with the xz command ' ...
            '(install xz or decompress the file to measurements.tsv): %s'], xzFile, strtrim(output))
    end
    file = tempFile;
elseif ~isfile(file)
    error('loadBrendaMeasurements:missing', ['No BRENDA measurement table (measurements.tsv) in %s. ' ...
        'It comes with the BRENDA files of GECKO (databases folder).'], folder)
end
fid = fopen(file, 'r', 'n', 'UTF-8');
first = fgetl(fid);
header = strsplit(fgetl(fid), '\t');
data = textscan(fid, repmat('%s', 1, numel(header)), 'Delimiter', '\t', 'Whitespace', '', ...
    'ReturnOnError', false);
fclose(fid);
if ~all(ismember({'ec', 'kcat'}, header))
    error('loadBrendaMeasurements:columns', '%s lacks the columns ec and kcat.', file)
end
measurements = struct();
for c = 1:numel(header)
    measurements.(matlab.lang.makeValidName(header{c})) = data{c};
end
measurements.kcat = str2double(measurements.kcat);
keep = measurements.kcat > 0 & ~cellfun(@isempty, measurements.ec);
for f = fieldnames(measurements)'
    measurements.(f{1}) = measurements.(f{1})(keep);
end
measurements.ec = strtrim(measurements.ec);
measurements.log10kcat = log10(measurements.kcat);
tok = regexp(first, '^#\s*(BRENDA release [^ ]+)', 'tokens', 'once');
measurements.release = '';
if ~isempty(tok)
    measurements.release = tok{1};
end
end

function deleteIfExists(file)
if isfile(file)
    delete(file)
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
