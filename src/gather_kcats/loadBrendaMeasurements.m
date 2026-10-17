function measurements = loadBrendaMeasurements(varargin)
% loadBrendaMeasurements  Load the BRENDA measurement table.
%
% Reads measurements.tsv from the BRENDA database folder of the model
% adapter (next to kcat.tsv, sa.tsv and mw.tsv): one row per wild-type
% BRENDA kcat with an EC number, after the outlier screen, with the
% substrate, organism, UniProt accession and sequence id of each
% measurement. geckopy reads the same file (load_brenda_measurements).
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% folder : char
%     folder with measurements.tsv, instead of the adapter's BRENDA folder.
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

checkNameValue(varargin, {'modelAdapter','folder'}, 'loadBrendaMeasurements');
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
file = fullfile(folder, 'measurements.tsv');
if ~isfile(file)
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
