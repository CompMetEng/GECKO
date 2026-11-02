function folder = fetchBrendaAsset(files)
% fetchBrendaAsset  Download (once) GECKO's BRENDA files and return their folder.
%
% GECKO does not ship BRENDA data. The files of the toolbox-assets release
% brenda-2026.1-1 (BRENDA 2026.1, CC BY 4.0) are downloaded on first use
% with RAVEN's fetchToolboxAsset, which checks them against the checksums
% of the toolbox-assets manifest, and are unpacked into
% databases/brenda-2026.1-1 in the GECKO folder. geckopy reads the same
% release (fetch_brenda_asset). A project can use its own BRENDA files
% instead, in data/brenda (see ModelAdapter.getBrendaDBFolder).
%
% Parameters
% ----------
% files : cell
%     file names (default: kcat.tsv, sa.tsv, mw.tsv and measurements.tsv);
%     the release also has sequences.fasta, substrates.tsv and
%     kcatCombineParams.tsv.
%
% Returns
% -------
% folder : char
%     the folder with the files, unpacked.
%
% Examples
% --------
%     folder = fetchBrendaAsset();
%     folder = fetchBrendaAsset({'kcatCombineParams.tsv'});
%
% See also
% --------
% loadBRENDAdata, loadBrendaMeasurements, ModelAdapter

release = 'brenda-2026.1-1';
if nargin < 1 || isempty(files)
    files = {'kcat.tsv', 'sa.tsv', 'mw.tsv', 'measurements.tsv'};
end
files = cellstr(files);
packed = struct('kcat_tsv', 1, 'sa_tsv', 1, 'mw_tsv', 1, 'measurements_tsv', 1, ...
    'sequences_fasta', 1, 'substrates_tsv', 1, 'kcatCombineParams_tsv', 0);
folder = fullfile(findGECKOroot(), 'databases', release);
if exist('fetchToolboxAsset', 'file') ~= 2
    error('fetchBrendaAsset:raven', ['GECKO downloads its BRENDA files with RAVEN''s ' ...
        'fetchToolboxAsset, which this RAVEN does not have. Update RAVEN to version 3.'])
end
for i = 1:numel(files)
    name = files{i};
    key = matlab.lang.makeValidName(name);
    if ~isfield(packed, key)
        error('fetchBrendaAsset:file', '%s is not a file of %s.', name, release)
    end
    target = fullfile(folder, name);
    if isfile(target)
        continue
    end
    if packed.(key)
        gz = fetchToolboxAsset(folder, release, [name '.gz']);
        gunzip(gz, folder);
    else
        fetchToolboxAsset(folder, release, name);
    end
end
end
