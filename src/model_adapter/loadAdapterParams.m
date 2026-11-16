function params = loadAdapterParams(tomlFile)
% loadAdapterParams  Read model adapter parameters from a TOML file.
%
% Reads model_adapter.toml, the configuration file shared by GECKO and
% geckopy, and returns the parameters with GECKO's field names, so that
% they can be used as ModelAdapter.params. Keys are translated with
% adapterParameters.tsv (for instance p_tot becomes Ptot); keys that are
% not given get the default listed there, which is geckopy's default.
%
% Parameters
% ----------
% tomlFile : char
%     path to model_adapter.toml, or to the project folder that holds it.
%
% Returns
% -------
% params : struct
%     the adapter parameters, as in ModelAdapter.params:
%     - path: the project folder (the folder of the TOML file);
%     - convGEM and spontaneous.file: absolute paths;
%     - evotune.kcatSources and evotune.sigma0logSource: built from the
%       evotune.source_groups tables, one entry per element of sources
%       (a prefix pattern keeps its trailing *), in group order, with the
%       sigma0_log of its group, and evotune.kcatGroups the group name of
%       each entry. Without [evotune.source_groups] the
%       groups are custom 0.1, brenda 0.2, ec_tree 0.3 and combined 0.4;
%     - fields under [matlab.params] are copied as they are, for
%       parameters that only GECKO uses.
%
% Examples
% --------
%     params = loadAdapterParams(fullfile(projectPath,'model_adapter.toml'));
%
% See also
% --------
% TomlAdapter, readAdapterToml, convertAdapterToToml

if isfolder(tomlFile)
    tomlFile = fullfile(tomlFile, 'model_adapter.toml');
end
if ~isfile(tomlFile)
    error('GECKO:loadAdapterParams', 'No model adapter file found at %s', tomlFile);
end
tomlFile = char(tomlFile);
[projectPath, ~, ~] = fileparts(getAbsolutePath(tomlFile));

data = readAdapterToml(tomlFile);
map = adapterParameterMap();

% Unknown keys are errors, as in geckopy, so that a typo is not ignored.
known = map.tomlKey;
checkKnownKeys(data, '', known, tomlFile);

params = struct();
params.path = projectPath;
for k = 1:numel(map.tomlKey)
    key = map.tomlKey{k};
    if strcmp(key, 'evotune.source_groups')
        continue
    end
    [found, value] = getDotted(data, key);
    if ~found
        if ~map.hasDefault(k)
            continue
        end
        value = map.default{k};
    end
    if isempty(map.matlabField{k})
        field = key; % a key that only geckopy reads keeps its name
    else
        field = map.matlabField{k};
    end
    params = setDotted(params, field, value);
end

if ~isfield(params, 'complex')
    params.complex.taxonomicID = []; % geckopy's default: no Complex Portal query
end

for f = {'convGEM'}
    if isfield(params, f{1})
        params.(f{1}) = projectFile(projectPath, params.(f{1}));
    end
end
if isfield(params, 'spontaneous') && isfield(params.spontaneous, 'file')
    params.spontaneous.file = projectFile(projectPath, params.spontaneous.file);
end

params = sourceGroups(params, data);

if isfield(data, 'matlab') && isfield(data.matlab, 'params')
    params = mergeStruct(params, data.matlab.params);
end
end

function map = adapterParameterMap()
% Reads adapterParameters.tsv into toml keys, MATLAB fields and defaults.
file = fullfile(fileparts(mfilename('fullpath')), 'adapterParameters.tsv');
fid = fopen(file, 'r');
cols = textscan(fid, '%s %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1, ...
    'Whitespace', '', 'ReturnOnError', false);
fclose(fid);
map.tomlKey = cols{1};
map.matlabField = cols{2};
map.default = cell(size(cols{3}));
% A key without a default is left out of params; a default of "" is kept.
map.hasDefault = ~cellfun(@isempty, cols{3});
for k = 1:numel(cols{3})
    if map.hasDefault(k)
        map.default{k} = parseDefault(cols{3}{k});
    end
end
end

function value = parseDefault(text)
% A default is written as a TOML value.
tmp = [tempname '.toml'];
fid = fopen(tmp, 'w');
fprintf(fid, 'x = %s\n', text);
fclose(fid);
cleanup = onCleanup(@() delete(tmp));
data = readAdapterToml(tmp);
value = data.x;
if ischar(value) && isempty(value)
    value = '';
end
end

function checkKnownKeys(s, prefix, known, tomlFile)
fields = fieldnames(s);
for k = 1:numel(fields)
    key = fields{k};
    if ~isempty(prefix)
        key = [prefix '.' key];
    end
    if strcmp(key, 'matlab.params')
        continue % free-form parameters for GECKO only
    end
    if any(strcmp(key, known)) || any(startsWith(known, [key '.']))
        if isstruct(s.(fields{k})) && ~any(strcmp(key, known))
            checkKnownKeys(s.(fields{k}), key, known, tomlFile);
        end
        continue
    end
    error('GECKO:loadAdapterParams', ['Unknown key "%s" in %s. Parameters that only ' ...
        'GECKO uses go under [matlab.params].'], key, tomlFile);
end
end

function params = sourceGroups(params, data)
% evotune.source_groups becomes the kcatSources and sigma0logSource lists
% that evotune reads: one entry per element of sources, in group order,
% with the sigma0_log of its group.
groups = struct();
if isfield(data, 'evotune') && isfield(data.evotune, 'source_groups')
    groups = data.evotune.source_groups;
else
    defaults = {'custom', 0.1; 'brenda', 0.2; 'ec_tree', 0.3; 'combined', 0.4};
    for k = 1:size(defaults, 1)
        groups.(defaults{k,1}).sources = defaults(k,1);
        groups.(defaults{k,1}).sigma0_log = defaults{k,2};
    end
end
names = fieldnames(groups);
kcatSources = {};
kcatGroups = {};
sigma0logSource = [];
for k = 1:numel(names)
    sources = groups.(names{k}).sources(:)';
    kcatSources = [kcatSources, sources]; %#ok<AGROW>
    kcatGroups = [kcatGroups, repmat(names(k), 1, numel(sources))]; %#ok<AGROW>
    sigma0logSource = [sigma0logSource; repmat(groups.(names{k}).sigma0_log, numel(sources), 1)]; %#ok<AGROW>
end
params.evotune.kcatSources = kcatSources;
params.evotune.kcatGroups = kcatGroups;
params.evotune.sigma0logSource = sigma0logSource;
end

function [found, value] = getDotted(s, key)
parts = strsplit(key, '.');
value = s;
found = true;
for k = 1:numel(parts)
    if ~isstruct(value) || ~isfield(value, parts{k})
        found = false;
        value = [];
        return
    end
    value = value.(parts{k});
end
end

function s = setDotted(s, key, value)
parts = strsplit(key, '.');
if numel(parts) == 1
    s.(parts{1}) = value;
    return
end
if ~isfield(s, parts{1}) || ~isstruct(s.(parts{1}))
    s.(parts{1}) = struct();
end
s.(parts{1}) = setDotted(s.(parts{1}), strjoin(parts(2:end), '.'), value);
end

function s = mergeStruct(s, extra)
fields = fieldnames(extra);
for k = 1:numel(fields)
    value = extra.(fields{k});
    if isstruct(value) && isscalar(value) && isfield(s, fields{k}) && isstruct(s.(fields{k}))
        s.(fields{k}) = mergeStruct(s.(fields{k}), value);
    else
        s.(fields{k}) = value;
    end
end
end

function p = projectFile(projectPath, p)
if isempty(p)
    return
end
if ~isAbsolute(p)
    p = fullfile(projectPath, p);
end
end

function tf = isAbsolute(p)
tf = startsWith(p, '/') || startsWith(p, '\') || ~isempty(regexp(p, '^[A-Za-z]:', 'once'));
end

function p = getAbsolutePath(p)
if ~isAbsolute(p)
    p = fullfile(pwd, p);
end
end
