function params = loadAdapterParams(tomlFile)
% loadAdapterParams  Read model adapter parameters from a TOML file.
%
% Reads model_adapter.toml, the configuration file shared by GECKO and
% geckopy, and returns the parameters with GECKO's field names, so that
% they can be used as ModelAdapter.params. Keys are translated with
% adapterParameters.tsv (for instance p_tot becomes Ptot); keys that are
% not given get the default listed there, which is geckopy's default.
% A key that is not known is an error, also inside [evotune.source_groups],
% [[kcat_assignment.stage]] and [kcat_sources.<name>]; the keys of a kcat
% source are those that updateKcatSources reads.
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
checkNested(data, tomlFile);

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
    checkMatlabParams(data.matlab.params, tomlFile);
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
    unknownKey(key, tomlFile);
end
end

function checkMatlabParams(extra, tomlFile)
% [matlab.params] is for parameters that only GECKO uses; the settings that
% have a key of their own are not set there.
fields = fieldnames(extra);
for k = 1:numel(fields)
    if strcmp(fields{k}, 'kcatAggregation')
        where = '[kcat_sources.brenda] aggregation';
    elseif strcmp(fields{k}, 'evotune')
        where = '[evotune]';
    else
        continue
    end
    error('GECKO:loadAdapterParams', '%s: matlab.params.%s is not allowed; set it in %s', ...
        tomlFile, fields{k}, where);
end
end

function unknownKey(key, tomlFile)
error('GECKO:loadAdapterParams', '%s: unknown key %s in model_adapter.toml', tomlFile, key);
end

function checkNested(data, tomlFile)
% Keys inside the tables that adapterParameters.tsv lists as one key.
if isfield(data, 'evotune') && isfield(data.evotune, 'source_groups')
    groups = data.evotune.source_groups;
    names = fieldnames(groups);
    for k = 1:numel(names)
        g = groups.(names{k});
        prefix = ['evotune.source_groups.' names{k}];
        if ~isstruct(g)
            error('GECKO:loadAdapterParams', ...
                '%s: %s must be a table with sources and sigma0_log', tomlFile, prefix);
        end
        checkKeys(g, {'sources', 'sigma0_log'}, prefix, tomlFile);
        if ~isfield(g, 'sources') || ~isfield(g, 'sigma0_log')
            error('GECKO:loadAdapterParams', ...
                '%s: %s needs both sources and sigma0_log', tomlFile, prefix);
        end
        if ~(iscell(g.sources) && ~isempty(g.sources) && all(cellfun(@ischar, g.sources)))
            error('GECKO:loadAdapterParams', ...
                '%s: %s.sources must be a non-empty list of strings', tomlFile, prefix);
        end
        w = g.sigma0_log;
        if ~(isnumeric(w) && isscalar(w) && isfinite(w) && w > 0)
            error('GECKO:loadAdapterParams', ...
                '%s: %s.sigma0_log must be a positive number', tomlFile, prefix);
        end
    end
end
if isfield(data, 'kcat_assignment') && isfield(data.kcat_assignment, 'stage')
    stages = data.kcat_assignment.stage;
    for k = 1:numel(stages)
        checkKeys(stages(k), {'action', 'sources', 'criteria', 'top_origin_limit', ...
            'bottom_origin_limit', 'wildcard_limit'}, 'kcat_assignment.stage', tomlFile);
    end
end
if isfield(data, 'kcat_sources') && isstruct(data.kcat_sources)
    % The sources that GECKO knows and the keys that updateKcatSources
    % reads for each; dlkcat, custom and the okp_<method> files take none.
    allowed = struct('brenda', {{'aggregation', 'force_wildcard_level'}}, ...
        'ec_tree', {{'k'}}, 'homologs', {{'k', 'min_identity', 'coverage'}}, ...
        'okp', {{'methods'}}, 'dlkcat', {{}}, 'custom', {{}});
    names = fieldnames(data.kcat_sources);
    for k = 1:numel(names)
        if isfield(allowed, names{k})
            keys = allowed.(names{k});
        elseif any(strcmp(names{k}, okpFileNames()))
            keys = {};
        else
            unknownKey(['kcat_sources.' names{k}], tomlFile);
        end
        source = data.kcat_sources.(names{k});
        if ~isstruct(source)
            error('GECKO:loadAdapterParams', '%s: kcat_sources.%s must be a table', tomlFile, names{k});
        end
        checkKeys(source, keys, ['kcat_sources.' names{k}], tomlFile);
        if strcmp(names{k}, 'brenda') && isfield(source, 'aggregation') && ...
                ~(ischar(source.aggregation) && any(strcmp(source.aggregation, {'max', 'median'})))
            error('GECKO:loadAdapterParams', ...
                '%s: kcat_sources.brenda.aggregation must be max or median', tomlFile);
        end
        if strcmp(names{k}, 'okp') && isfield(source, 'methods')
            checkOkpMethods(source.methods, tomlFile);
        end
    end
end
end

function checkOkpMethods(methods, tomlFile)
% [kcat_sources.okp] methods: a list of the OKP methods that kcatsFromOKP accepts.
if ~(iscell(methods) && all(cellfun(@ischar, methods)))
    error('GECKO:loadAdapterParams', '%s: kcat_sources.okp.methods must be a list of strings', tomlFile);
end
known = okpMethodNames();
unknown = unique(methods(~ismember(methods, known)), 'stable');
if ~isempty(unknown)
    error('kcatAssignment:config', 'Unknown OKP methods in [kcat_sources.okp] methods: %s; methods: %s.', ...
        strjoin(unknown, ', '), strjoin(known, ', '))
end
end

function known = okpMethodNames()
% The methods that predict from one enzyme and one substrate (CatPred and
% TurNup read whole reactions).
known = {'CataPro','CatRange','DLKcat','EITLEM','KinForm-H','KinForm-L', ...
    'OmniESI','OmniESI-O2DENet','UniKP'};
end

function files = okpFileNames()
% The files of the OKP methods that kcatsFromOKP accepts: okp_ and the
% lower-case method name, with runs of characters other than letters and
% digits replaced by _.
known = okpMethodNames();
files = cellfun(@(m) ['okp_' lower(regexprep(regexprep(m, '[^A-Za-z0-9]+', '_'), '^_+|_+$', ''))], ...
    known, 'UniformOutput', false);
end

function checkKeys(s, allowed, prefix, tomlFile)
fields = fieldnames(s);
for k = 1:numel(fields)
    if ~any(strcmp(fields{k}, allowed))
        unknownKey([prefix '.' fields{k}], tomlFile);
    end
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
