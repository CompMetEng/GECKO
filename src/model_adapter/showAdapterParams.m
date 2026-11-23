function T = showAdapterParams(modelAdapter)
% showAdapterParams  Print the parameters of a model adapter and where each comes from.
%
% Prints one line per parameter of adapterParameters.tsv, in that order:
%
%     <key> = <value>  [<origin>]
%
% The origin is file (set in model_adapter.toml, or by an adapter class),
% default (the default of adapterParameters.tsv or of the code), preset
% <name> (given by the kcat assignment preset) or not set (no value; the
% line is "<key> =  [not set]"). Strings are shown bare, booleans as true or
% false, lists comma-joined, integer-valued numbers whole and other
% numbers with %g. The resolved kcat
% assignment is listed as kcat_assignment.stage[<i>] lines (the stages of the
% preset or the file, with the settings that the file sets or that differ
% from the defaults criteria max, top_origin_limit 4, bottom_origin_limit 6
% and wildcard_limit 3), the settings of each kcat source as
% kcat_sources.<name>.<key> lines and the evotune source groups, in group
% order, as evotune.source_groups.<group> lines with their sigma0_log.
% geckopy's "geckopy adapter show" prints the same lines.
%
% Parameters
% ----------
% modelAdapter : ModelAdapter, optional
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% T : table, optional
%     the listing, with the columns key, value and origin. Without an
%     output argument the lines are printed.
%
% Examples
% --------
%     showAdapterParams(ModelAdapterManager.getDefault())
%
% See also
% --------
% loadAdapterParams, kcatAssignmentConfig

if nargin < 1 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();
tomlFile = fullfile(params.path, 'model_adapter.toml');
hasFile = isfile(tomlFile);
data = struct();
if hasFile
    data = readAdapterToml(tomlFile);
end
config = kcatAssignmentConfig(params);

file = fullfile(fileparts(mfilename('fullpath')), 'adapterParameters.tsv');
fid = fopen(file, 'r');
cols = textscan(fid, '%s %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1, 'Whitespace', '');
fclose(fid);

key = {};
value = {};
origin = {};
    function add(k, v, o)
        key{end+1, 1} = k;
        value{end+1, 1} = v;
        origin{end+1, 1} = o;
    end

fromPreset = 'file';
if ~isempty(config.preset)
    fromPreset = ['preset ' config.preset];
end
for k = 1:numel(cols{1})
    tomlKey = cols{1}{k};
    field = cols{2}{k};
    if isempty(field)
        field = tomlKey;
    end
    switch tomlKey
        case 'evotune.source_groups'
            if ~isfield(params, 'evotune') || ~isfield(params.evotune, 'kcatGroups')
                add(tomlKey, '', 'not set');
                continue
            end
            entries = params.evotune.kcatSources;
            groups = params.evotune.kcatGroups;
            o = originOf(data, tomlKey, hasFile);
            [names, first] = unique(groups, 'stable');
            for j = 1:numel(names)
                in = strcmp(groups, names{j});
                add(['evotune.source_groups.' names{j}], sprintf('%s; sigma0_log %s', ...
                    strjoin(entries(in), ', '), formatValue(params.evotune.sigma0logSource(first(j)))), o);
            end
        case 'kcat_assignment.preset'
            if isempty(config.preset)
                add(tomlKey, 'none (stages are given)', 'file');
            else
                add(tomlKey, config.preset, originOf(data, tomlKey, hasFile));
            end
        case 'kcat_assignment.stage'
            given = struct();
            if isempty(config.preset)
                given = params.kcat_assignment.stage;
            end
            % A stage shows the settings that the file sets, or, for a preset,
            % those that are not the default.
            limits = {'criteria', 'criteria', 'max'; 'top_origin_limit', 'topOriginLimit', 4; ...
                'bottom_origin_limit', 'bottomOriginLimit', 6; 'wildcard_limit', 'wildcardLimit', 3};
            for j = 1:numel(config.stages)
                st = config.stages(j);
                text = sprintf('%s: %s', st.action, strjoin(st.sources, ', '));
                for i = 1:size(limits, 1)
                    if isempty(config.preset)
                        shown = isfield(given(j), limits{i,1}) && ~isempty(given(j).(limits{i,1}));
                    else
                        shown = ~isequal(st.(limits{i,2}), limits{i,3});
                    end
                    if shown
                        text = [text '; ' limits{i,1} '=' formatValue(st.(limits{i,2}))]; %#ok<AGROW>
                    end
                end
                add(sprintf('%s[%d]', tomlKey, j), text, fromPreset);
            end
        case 'kcat_sources'
            sources = config.sources;
            if ~(isfield(sources, 'okp') && isfield(sources.okp, 'method'))
                sources.okp.method = config.okpMethod;
            end
            % brenda, ec_tree, homologs and okp first, other sources and the
            % keys of each source alphabetically
            names = fieldnames(sources);
            first = {'brenda', 'ec_tree', 'homologs', 'okp'};
            names = [first(ismember(first, names)), sort(names(~ismember(names, first)))'];
            for j = 1:numel(names)
                keys = sort(fieldnames(sources.(names{j})));
                for i = 1:numel(keys)
                    given = isfield(params, 'kcat_sources') && isfield(params.kcat_sources, names{j}) ...
                        && isfield(params.kcat_sources.(names{j}), keys{i});
                    o = fromPreset;
                    if given
                        o = 'file';
                    elseif strcmp(names{j}, 'okp') && strcmp(keys{i}, 'method')
                        o = 'default';
                    end
                    add(sprintf('%s.%s.%s', tomlKey, names{j}, keys{i}), ...
                        formatValue(sources.(names{j}).(keys{i})), o);
                end
            end
        otherwise
            found = false;
            v = params;
            for part = strsplit(field, '.')
                if isstruct(v) && isfield(v, part{1})
                    v = v.(part{1});
                    found = true;
                else
                    found = false;
                    break
                end
            end
            if ~found || (isnumeric(v) && isempty(v))
                add(tomlKey, '', 'not set');
                continue
            end
            o = originOf(data, tomlKey, hasFile);
            add(tomlKey, formatValue(v), o);
    end
end

T = table(key, value, origin, 'VariableNames', {'key', 'value', 'origin'});
if nargout == 0
    for k = 1:height(T)
        if strcmp(origin{k}, 'not set')
            fprintf('%s =  [not set]\n', key{k});
        else
            fprintf('%s = %s  [%s]\n', key{k}, value{k}, origin{k});
        end
    end
    clear T
end
end

function o = originOf(data, tomlKey, hasFile)
% file when the TOML sets the key, default otherwise, file for an adapter
% class without a TOML file (it sets the value).
inFile = true;
w = data;
for part = strsplit(tomlKey, '.')
    inFile = isstruct(w) && isfield(w, part{1});
    if ~inFile
        break
    end
    w = w.(part{1});
end
if inFile || ~hasFile
    o = 'file';
else
    o = 'default';
end
end

function text = formatValue(v)
if ischar(v) || isstring(v)
    text = char(v);
elseif islogical(v)
    text = char(string(v));
elseif isnumeric(v)
    % integer-valued numbers whole, others with %g
    text = strjoin(arrayfun(@(x) numberText(x), v(:)', 'UniformOutput', false), ', ');
elseif iscell(v)
    text = strjoin(cellfun(@formatValue, v(:)', 'UniformOutput', false), ', ');
elseif isstruct(v)
    text = sprintf('struct with fields %s', strjoin(fieldnames(v), ', '));
else
    text = class(v);
end
end

function text = numberText(x)
if isfinite(x) && x == round(x)
    text = sprintf('%.0f', x);
else
    text = sprintf('%g', x);
end
end
