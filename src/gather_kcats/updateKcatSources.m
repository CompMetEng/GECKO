function report = updateKcatSources(model, varargin)
% updateKcatSources  Gather the kcat sources that the configuration needs.
%
% Lists the kcat source files that [kcat_assignment] in model_adapter.toml
% uses (see collectKcats), and gathers those that are missing or outdated
% with kcatsFromBrenda, kcatsFromECTree, kcatsFromDLKcat and kcatsFromOKP:
%
% - a file is outdated when its header records other settings than
%   [kcat_sources.<name>] asks for (presets set these for brenda and
%   ec_tree), or, for brenda and ec_tree, another BRENDA release than the
%   BRENDA files in use;
% - brenda and ec_tree are made again; dlkcat reads data/DLKcat.tsv if it exists and
%   runs DLKcat (Docker) otherwise; okp_<method> submits or fetches a job,
%   and reports it as pending while the job runs;
% - custom.tsv and other files it cannot gather are left as they are.
%
% Requirements are checked first: Docker for a DLKcat run and an
% OpenKineticsPredictor API key for OKP. Calling it again is safe: files
% that are current are not gathered again.
%
% Per-source settings are read from model_adapter.toml, for instance:
%
%     [kcat_sources.brenda]
%     aggregation = "median"
%
%     [kcat_sources.okp_catapro]
%     method = "CataPro"
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% dryRun : logical
%     only report what would be gathered (default false).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
%
% Returns
% -------
% report : table
%     one row per source: source, status (current, gathered, pending,
%     outdated, missing, failed or skipped) and note.
%
% Examples
% --------
%     updateKcatSources(ecModel, 'dryRun', true)
%     updateKcatSources(ecModel);
%
% See also
% --------
% collectKcats, assignKcats, kcatsFromBrenda, kcatsFromDLKcat, kcatsFromOKP

checkOptions(varargin, {'dryRun','modelAdapter'}, 'updateKcatSources');
p = parseGECKOargs(varargin, {'dryRun', false; 'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
config = kcatAssignmentConfig(params);
kcatDir = fullfile(params.path, 'data', 'kcats');

names = {};
for s = 1:numel(config.stages)
    for t = 1:numel(config.stages(s).files)
        token = config.stages(s).files{t};
        if startsWith(token, 'database_')
            token = 'brenda';
        end
        if ~any(strcmp(names, token))
            names{end+1} = token; %#ok<AGROW>
        end
    end
end

status = cell(numel(names), 1);
note = cell(numel(names), 1);
for k = 1:numel(names)
    name = names{k};
    file = fullfile(kcatDir, [name '.tsv']);
    wanted = struct();
    if isfield(config.sources, name)
        wanted = config.sources.(name);
    end
    [status{k}, note{k}] = check(name, file, wanted);
    if any(strcmp(status{k}, {'missing','outdated'})) && ~p.dryRun
        [status{k}, note{k}] = gather(model, name, wanted, modelAdapter);
    end
end
report = table(names(:), status, note, 'VariableNames', {'source','status','note'});
disp(report)

    function [st, nt] = check(name, file, wanted)
        st = 'current'; nt = '';
        if strcmp(name, 'custom')
            st = 'skipped'; nt = 'curated by hand';
            return
        elseif ~any(strcmp(name, {'brenda','ec_tree','dlkcat'})) && ~startsWith(name, 'okp_')
            st = 'skipped'; nt = 'not a source GECKO gathers';
            return
        end
        if ~isfile(file)
            st = 'missing';
            return
        end
        meta = readHeaderLines(file);
        recorded = struct();
        if isfield(meta, 'parameters') && ~isempty(meta.parameters)
            recorded = jsondecode(meta.parameters);
        end
        f = fieldnames(wanted);
        for i = 1:numel(f)
            if ~isfield(recorded, f{i}) || ~isequal(normalizeValue(recorded.(f{i})), normalizeValue(wanted.(f{i})))
                st = 'outdated';
                nt = sprintf('%s differs from model_adapter.toml', f{i});
                return
            end
        end
        if any(strcmp(name, {'brenda','ec_tree'}))
            release = brendaRelease(modelAdapter, name);
            if isfield(meta, 'data') && ~isempty(release) && ~strcmp(meta.data, release)
                st = 'outdated';
                nt = sprintf('made from %s; GECKO has %s', meta.data, release);
            end
        end
    end
end

function [st, nt] = gather(model, name, wanted, modelAdapter)
st = 'gathered'; nt = '';
try
    if strcmp(name, 'brenda')
        args = {};
        if isfield(wanted, 'aggregation'), args = [args, {'aggregation', wanted.aggregation}]; end
        if isfield(wanted, 'force_wildcard_level'), args = [args, {'forceWClvl', wanted.force_wildcard_level}]; end
        kcatsFromBrenda(model, args{:}, 'modelAdapter', modelAdapter);
    elseif strcmp(name, 'ec_tree')
        args = {};
        if isfield(wanted, 'k'), args = {'k', double(wanted.k)}; end
        kcatsFromECTree(model, args{:}, 'modelAdapter', modelAdapter);
    elseif strcmp(name, 'dlkcat')
        run = ~isfile(fullfile(modelAdapter.params.path, 'data', 'DLKcat.tsv'));
        if run
            [status, ~] = system('docker --version');
            if status ~= 0
                st = 'failed'; nt = 'DLKcat needs Docker to run';
                return
            end
        end
        kcatsFromDLKcat(model, 'run', run, 'modelAdapter', modelAdapter);
    else % okp_<method>
        method = okpMethod(name, wanted);
        resolveOkpApiKey('', modelAdapter.params.path); % errors without a key
        [~, done] = kcatsFromOKP(model, 'method', method, 'saveAs', name, 'modelAdapter', modelAdapter);
        if ~done
            st = 'pending'; nt = sprintf('OKP job for %s still running; call again later', method);
        end
    end
catch err
    st = 'failed'; nt = err.message;
end
end

function method = okpMethod(name, wanted)
% The OKP method of okp_<method>: from [kcat_sources.okp_<method>] method,
% or the known method whose name gives this file name.
if isfield(wanted, 'method')
    method = char(wanted.method);
    return
end
known = {'CataPro','CatPred','CatRange','DLKcat','EITLEM','KinForm-H','KinForm-L', ...
    'OmniESI','OmniESI-O2DENet','TurNuP','UniKP'};
tags = cellfun(@(m) ['okp_' lower(regexprep(m, '[^A-Za-z0-9]+', '_'))], known, 'UniformOutput', false);
hit = strcmp(tags, name);
if ~any(hit)
    error('updateKcatSources:okpMethod', ['Set the OKP method of %s in model_adapter.toml: ' ...
        '[kcat_sources.%s] method = "..."'], name, name)
end
method = known{hit};
end

function v = normalizeValue(v)
if isstring(v)
    v = char(v);
end
if ischar(v)
    v = lower(v);
elseif islogical(v)
    v = double(v);
end
end

function release = brendaRelease(modelAdapter, name)
% The release line of the BRENDA file a source is made from.
release = '';
files = struct('brenda', 'kcat.tsv', 'ec_tree', 'measurements.tsv');
fid = fopen(fullfile(modelAdapter.getBrendaDBFolder(), files.(name)), 'r');
if fid < 0
    return
end
line = fgetl(fid);
fclose(fid);
tok = regexp(line, '^#\s*(BRENDA release [^ ]+)', 'tokens', 'once');
if ~isempty(tok)
    release = tok{1};
end
end

function meta = readHeaderLines(filename)
% Only the header lines at the top, so that a damaged table does not stop
% a gathered file from being written again.
meta = struct();
fid = fopen(filename, 'r', 'n', 'UTF-8');
if fid < 0
    return
end
cleanup = onCleanup(@() fclose(fid));
line = fgetl(fid);
if ischar(line)
    line = regexprep(line, ['^' char(65279)], '');
end
while ischar(line) && (isempty(strtrim(line)) || startsWith(unquoteField(line), '#'))
    [key, value] = parseHeaderLine(line);
    if ~isempty(key)
        meta.(key) = value;
    end
    line = fgetl(fid);
end
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
