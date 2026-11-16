function model = loadConventionalGEM(varargin)
% loadConventionalGEM  Load the conventional GEM (non-ecModel).
%
% Loads the conventional GEM (non-ecModel) from the location specified in
% the modelAdapter. By default, it looks in the models/ subdirectory of the
% param.path specified in the modelAdapter. When loading conventional GEMs
% from other locations, one can directly use importModel.
%
% Name-Value Arguments
% --------------------
% filename : char
%     name of the model file, located in the models/ subfolder of
%     param.path as specified in the modelAdapter. A name ending in .xml is
%     loaded as SBML and a name ending in .yml as YAML; any other name
%     gets .xml appended (default: the value specified as param.convGEM in
%     the modelAdapter, which must end in .xml or .yml, or the error
%     loadConventionalGEM:extension is raised).
% modelAdapter : ModelAdapter
%     a loaded model adapter, from where the model folder is read (default:
%     the current default model adapter).
%
% Returns
% -------
% model : struct
%     model in RAVEN format.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     model = loadConventionalGEM(filename, modelAdapter);
%     model = loadConventionalGEM('filename', filename);
%
% See also
% --------
% loadEcModel

p = parseGECKOargs(varargin, { ...
    'filename',     []; ...
    'modelAdapter', []});
filename     = p.filename;
modelAdapter = p.modelAdapter;

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();
if isempty(filename)
    filename = params.convGEM;
else
    filename = char(filename);
    % A name like yeast-GEM-9.1.1 has dots, so only .xml and .yml count as extension
    if ~endsWith(filename, {'.xml', '.yml'})
        filename = [filename, '.xml'];
    end
    filename = fullfile(params.path,'models',filename);
end

[~, ~, ext] = fileparts(filename);
switch ext
    case '.xml'
        model = importModel(filename);
    case '.yml'
        model = readYAMLmodel(filename);
    otherwise
        error('loadConventionalGEM:extension', ...
            'The model file %s is neither SBML (.xml) nor YAML (.yml).', filename)
end
