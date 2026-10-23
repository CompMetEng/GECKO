function model = loadEcModel(varargin)
% loadEcModel  Load the ecModel that matches the modelAdapter.
%
% Loads the ecModel that matches the modelAdapter. By default, it loads the
% models/ecModel.yml in the directory specified as param.path in the
% modelAdapter. Alternative files in the same folder can be loaded by
% providing the appropriate filename. If loading models from other
% locations, one can directly use readYAMLmodel.
%
% An ecModel whose protein usage and pool reactions run in the GECKO 3
% direction is upgraded in the loaded model with upgradeEcModel, with a
% warning. An ecModel format (gecko_format in the file's metaData, see
% saveEcModel) above 4 is an error.
%
% Name-Value Arguments
% --------------------
% filename : char
%     name of the ecModel file (default 'ecModel.yml'). File should be
%     located in the models/ subfolder of param.path as specified in the
%     modelAdapter.
% modelAdapter : ModelAdapter
%     a loaded model adapter, from where the model folder is read (default:
%     the current default model adapter).
%
% Returns
% -------
% model : struct
%     an ecModel (with ecModel.ec structure).
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     model = loadEcModel(filename, modelAdapter);
%     model = loadEcModel('filename', filename);
%
% See also
% --------
% loadConventionalGEM, saveEcModel

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
    filename = 'ecModel.yml';
end

if ~endsWith(filename,{'yml','yaml'})
    error(['ecModels should preferably be distributed in YAML file format, ', ...
           'as otherwise model content will be lost. If the model is in ', ...
           'SBML format, you can use the usual importModel function to load ', ...
           'this as any other genome-scale model.'])
end

filename = fullfile(params.path,'models',filename); 
if endsWith(filename,{'yml','yaml'})
    model = readYAMLmodel(filename);
elseif endsWith(filename,{'xml','sbml'})
    model = importModel(filename);
end

if isfield(model, 'annotation') && isfield(model.annotation, 'gecko_format') ...
        && str2double(model.annotation.gecko_format) > 4
    writer = 'an unknown release';
    if isfield(model.annotation, 'generator')
        writer = model.annotation.generator;
    end
    error('loadEcModel:format', ['%s has ecModel format %s (written by %s); ' ...
        'this GECKO reads format 4 and older. Update GECKO to load it.'], ...
        filename, model.annotation.gecko_format, writer)
end
[model, changed] = upgradeEcModel(model);
if ~isempty(changed)
    warning('loadEcModel:upgraded', ['%s has %d protein usage and pool reactions in the GECKO 3 ' ...
        'direction; they are turned forward in the loaded model. upgradeGECKOproject upgrades the ' ...
        'files of the project.'], filename, numel(changed))
end
end
