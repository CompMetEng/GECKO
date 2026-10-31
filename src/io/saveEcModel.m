function saveEcModel(ecModel,varargin)
% saveEcModel  Save an ecModel in YAML and/or SBML format.
%
% Saves the ecModel in either YAML format (= default and preferred, as all
% ecModel content is reserved) and/or SBML format (= more widely compatible
% with other constraint-based modelling tools, can be used for running
% simulations like FBA etc., but this model cannot be loaded back into
% MATLAB for applying further GECKO functions, as model.ec is lost).
%
% A YAML file records its ecModel format (gecko_format, see loadEcModel)
% and the GECKO release that wrote it (generator) in its metaData, in
% place of the geckopy_version that earlier geckopy releases wrote.
% A model without a description (ecModel.annotation.description) gets
% 'Enzyme-constrained model of <id>' in its metaData; an existing description
% is kept.
%
% The ecModel is saved to the ecModel-specific models/ subfolder. For
% saving to other locations, writeYAMLmodel or exportModel are more
% suitable.
%
% Parameters
% ----------
% ecModel : struct
%     an ecModel (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% filename : char
%     ending with either .yml or .xml, specifying if the ecModel should be
%     saved in YAML or SBML file format. If no file extension is given, the
%     ecModel will be saved in YAML format. If no filename is given,
%     'ecModel.yml' is used.
% modelAdapter : ModelAdapter
%     a loaded model adapter, from where the model folder is read (default:
%     the current default model adapter).
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     saveEcModel(ecModel,filename,modelAdapter);
%     saveEcModel(ecModel,'filename',filename);
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
    filename = 'ecModel';
end
filename = fullfile(params.path,'models',filename);

switch filename(end-3:end)
    case {'.xml','sbml'}
        exportModel(ecModel, filename);
        %Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
        %inconsistencies between Windows and MAC:
        copyfile(filename,'backup.xml')
        fin  = fopen('backup.xml', 'r');
        fout = fopen(filename, 'w');
        still_reading = true;
        while still_reading
            inline = fgets(fin);
            if ~ischar(inline)
                still_reading = false;
            else
                if ~isempty(regexp(inline,'-00[0-9]','once'))
                    inline = strrep(inline,'-00','-0');
                elseif ~isempty(regexp(inline,'-01[0-9]','once'))
                    inline = strrep(inline,'-01','-1');
                end
                fwrite(fout, inline);
            end
        end
        fclose('all');
        delete('backup.xml');
    otherwise % Assume YAML
        % The ecModel format the file follows (as project_format in
        % model_adapter.toml) and the release that wrote it.
        ecModel.annotation.gecko_format = '4';
        ecModel.annotation.generator = ['GECKO ' geckoVersion()];
        if ~isfield(ecModel.annotation, 'description') || isempty(ecModel.annotation.description)
            ecModel.annotation.description = ['Enzyme-constrained model of ' ecModel.id];
        end
        if isfield(ecModel.annotation, 'geckopy_version')
            ecModel.annotation = rmfield(ecModel.annotation, 'geckopy_version');
        end
        writeYAMLmodel(ecModel,filename);
end
end

function v = geckoVersion()
file = fullfile(findGECKOroot, 'version.txt');
if isfile(file)
    v = strtrim(fileread(file));
else
    v = 'unknown';
end
end
