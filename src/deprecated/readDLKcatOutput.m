function kcatList = readDLKcatOutput(model, varargin)
% readDLKcatOutput  Deprecated: use kcatsFromDLKcat(model, 'via', 'file').
%
% Reads a DLKcat output file into a kcatList through kcatsFromDLKcat, which
% reads data/DLKcat.tsv of the project. A file elsewhere (outFile) is read
% from a temporary copy. Will be removed in a future release.
%
% Name-Value Arguments
% --------------------
% outFile : char
%     the DLKcat output file (default: data/DLKcat.tsv of the project).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% strict : logical
%     ignored, with a warning: rows with a substrate that is not in
%     model.metNames are always left out with a warning.
%
% See also
% --------
% kcatsFromDLKcat

warning('GECKO:deprecatedFunction', ['readDLKcatOutput is deprecated; use ' ...
    'kcatsFromDLKcat(model, ''via'', ''file''). It will be removed in a future release.'])
p = parseGECKOargs(varargin, {'outFile', []; 'modelAdapter', []; 'strict', []});
if ~isempty(p.strict)
    warning('GECKO:deprecatedOption', ['The strict option of readDLKcatOutput is ignored: rows ' ...
        'with a substrate that is not in model.metNames are left out with a warning.'])
end
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
if ~isempty(p.outFile)
    project = tempname;
    mkdir(fullfile(project, 'data'));
    cleanup = onCleanup(@() rmdir(project, 's'));
    copyfile(p.outFile, fullfile(project, 'data', 'DLKcat.tsv'));
    modelAdapter.params.path = project;
end
kcatList = kcatsFromDLKcat(model, 'via', 'file', 'save', false, 'modelAdapter', modelAdapter);
end
