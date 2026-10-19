function writtenTable = writeDLKcatInput(model, varargin)
% writeDLKcatInput  Prepare and write the input file for DLKcat.
%
% Prepares the input for DLKcat, and writes it to data/DLKcat.tsv in the
% obj.params.path specified in the ModelAdapter.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure).
%
% Name-Value Arguments
% --------------------
% ecRxns : logical
%     for which reactions (from model.ec.rxns) DLKcat should predict kcat
%     values, provided as logical vector with same length as model.ec.rxns
%     (default: all reactions).
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% onlyWithSmiles : logical
%     whether to only include metabolites with SMILES (default true).
% filename : char
%     path to the input file, including the filename and .tsv extension
%     (default: data/DLKcat.tsv from the obj.params.path folder specified in
%     the modelAdapter).
% overwrite : logical
%     whether existing file should be overwritten (default false, to prevent
%     overwriting a file that already contains DLKcat-predicted kcat values).
%
% Returns
% -------
% writtenTable : cell
%     the table written, mainly to be used for testing purposes.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     writtenTable = writeDLKcatInput(model);
%     writtenTable = writeDLKcatInput(model, 'onlyWithSmiles', false);
%
% See also
% --------
% runDLKcat, readDLKcatOutput, extractEnzymeSubstratePairs

p = parseGECKOargs(varargin, { ...
    'ecRxns',         []; ...
    'modelAdapter',   []; ...
    'onlyWithSmiles', []; ...
    'filename',       []; ...
    'overwrite',      []});
ecRxnsSel      = p.ecRxns;
modelAdapter   = p.modelAdapter;
onlyWithSmiles = p.onlyWithSmiles;
filename       = p.filename;
overwrite      = p.overwrite;

if isempty(ecRxnsSel)
    ecRxnsSel = true(numel(model.ec.rxns),1);
elseif ~logical(ecRxnsSel)
    error('ecRxns should be provided as logical vector')
elseif numel(ecRxnsSel)~=numel(model.ec.rxns)
    error('Length of ecRxns is not the same as model.ec.rxns')
end
ecRxnsSel = find(ecRxnsSel); % Change to indices, into the full model.ec.rxns

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;

if isempty(onlyWithSmiles)
    onlyWithSmiles=true;
end

if isempty(filename)
    filename = fullfile(params.path,'data','DLKcat.tsv');
elseif ~endsWith(filename,'.tsv')
    error('If filename is provided, it should include the .tsv extension.')
end

if isempty(overwrite) || ~overwrite % If is true
    if exist(filename,'file')
        error([filename ' already exists, either delete it first, or set the overwrite input argument as true'])
    end
end

pairs = extractEnzymeSubstratePairs(model, 'ecRxns', ismember((1:numel(model.ec.rxns))', ecRxnsSel), ...
    'modelAdapter', modelAdapter);

% Prepare output
out(1,:) = model.ec.rxns(pairs.ecRxn);
out(2,:) = model.ec.genes(pairs.protein);
out(3,:) = model.metNames(pairs.met);
if isfield(model,'metSmiles')
    out(4,:) = model.metSmiles(pairs.met);
else
    out(4,:) = cell(numel(pairs.met),1);
end

out(5,:) = model.ec.sequence(pairs.protein);
if onlyWithSmiles
    out(:,cellfun(@isempty,out(4,:))) = [];
else
    out(4,cellfun(@isempty,out(4,:))) = {'None'};
end
out(6,:) = cell(numel(out(1,:)),1);
out(6,:) = {'NA'};

% Write file
fID = fopen(filename,'w');
fprintf(fID,'%s\t%s\t%s\t%s\t%s\t%s\n',out{:});
fclose(fID);
fprintf('Model-specific DLKcat input stored at %s\n',filename);

writtenTable = out;
end
