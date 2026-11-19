function [evidence, summary] = kcatsForReactions(model, rxns, varargin)
% kcatsForReactions  Get the kcats of every source for a few reactions.
%
% Reads the kcats that the source files of the project give for the
% reactions, computes the EC-tree estimate for them, and optionally
% predicts their kcats with one or several OpenKineticsPredictor (OKP)
% methods. The result is shown next to the kcat of the model (see
% kcatEvidence).
%
% The OKP jobs have their own tag (see the jobTag option of kcatsFromOKP):
% they submit the enzyme-substrate pairs of the selected reactions only, and
% leave the job and the result that kcatsFromOKP keeps for all reactions as
% they are. The predictions are written to data/kcats/okp_<method>_query.tsv
% and read back as a source of kcatEvidence. The function waits for each OKP
% job to finish and submits a new job on every call. An OKP method that fails is
% reported and the other methods continue. Each OKP method submits a job and
% counts against the daily prediction quota of the API key.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure), with metabolite SMILES
%     (findMetSmiles) and protein sequences when OKP is used.
% rxns : cell, char, double or logical
%     the reactions, as in kcatEvidence: names of model.ec.rxns or
%     conventional reaction names (a reverse direction ends with _REV),
%     indices or a logical vector.
%
% Name-Value Arguments
% --------------------
% okp : logical, char or cell
%     OKP methods to predict with: false for none (default), true for the
%     method of the model adapter, 'all' for every method that kcatsFromOKP
%     accepts, or a method name or cell array of method names.
% sources, compute, measurements, modelAdapter
%     as in kcatEvidence; compute is {'ec_tree'} by default.
%
% Returns
% -------
% evidence, summary : table
%     as returned by kcatEvidence for the reactions.
%
% Examples
% --------
%     [~, summary] = kcatsForReactions(ecModel, 'MAR03905');
%     [~, summary] = kcatsForReactions(ecModel, {'MAR03905', 'MAR03905_REV'}, 'okp', 'all');
%     [~, summary] = kcatsForReactions(ecModel, 12, 'okp', {'KinForm-H', 'UniKP'});
%
% See also
% --------
% kcatEvidence, kcatsFromOKP, kcatsFromECTree

p = parseGECKOargs(varargin, { ...
    'okp',          false; ...
    'sources',      []; ...
    'compute',      {'ec_tree'}; ...
    'measurements', []; ...
    'modelAdapter', []});
modelAdapter = p.modelAdapter;
if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
[~, ~, rows] = kcatEvidence(model, 'rxns', rxns, 'sources', {}, 'modelAdapter', modelAdapter);
if isempty(rows)
    error('kcatsForReactions:none', 'No reaction of model.ec is selected.')
end

methods = okpMethods(p.okp, modelAdapter.params);
selected = false(numel(model.ec.rxns), 1);
selected(rows) = true;
queried = cell(0, 1);
for i = 1:numel(methods)
    try
        [~, done] = kcatsFromOKP(model, 'method', methods{i}, 'ecRxns', selected, ...
            'jobTag', 'query', 'resubmit', true, 'wait', true, 'modelAdapter', modelAdapter);
        if done
            queried{end+1, 1} = ['okp_' okpTag(methods{i}) '_query']; %#ok<AGROW>
        end
    catch err
        printOrange(sprintf('OKP %s failed: %s\n', methods{i}, err.message));
    end
end

sources = p.sources;
if isnumeric(sources) && isempty(sources)
    files = dir(fullfile(modelAdapter.params.path, 'data', 'kcats', '*.tsv'));
    sources = sort(cellfun(@(f) f(1:end-4), {files.name}, 'UniformOutput', false));
    % Predictions of an earlier query that were not asked for now are left out.
    old = startsWith(sources, 'okp_') & endsWith(sources, '_query') & ~ismember(sources, queried');
    sources = sources(~old);
end
[evidence, summary] = kcatEvidence(model, 'rxns', rows, 'sources', sources, 'compute', p.compute, ...
    'measurements', p.measurements, 'modelAdapter', modelAdapter);
end

function methods = okpMethods(okp, params)
% The OKP methods that okp asks for.
known = {'CataPro','CatRange','DLKcat','EITLEM','KinForm-H','KinForm-L', ...
    'OmniESI','OmniESI-O2DENet','UniKP'};
if islogical(okp) || isnumeric(okp)
    if ~okp
        methods = {};
    else
        config = kcatAssignmentConfig(params);
        methods = config.okpMethods(1);
    end
elseif (ischar(okp) || isstring(okp)) && strcmpi(okp, 'all')
    methods = known;
else
    methods = cellstr(okp);
    unknown = setdiff(methods, known);
    if ~isempty(unknown)
        error('kcatsForReactions:okp', 'Unknown OKP methods: %s. Methods: %s.', ...
            strjoin(unknown, ', '), strjoin(known, ', '))
    end
end
methods = methods(:)';
end

function tag = okpTag(method)
tag = lower(regexprep(regexprep(method, '[^A-Za-z0-9]+', '_'), '^_+|_+$', ''));
end
