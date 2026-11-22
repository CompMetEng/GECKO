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
% they are. The predictions are written to data/kcats/okp_<method>__query.tsv
% and read back as a source of kcatEvidence. The function waits for each OKP
% job to finish and submits a new job on every call. An OKP method that fails is
% reported and the other methods continue; the methods without predictions are
% listed at the end. Each OKP method submits a job and
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
%     OKP methods to predict with: false or empty for none (default), true for the
%     method of the model adapter, 'all' for every method that kcatsFromOKP
%     accepts, or a method name or cell array of method names.
% sources, compute, measurements, modelAdapter
%     as in kcatEvidence; compute is {'ec_tree'} by default, so the EC-tree
%     estimate is computed from the BRENDA measurements, also when the project
%     has an ec_tree.tsv. With sources, the files of the OKP methods that were
%     asked for are added to it. The default sources include the okp_<method>__query.tsv
%     files of earlier queries, which show the predictions of an earlier query for the same
%     reactions.
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

checkOptions(varargin, {'okp','sources','compute','measurements','modelAdapter'}, 'kcatsForReactions');
p = parseRAVENargs(varargin, { ...
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
queried = cell(1, 0);
failed = cell(1, 0);
for i = 1:numel(methods)
    try
        [~, done] = kcatsFromOKP(model, 'method', methods{i}, 'ecRxns', selected, ...
            'jobTag', 'query', 'resubmit', true, 'wait', true, 'modelAdapter', modelAdapter);
    catch err
        printOrange(sprintf('OKP %s failed: %s\n', methods{i}, err.message));
        done = false;
    end
    if done
        queried{end+1} = ['okp_' okpTag(methods{i}) '__query']; %#ok<AGROW>
    else
        failed{end+1} = methods{i}; %#ok<AGROW>
    end
end
if ~isempty(failed)
    printOrange(sprintf('No OKP predictions from: %s.\n', strjoin(failed, ', ')));
end

sources = p.sources;
if ~(isnumeric(sources) && isempty(sources))
    sources = [cellstr(sources(:))', queried];
end
[evidence, summary] = kcatEvidence(model, 'rxns', rows, 'sources', sources, 'compute', p.compute, ...
    'measurements', p.measurements, 'modelAdapter', modelAdapter);
end

function methods = okpMethods(okp, params)
% The OKP methods that okp asks for.
known = {'CataPro','CatRange','DLKcat','EITLEM','KinForm-H','KinForm-L', ...
    'OmniESI','OmniESI-O2DENet','UniKP'};
if isempty(okp) || islogical(okp) || isnumeric(okp)
    if isempty(okp) || ~okp
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
