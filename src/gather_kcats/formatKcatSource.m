function text = formatKcatSource(origin, detail, tuned, prior)
% formatKcatSource  Write the provenance text of a kcat in model.ec.source.
%
% The text has the form
%
%     <origin>[ (<detail>)][ | tuned=<method>[+<method>...] (prior=<value>)]
%
% for instance 'brenda (wc=0, origin=1)', 'combined (sd=0.61, n=6)' or
% 'brenda (wc=0, origin=2) | tuned=cmaes (prior=6.55)'. The origin and the
% detail keys are written in lower case, with every run of characters other
% than letters and digits replaced by one underscore, so 'OKP-CataPro' gives
% 'okp_catapro'. parseKcatSource reads the text back.
%
% Parameters
% ----------
% origin : char
%     where the kcat comes from, e.g. 'brenda', 'custom' or 'standard'. An
%     empty origin gives '' when there is no detail and no tuning status,
%     and 'unknown' otherwise.
% detail : struct
%     scalar struct whose fields are the detail keys, in output order, e.g.
%     struct('wc', 0, 'origin', 1). A number is written with '%.6g' and
%     text as it is, so a fixed number of decimals is given as text, e.g.
%     sprintf('%.2f', sd). A field whose value is NaN, Inf or [] is left out.
%     A field is also left out when its name does not start with a letter
%     after the normalisation. The brackets are left out when no field is
%     left (default []).
% tuned : cell or char
%     the tuning methods that changed the kcat, in the order of application,
%     e.g. {'sensitivity','cmaes'} (default {}).
% prior : double
%     the kcat before the first tuning, written when tuned is not empty and
%     the prior is a finite number; NaN, Inf or [] writes none (default []).
%
% Returns
% -------
% text : char
%     the provenance text.
%
% Examples
% --------
%     formatKcatSource('standard')
%     formatKcatSource('brenda', struct('wc', 0, 'origin', 1))
%     formatKcatSource('combined', struct('sd', sprintf('%.2f', 0.61), 'n', 6))
%     formatKcatSource('brenda', [], {'cmaes'}, 6.55)
%
% See also
% --------
% parseKcatSource, assignKcatValues, sensitivityTuning
if nargin < 2
    detail = [];
end
if nargin < 3 || isempty(tuned)
    tuned = {};
end
if nargin < 4
    prior = [];
end
tuned = cellstr(tuned);

pairs = {};
if isstruct(detail)
    keys = fieldnames(detail);
    for i = 1:numel(keys)
        value = detail.(keys{i});
        key = normaliseToken(keys{i});
        if isempty(regexp(key, '^[a-z][a-z0-9_]*$', 'once'))
            continue
        end
        if isnumeric(value) || islogical(value)
            if isempty(value) || ~isfinite(value(1))
                continue
            end
            value = sprintf('%.6g', value(1));
        elseif isstring(value)
            value = char(value);
        end
        pairs{end+1} = [key '=' value]; %#ok<AGROW>
    end
end

text = normaliseToken(origin);
if isempty(text) && (~isempty(pairs) || ~isempty(tuned))
    text = 'unknown';
end
if ~isempty(pairs)
    text = [text ' (' strjoin(pairs, ', ') ')'];
end
if ~isempty(tuned)
    text = [text ' | tuned=' strjoin(tuned(:)', '+')];
    if ~isempty(prior) && isfinite(prior)
        text = [text sprintf(' (prior=%.6g)', prior)];
    end
end
end

function token = normaliseToken(token)
% Lower case, every run of characters outside [0-9a-z] one underscore, no
% leading or trailing underscore.
token = lower(strtrim(char(token)));
token = regexprep(token, '[^0-9a-z]+', '_');
token = regexprep(token, '^_+|_+$', '');
end
