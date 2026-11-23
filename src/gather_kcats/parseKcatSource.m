function info = parseKcatSource(source)
% parseKcatSource  Read the provenance text of a kcat in model.ec.source.
%
% Reads the form written by formatKcatSource,
%
%     <origin>[ (<detail>)][ | tuned=<method>[+<method>...] (prior=<value>)]
%
% A text that does not follow the form is read as an origin without detail
% and without tuning status, so a text that GECKO wrote in the form of an
% earlier release stays readable. A detail pair is ignored when its key does
% not start with a letter, a value is a number when it is in plain decimal
% form (for example 0.61, -2, .5 or 1e3) and text otherwise, and the prior
% needs exactly one space before '(prior='.
%
% Parameters
% ----------
% source : char, string or cell
%     one provenance text, or a cell array of them, e.g. model.ec.source.
%     Missing, NaN and other non-text input gives an empty origin.
%
% Returns
% -------
% info : struct
%     scalar struct for a char or string, and a struct array of the size of
%     the cell array for a cell array, with the fields
%
%     - origin : the origin as written, in the case it was written in.
%     - detail : scalar struct of the detail pairs; a value is a double when
%       it parses as a number and char otherwise.
%     - tuned : 1xN cell of the tuning methods, {} when none.
%     - prior : the kcat before the first tuning, NaN when none.
%
% Examples
% --------
%     info = parseKcatSource('brenda (wc=0, origin=2) | tuned=cmaes (prior=6.55)');
%     info.origin       % 'brenda'
%     info.detail.wc    % 0
%     info.tuned        % {'cmaes'}
%     info.prior        % 6.55
%
% See also
% --------
% formatKcatSource, getStandardKcat, resolveTunableContext
if iscell(source)
    info = repmat(emptyInfo(), size(source));
    for i = 1:numel(source)
        info(i) = parseText(source{i});
    end
else
    info = parseText(source);
end
end

function info = emptyInfo()
info = struct('origin', '', 'detail', struct(), 'tuned', {{}}, 'prior', NaN);
end

function info = parseText(text)
info = emptyInfo();
if isstring(text) && isscalar(text) && ~ismissing(text)
    text = char(text);
end
if ~ischar(text) || (~isempty(text) && ~isrow(text))
    return
end
text = strtrim(text);

tail = '';
cut = strfind(text, ' | tuned=');
if ~isempty(cut)
    tail = text(cut(1)+9:end);
    text = strtrim(text(1:cut(1)-1));
end

open = strfind(text, ' (');
if ~isempty(open) && endsWith(text, ')')
    info.origin = strtrim(text(1:open(1)-1));
    parts = strsplit(text(open(1)+2:end-1), ',');
    for i = 1:numel(parts)
        eq = find(parts{i} == '=', 1);
        if isempty(eq)
            continue
        end
        key = strtrim(parts{i}(1:eq-1));
        if isempty(regexp(key, '^[A-Za-z][A-Za-z0-9_]*$', 'once'))
            continue
        end
        info.detail.(key) = parseValue(strtrim(parts{i}(eq+1:end)));
    end
else
    info.origin = strtrim(text);
end

tail = strtrim(tail);
if ~isempty(tail)
    at = strfind(tail, ' (prior=');
    if ~isempty(at)
        names = tail(1:at(1)-1);
        rest = tail(at(1)+8:end);
        if endsWith(rest, ')') && ~any(rest(1:end-1) == ')') && ~(~isempty(names) && isspace(names(end)))
            tail = names;
            prior = parseValue(strtrim(rest(1:end-1)));
            if isnumeric(prior)
                info.prior = prior;
            end
        end
    end
    methods = strtrim(strsplit(tail, '+'));
    methods(cellfun(@isempty, methods)) = [];
    if ~isempty(methods)
        info.tuned = methods;
    end
end
end

function value = parseValue(text)
% A number when the text is a plain decimal number, the text otherwise.
if ~isempty(regexp(text, '^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$', 'once'))
    value = str2double(text);
else
    value = text;
end
end
