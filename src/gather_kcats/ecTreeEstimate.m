function [log10kcat, n, depth, ecUsed] = ecTreeEstimate(eccodes, measurements, k)
% ecTreeEstimate  EC-tree kcat estimates for EC numbers.
%
% For each EC number, the median log10 kcat of its measurements is pulled
% toward the estimate of its parent class with weight n / (n + k), where n
% is the number of measurements: the global median, then the 2-digit
% class, the 3-digit class and the full EC number, in that order. A
% well-measured EC keeps its own median; a sparsely measured one leans on
% its class. geckopy has the same (geckopy.gather_kcats.ec_tree.ECTree).
%
% Parameters
% ----------
% eccodes : cell
%     per query one or more EC numbers separated by ';'. The one with most
%     measurements of its full EC is used (the first on a tie).
% measurements : struct
%     BRENDA measurement table with fields ec and log10kcat (see
%     loadBrendaMeasurements).
% k : double
%     shrinkage constant (default 5): the number of measurements at which
%     an EC's own median and its class count equally.
%
% Returns
% -------
% log10kcat : double
%     the estimates; the global median for queries without measurements.
% n : double
%     number of measurements of the full EC number used.
% depth : double
%     deepest EC level with measurements: 4, 3, 2, or 0 when there are
%     none (or no EC number).
% ecUsed : cell
%     the EC number used.
%
% Examples
% --------
%     m = loadBrendaMeasurements();
%     [v, n, depth] = ecTreeEstimate({'1.1.1.1'; '2.7.1.1;2.7.1.2'}, m);
%     kcat = 10.^v;
%
% See also
% --------
% kcatsFromECTree, loadBrendaMeasurements

if nargin < 3 || isempty(k)
    k = 5;
end
if ~(isnumeric(k) && isscalar(k) && k > 0)
    error('ecTreeEstimate:k', 'k must be a positive number.')
end
if ischar(eccodes) || isstring(eccodes)
    eccodes = cellstr(eccodes);
end
ec = measurements.ec(:);
v = measurements.log10kcat(:);
globalMedian = median(v);
lv = cell(4, 1);
for d = [2 3 4]
    if d == 4
        key = ec;
    else
        key = cellfun(@(e) parentEC(e, d), ec, 'UniformOutput', false);
    end
    [keys, ~, g] = unique(key);
    lv{d} = struct('keys', {keys}, 'median', accumarray(g, v, [], @median), ...
        'n', accumarray(g, 1));
end

nq = numel(eccodes);
log10kcat = repmat(globalMedian, nq, 1);
n = zeros(nq, 1);
depth = zeros(nq, 1);
ecUsed = repmat({''}, nq, 1);
for q = 1:nq
    ecs = strtrim(strsplit(char(eccodes{q}), ';'));
    ecs = ecs(~cellfun(@isempty, ecs));
    if isempty(ecs)
        continue
    end
    best = ecs{1};
    bestN = -1;
    for e = 1:numel(ecs)
        c = countOf(lv{4}, ecs{e});
        if c > bestN
            best = ecs{e};
            bestN = c;
        end
    end
    value = globalMedian;
    for d = [2 3 4]
        if d == 4
            key = best;
        else
            key = parentEC(best, d);
        end
        [c, med] = countOf(lv{d}, key);
        if c > 0
            w = c / (c + k);
            value = w * med + (1 - w) * value;
            depth(q) = d;
        end
    end
    log10kcat(q) = value;
    n(q) = bestN;
    ecUsed{q} = best;
end
end

function p = parentEC(ec, d)
parts = strsplit(ec, '.');
p = strjoin(parts(1:min(d, numel(parts))), '.');
end

function [c, med] = countOf(level, key)
[found, i] = ismember(key, level.keys);
c = 0;
med = NaN;
if found
    c = level.n(i);
    med = level.median(i);
end
end
