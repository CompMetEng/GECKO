function families = getComplexFamilies(model, varargin)
% getComplexFamilies  Group the enzymes of multi-subunit reactions into complexes.
%
% A complex is a set of enzymes that catalyses a reaction together, taken
% from model.ec.rxnEnzMat. Two enzyme sets are one complex when they share
% at least minShared subunits and at least a fraction minOverlap of the
% larger set (isozyme variants that swap a subunit, the same complex on
% different reactions). A single shared protein never merges two sets, so
% electron carriers and shared subunits (cytochrome c, thioredoxin, the E3
% subunit of the 2-oxo-acid dehydrogenases) do not chain unrelated
% complexes together.
%
% A subunit in exactly two complexes belongs to the one whose enzyme sets
% catalyse more reactions (the core subunits COX1-3 of complex IV also
% appear in one assembly reaction). A subunit in three or more complexes,
% in two with equal weight, or that also catalyses a reaction on its own is
% shared: it is left as measured and does not take part in the summary of a
% complex.
%
% When isozyme variants of a complex disagree on the number of copies of a
% subunit, the first variant is used and a warning is given.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure), with
%     ec.rxnEnzMat (copies per subunit).
%
% Name-Value Arguments
% --------------------
% minOverlap : double
%     two enzyme sets are one complex when they share at least this
%     fraction of the larger set's subunits (default 0.5)...
% minShared : double
%     ...and at least this many subunits (default 2).
%
% Returns
% -------
% families : struct
%     the complexes of the model, with the fields:
%
%     - complexes : cell array with one struct per complex, with the fields
%       enzymes (cell array with every subunit that appears in any of its
%       enzyme sets), copies (copies of each of these subunits) and owned
%       (logical, true for the subunits that belong to this complex).
%     - owner : numeric vector with one entry per model.ec.enzymes, the
%       index in complexes of the complex that the subunit belongs to, or 0
%       if it belongs to none.
%     - shared : cell array with the subunits of a complex that belong to
%       none: carriers and subunits shared between complexes.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     families = getComplexFamilies(model);
%     families = getComplexFamilies(model, 'minShared', 3);
%
% See also
% --------
% levelComplexConcs, flexibilizeEnzConcs

p = parseGECKOargs(varargin, { ...
    'minOverlap', []; ...
    'minShared',  []});
minOverlap = p.minOverlap;
minShared  = p.minShared;
if isempty(minOverlap); minOverlap = 0.5; end
if isempty(minShared);  minShared  = 2;   end

if ~isfield(model, 'ec') || ~isfield(model.ec, 'rxnEnzMat') || isempty(model.ec.rxnEnzMat)
    error('getComplexFamilies:noRxnEnzMat', ...
        'The model has no ec.rxnEnzMat. Run makeEcModel first.')
end
enzymes = model.ec.enzymes;
matrix  = sparse(model.ec.rxnEnzMat);

% Distinct enzyme sets of the reactions with more than one enzyme, in order
% of first appearance, with the copies of the first reaction that has the set
setIdx     = {};
setCopies  = {};
setCount   = zeros(0, 1);
setLookup  = containers.Map('KeyType', 'char', 'ValueType', 'double');
alone      = false(numel(enzymes), 1);
for row = 1:size(matrix, 1)
    [~, cols, vals] = find(matrix(row, :));
    if numel(cols) == 1
        alone(cols) = true;
    elseif ~isempty(cols)
        key = sprintf('%d,', cols);
        if isKey(setLookup, key)
            k = setLookup(key);
        else
            k = numel(setIdx) + 1;
            setLookup(key) = k;
            setIdx{k, 1}    = cols(:);
            setCopies{k, 1} = full(vals(:));
            setCount(k, 1)  = 0;
        end
        setCount(k) = setCount(k) + 1;
    end
end

% Merge the enzyme sets that overlap enough
nSets  = numel(setIdx);
parent = 1:nSets;
byEnzyme = cell(numel(enzymes), 1);
for i = 1:nSets
    for e = setIdx{i}'
        byEnzyme{e}(end+1) = i;
    end
end
for e = 1:numel(byEnzyme)
    indices = byEnzyme{e};
    for x = 1:numel(indices)
        for y = x+1:numel(indices)
            i = indices(x);
            j = indices(y);
            common = numel(intersect(setIdx{i}, setIdx{j}));
            if common >= minShared && common >= minOverlap * max(numel(setIdx{i}), numel(setIdx{j}))
                parent(findRoot(parent, i)) = findRoot(parent, j);
            end
        end
    end
end
roots = zeros(1, nSets);
for i = 1:nSets
    roots(i) = findRoot(parent, i);
end

% One complex per group of sets, numbered by the first set of the group
[groupRoots, ~, groupOf] = unique(roots, 'stable');
nComplex  = numel(groupRoots);
complexes = cell(nComplex, 1);
weight    = cell(nComplex, 1);
for g = 1:nComplex
    idx = zeros(0, 1);
    copies = zeros(0, 1);
    w = zeros(numel(enzymes), 1);
    for i = find(groupOf(:)' == g)
        w(setIdx{i}) = w(setIdx{i}) + setCount(i);
        for n = 1:numel(setIdx{i})
            e = setIdx{i}(n);
            pos = find(idx == e, 1);
            if isempty(pos)
                idx(end+1, 1)    = e; %#ok<AGROW>
                copies(end+1, 1) = setCopies{i}(n); %#ok<AGROW>
            elseif copies(pos) ~= setCopies{i}(n)
                warning('getComplexFamilies:copiesDisagree', ...
                    '%s has %g copies in one variant and %g in another; the first is used.', ...
                    enzymes{e}, copies(pos), setCopies{i}(n))
            end
        end
    end
    complexes{g} = struct('enzymes', {enzymes(idx)}, 'copies', copies, 'idx', idx);
    weight{g} = w;
end

% A subunit belongs to the complex with most reactions when it is in two
owner = zeros(numel(enzymes), 1);
inComplex = false(numel(enzymes), 1);
for e = 1:numel(enzymes)
    hits = [];
    for g = 1:nComplex
        if weight{g}(e) > 0
            hits(end+1, :) = [weight{g}(e), g]; %#ok<AGROW>
        end
    end
    if isempty(hits)
        continue
    end
    inComplex(e) = true;
    if alone(e) || size(hits, 1) > 2
        continue
    end
    if size(hits, 1) == 1
        owner(e) = hits(1, 2);
    elseif hits(1, 1) ~= hits(2, 1)
        [~, best] = max(hits(:, 1));
        owner(e) = hits(best, 2);
    end
end
for g = 1:nComplex
    complexes{g}.owned = owner(complexes{g}.idx) == g;
    complexes{g} = rmfield(complexes{g}, 'idx');
end

families.complexes = complexes;
families.owner     = owner;
families.shared    = enzymes(inComplex & owner == 0);
end

function r = findRoot(parent, i)
% Root of the union-find tree that contains i
r = i;
while parent(r) ~= r
    r = parent(r);
end
end
