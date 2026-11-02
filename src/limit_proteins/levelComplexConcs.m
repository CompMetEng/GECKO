function [model, result] = levelComplexConcs(model, varargin)
% levelComplexConcs  Give every enzyme complex one abundance.
%
% Absolute proteomics recovers the subunits of a complex unevenly, so after
% fillEnzConcs the subunits of one complex can disagree several fold about
% how much of it there is. Capped one by one, the complex can then carry no
% more flux than its least abundant subunit allows. This function groups
% the enzymes of every multi-subunit reaction into complexes (see
% getComplexFamilies) and gives each complex one abundance, scaled to every
% subunit by its stoichiometry.
%
% Run after fillEnzConcs and before constrainEnzConcs. The concentrations
% (mg/gDCW) in model.ec.concs are converted to molar amounts with model.ec.mw
% for the levelling and back afterwards. Per complex, the statistic over the
% measured subunits that belong to the complex of the molar amount divided
% by the copies of the subunit is written back, times the copies and the
% molecular weight, to every subunit that belongs to the complex. Subunits
% that are shared between complexes, or that also catalyse a reaction on
% their own, are left as measured.
%
% Parameters
% ----------
% model : struct
%     an ecModel in GECKO 3 format (with ecModel.ec structure), with
%     ec.concs filled.
%
% Name-Value Arguments
% --------------------
% statistic : char
%     the summary over the measured subunits, 'median' or 'mean' (default
%     'median').
% fillMissing : logical
%     whether subunits without a measurement, that belong to a complex, are
%     given the concentration of the complex (default true).
% families : struct
%     the complexes of the model, as returned by getComplexFamilies
%     (default: calculated from the model).
%
% Returns
% -------
% model : struct
%     ecModel where model.ec.concs has one abundance per complex.
% result : struct
%     what was changed, with the fields:
%
%     - adjusted : cell array with the measured subunits whose
%       concentration was replaced by that of the complex.
%     - added : cell array with the unmeasured subunits that were given the
%       concentration of the complex.
%     - shared : cell array with the measured subunits that were left as
%       measured because they belong to more than one complex or also
%       catalyse a reaction on their own.
%     - skipped : cell array with one cell array of subunits per complex
%       that has no measured subunit that belongs to it.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     model = levelComplexConcs(model);
%     model = levelComplexConcs(model, 'statistic', 'mean');
%     [model, result] = levelComplexConcs(model, 'fillMissing', false);
%
% See also
% --------
% getComplexFamilies, fillEnzConcs, constrainEnzConcs, flexibilizeEnzConcs

p = parseGECKOargs(varargin, { ...
    'statistic',   []; ...
    'fillMissing', []; ...
    'families',    []});
statistic   = p.statistic;
fillMissing = p.fillMissing;
families    = p.families;
if isempty(statistic);   statistic   = 'median'; end
if isempty(fillMissing); fillMissing = true;     end
statistic = char(statistic);
if ~ismember(statistic, {'median','mean'})
    error('levelComplexConcs:statistic', ...
        'statistic must be ''mean'' or ''median'', not ''%s''.', statistic)
end
summarise = str2func(statistic);
if isempty(families)
    families = getComplexFamilies(model);
end

if ~isfield(model.ec, 'concs') || numel(model.ec.concs) ~= numel(model.ec.enzymes)
    error('levelComplexConcs:noConcs', ...
        'model.ec.concs has no entry per enzyme. Run fillEnzConcs first.')
end
enzymes = model.ec.enzymes;
mw      =model.ec.mw(:);
concs   = model.ec.concs(:);
molar   = concs ./ mw;

% Measured amounts by enzyme: finite entries, in order of the enzymes
values   = molar;
isKnown  = isfinite(molar);
isLevel  = isKnown & molar > 0;
result.adjusted = cell(0, 1);
result.added    = cell(0, 1);
result.shared   = sort(enzymes(ismember(enzymes, families.shared) & isKnown));
result.skipped  = cell(0, 1);

for c = 1:numel(families.complexes)
    cplx = families.complexes{c};
    [~, pos] = ismember(cplx.enzymes, enzymes);
    owned = find(cplx.owned);
    ownedPos = pos(owned);
    measured = isLevel(ownedPos);
    if ~any(measured)
        result.skipped{end+1, 1} = sort(cplx.enzymes(:))'; %#ok<AGROW>
        continue
    end
    perCopy = summarise(values(ownedPos(measured)) ./ cplx.copies(owned(measured)));
    for n = 1:numel(owned)
        e = ownedPos(n);
        level = perCopy * cplx.copies(owned(n));
        if isLevel(e)
            if abs(values(e) - level) > 1e-8 + 1e-5 * abs(level)
                result.adjusted{end+1, 1} = enzymes{e}; %#ok<AGROW>
            end
            values(e) = level;
        elseif fillMissing
            values(e)  = level;
            isKnown(e) = true;
            result.added{end+1, 1} = enzymes{e}; %#ok<AGROW>
        end
    end
end

model.ec.concs(isKnown) = values(isKnown) .* mw(isKnown);
end
