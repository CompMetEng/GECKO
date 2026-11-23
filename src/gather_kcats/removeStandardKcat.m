function model = removeStandardKcat(model)
% removeStandardKcat  Remove the standard pseudoenzyme and standard kcats.
%
% Low-level: it undoes getStandardKcat, which collectKcats and
% assignKcats run; call it directly only to change kcats by hand.
%
% Removes the "standard" pseudoenzyme and standard kcat and MW values, as
% they were introduced by getStandardKcat. Also standard kcat values that
% were assigned by getStandardKcat if fillZeroKcat was set to true are
% removed from model.ec.kcat. Both the model.ec and model.S structures are
% modified. A standard kcat value is recognised by the origin 'standard' of
% model.ec.source, compared without regard to letter case; a row with a
% tuning status is kept.
%
% Parameters
% ----------
% model : struct
%     an ecModel (with ecModel.ec structure) that has
%     standard kcat values implemented by getStandardKcat.
%
% Returns
% -------
% model : struct
%     ecModel without standard pseudoprotein and standard kcat values.
%
% Examples
% --------
%     model = removeStandardKcat(model);
%
% See also
% --------
% getStandardKcat, applyKcatConstraints

% Remove standard enzyme from ec structure
stdEnzIdx = find(strcmpi(model.ec.enzymes, 'standard'));
if ~isempty(stdEnzIdx)
    model.ec.genes(stdEnzIdx)       = [];
    model.ec.enzymes(stdEnzIdx)     = [];
    model.ec.mw(stdEnzIdx)          = [];
    model.ec.sequence(stdEnzIdx)    = [];
    if isfield(model.ec,'concs')
        model.ec.concs(stdEnzIdx)   = [];
    end
    rxnEnzIdx = find(model.ec.rxnEnzMat(:,stdEnzIdx));
    model.ec.rxns(rxnEnzIdx)        = [];
    model.ec.kcat(rxnEnzIdx)        = [];
    model.ec.source(rxnEnzIdx)      = [];
    model.ec.notes(rxnEnzIdx)       = [];
    model.ec.eccodes(rxnEnzIdx)     = [];
    model.ec.rxnEnzMat(:,stdEnzIdx) = [];
    model.ec.rxnEnzMat(rxnEnzIdx,:) = [];
end

% Remove standard kcat values and reapply kcat constraints for those
% specific reactions
% A row with a tuning status was changed on purpose and is kept.
stdSource  = parseKcatSource(model.ec.source);
stdKcatIdx = find(strcmpi({stdSource.origin}, 'standard') & cellfun(@isempty, {stdSource.tuned}));
stdKcatIdx = stdKcatIdx(:);
if ~isempty(stdKcatIdx)
    model.ec.source(stdKcatIdx)     = {''};
    model.ec.kcat(stdKcatIdx)       = 0;
    model = applyKcatConstraints(model,stdKcatIdx);
end

% Remove standard protein, usage and gene from model structure
stdMetIdx = find(strcmpi(model.mets, 'prot_standard'));
if ~isempty(stdMetIdx)
    model       = removeMets(model,stdMetIdx,false,false,false,false);
end
stdProtEx = find(strcmpi(model.rxns, 'usage_prot_standard'));
if ~isempty(stdProtEx)
    model       = removeReactions(model,stdProtEx,false,false,false);
end
stdProtGene = find(strcmpi(model.genes, 'standard'));
if ~isempty(stdProtGene)
    model       = removeGenes(model,stdProtGene,false,false,false);
end
end
