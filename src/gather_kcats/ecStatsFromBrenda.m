function ecStats = ecStatsFromBrenda(KCATcell,SAcell,KCATmedianCell)
% ecStatsFromBrenda  Summarise loadBRENDAdata's output per EC number, for
% reviewKcatAssignment's ecStats argument.
%
% loadBRENDAdata returns one aggregate of kcat.tsv per call, so the maximum
% view and the median view come from two calls:
%
%     [KCATcell, SAcell] = loadBRENDAdata('kcatAggregation', 'max');
%     KCATmedianCell     = loadBRENDAdata('kcatAggregation', 'median');
%     ecStats = ecStatsFromBrenda(KCATcell, SAcell, KCATmedianCell);
%
% geckopy's ec_stats_from_brenda gives the same statistics. The EC numbers
% keep the case they have in the BRENDA files.
%
% Parameters
% ----------
% KCATcell : cell
%     first output of loadBRENDAdata with kcatAggregation 'max': {ecCode,
%     substrate, organism, kcatMax}.
% SAcell : cell
%     second output of loadBRENDAdata with kcatAggregation 'max': {ecCode,
%     organism, kcatFromSA, mw}.
% KCATmedianCell : cell
%     first output of loadBRENDAdata with kcatAggregation 'median'
%     (optional). Without it kcatMedian is NaN, so reviewKcatAssignment
%     cannot flag 'ec-maximum', and a warning is given.
%
% Returns
% -------
% ecStats : containers.Map
%     EC code -> struct(nKcat,kcatMax,kcatMedian,values,nSa). nKcat and
%     kcatMax are the number and the maximum of the kcat rows of the
%     maximum view, kcatMedian is the median over the rows of the median
%     view, values holds the kcat values of both views rounded to six
%     decimals, and nSa is the number of specific-activity rows.
%
% See also
% --------
% reviewKcatAssignment, loadBRENDAdata
ecCol  = KCATcell{1};
maxCol = KCATcell{4};

hasMedian = nargin >= 3 && ~isempty(KCATmedianCell) && ~isempty(KCATmedianCell{1});
if nargin < 3
    warning('ecStatsFromBrenda:noMedianView', ...
        ['No median view of the kcat rows was given, so kcatMedian is NaN ' ...
         'and reviewKcatAssignment cannot flag ec-maximum.'])
end
medEcCol = {};
medCol   = [];
if hasMedian
    medEcCol = KCATmedianCell{1};
    medCol   = KCATmedianCell{4};
end

saEcCol = {};
if ~isempty(SAcell) && numel(SAcell) >= 1 && ~isempty(SAcell{1})
    saEcCol = SAcell{1};
end

ecList = unique([ecCol(:); saEcCol(:)]);
ecStats = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(ecList)
    ec = char(ecList{i});
    rows = strcmp(ecCol,ec);
    nKcat = sum(rows);
    if nKcat > 0
        kcatMax = max(maxCol(rows));
        values  = round(maxCol(rows),6);
    else
        kcatMax = NaN;
        values  = [];
    end
    kcatMedian = NaN;
    if hasMedian
        medRows = strcmp(medEcCol,ec);
        if any(medRows)
            kcatMedian = median(medCol(medRows),'omitnan');
            values = [values(:); round(medCol(medRows),6)];
        end
    end
    values = unique(values(:));
    nSa = sum(strcmp(saEcCol,ec));
    ecStats(ec) = struct('nKcat',nKcat,'kcatMax',kcatMax,'kcatMedian',kcatMedian,'values',values,'nSa',nSa);
end
end
