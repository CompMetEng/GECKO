%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KCATcell, SAcell] = loadBRENDAdata(varargin)
% loadBRENDAdata  Load kcat and specific activity data from BRENDA files.
%
% Reads the BRENDA data files (kcat.tsv, sa.tsv and mw.tsv) from the BRENDA
% database folder defined by the model adapter, and returns the kcat values
% and specific activities (the latter converted to kcat values using the
% molecular weights).
%
% kcat.tsv and sa.tsv each carry both a max and a median aggregate per
% (EC, substrate, organism) triple; kcatAggregation selects which of the
% two is read. Each file starts with a `#`-prefixed release-version line
% followed by a tab-delimited column header, both skipped on read.
%
% Name-Value Arguments
% --------------------
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% kcatAggregation : char
%     which aggregate to read from kcat.tsv and sa.tsv, 'max' or
%     'median' (default: params.kcatAggregation from the model adapter,
%     or 'max' if the adapter does not set it). A triple's max is the
%     highest turnover reported for it; its median is the middle one.
% filterOutlierKcats : logical
%     whether to drop a row whose kcat is implausibly low for its EC
%     (default: true). A row is dropped when its kcat is more than four
%     orders of magnitude below the median kcat of its EC. The reference
%     median of an EC is taken over the positive kcat values of the max and
%     median views of kcat.tsv and of the kcat values derived from sa.tsv
%     and mw.tsv, before any row is dropped, so it does not depend on
%     kcatAggregation. Wildcard ECs (containing '-') do not contribute to a
%     reference. The EC level widens from w.x.y.z to w.x.y, w.x and w while
%     fewer than five values are available. Set to false to keep every row.
%
% Returns
% -------
% KCATcell : cell
%     kcat data extracted from BRENDA.
% SAcell : cell
%     specific activity data extracted from BRENDA.
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     [KCATcell, SAcell] = loadBRENDAdata();
%     [KCATcell, SAcell] = loadBRENDAdata('modelAdapter', adapter);
%     [KCATcell, SAcell] = loadBRENDAdata('kcatAggregation', 'median');
%     [KCATcell, SAcell] = loadBRENDAdata('filterOutlierKcats', false);
%
% See also
% --------
% loadDatabases, getECfromDatabase

p = parseGECKOargs(varargin, { ...
    'modelAdapter',    [], @(x) assert(isempty(x) || isa(x,'ModelAdapter'), ...
                                'loadBRENDAdata:invalidModelAdapter', ...
                                'modelAdapter must be a ModelAdapter.'); ...
    'kcatAggregation', [], @(x) assert(isempty(x) || any(strcmpi(x, ...
                                {'max','median'})), ...
                                'loadBRENDAdata:invalidKcatAggregation', ...
                                'kcatAggregation must be ''max'' or ''median''.'); ...
    'filterOutlierKcats', true, @(x) assert((islogical(x) || isnumeric(x)) && ...
                                isscalar(x), ...
                                'loadBRENDAdata:invalidFilterOutlierKcats', ...
                                'filterOutlierKcats must be a logical scalar.')});
modelAdapter       = p.modelAdapter;
kcatAggregation    = p.kcatAggregation;
filterOutlierKcats = logical(p.filterOutlierKcats);

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end

kcatAggregation = resolveKcatAggregation(kcatAggregation, modelAdapter);
%Column 4 holds the per-triple maximum and column 5 the median; both are
%read, so that the outlier filter can use the max and the median views.
valueCol = 4 + strcmpi(kcatAggregation,'median');

basePath      = modelAdapter.getBrendaDBFolder();
KCAT_file      = fullfile(basePath,'kcat.tsv');
SA_file        = fullfile(basePath,'sa.tsv');
MW_file        = fullfile(basePath,'mw.tsv');

%Extract BRENDA DATA from files information. kcat.tsv/sa.tsv are seven
%columns (ec_code, substrate, organism, value_max, value_median, n,
%references); the data cells hold ec_code, substrate, organism, the max and
%the median. mw.tsv is six columns (no substrate-level aggregation choice):
%ec_code, substrate, organism, value, n, references.
KCAT          = openDataFile(KCAT_file,1,'%q %q %q %f %f %f %q',[4 5]);
scalingFactor = 1/60;    %[umol/min/mg] -> [mmol/s/g]
SA            = openDataFile(SA_file,scalingFactor,'%q %q %q %f %f %f %q',[4 5]);
scalingFactor = 1/1000;  %[g/mol] -> [g/mmol]
MW            = openDataFile(MW_file,scalingFactor,'%q %q %q %f %f %q',4);

%build an index on MW{1} to speed things up a bit
%first just extract the genus (i.e. the first part of the name)
MWECNum = unique(upper(MW{1}));
MWECNumIndices = cell(length(MWECNum),1);
MWECNumHashMap = containers.Map(MWECNum,1:length(MWECNum));
for i = 1:length(MW{1})
    matchInd = cell2mat(values(MWECNumHashMap, upper(MW{1}(i))));
    MWECNumIndices{matchInd} = [MWECNumIndices{matchInd};i];
end

%Match each specific activity to the molecular weight of its EC and
%organism. The first match is used.
previousEC = []; EC_indexes = [];
saRows = zeros(length(SA{1}),1);
mwVals = zeros(length(SA{1}),1);
nMatched = 0;
for i=1:length(SA{1})
    %Gets the indexes of the EC repetitions in the MW cell for every
    %new (different) EC
    if ~strcmpi(SA{1}(i), previousEC)
        key = upper(SA{1}(i));
        if isKey(MWECNumHashMap,key) %annoyingly, this seems to be needed
            matchInd = cell2mat(values(MWECNumHashMap,key));
            EC_indexes = MWECNumIndices{matchInd};
        else
            EC_indexes = [];
        end
    end
    mwEC{1} = MW{3}(EC_indexes); mwEC{2} = MW{4}(EC_indexes);
    % just looks for the first match because just the maximal value for
    % each EC# / Orgaism is reported on the file
    org_index = find(strcmpi(SA{3}(i),mwEC{1}),1);
    if ~isempty(org_index)
        nMatched = nMatched + 1;
        saRows(nMatched) = i;
        mwVals(nMatched) = mwEC{2}(org_index);
    end
    previousEC = SA{1}(i);
end
saRows = saRows(1:nMatched);
mwVals = mwVals(1:nMatched);

%Specific activities as kcat values [1/s], for the max and the median
saEC     = SA{1}(saRows);
saMax    = SA{4}(saRows).*mwVals;
saMedian = SA{5}(saRows).*mwVals;

if valueCol == 4
    saSelected = saMax;
else
    saSelected = saMedian;
end

%Drop kcats that are implausibly low for their EC
keepKcat = true(size(KCAT{1}));
keepSA   = true(size(saEC));
if filterOutlierKcats
    ecPool  = [KCAT{1}; KCAT{1}; saEC; saEC];
    valPool = [KCAT{4}; KCAT{5}; saMax; saMedian];
    keepKcat = ~isOutlier(KCAT{1},KCAT{valueCol},ecPool,valPool);
    keepSA   = ~isOutlier(saEC,saSelected,ecPool,valPool);
end

KCATcell = {KCAT{1}(keepKcat), KCAT{2}(keepKcat), KCAT{3}(keepKcat), ...
            KCAT{valueCol}(keepKcat)};
if nMatched > 0
    SAcell = {saEC(keepSA), SA{3}(saRows(keepSA)), saSelected(keepSA), ...
              mwVals(keepSA)};
else
    SAcell = {[],[],[],[]};
end
end

function data_cell = openDataFile(fileName,scalingFactor,formatSpec,valueCols)
fID          = fopen(fileName);
raw          = textscan(fID,formatSpec,'delimiter','\t','HeaderLines',2);
fclose(fID);
data_cell    = [raw(1:3), cellfun(@(c) raw{c}*scalingFactor, num2cell(valueCols), ...
                                  'UniformOutput', false)];
end

function outlier = isOutlier(ecs,vals,ecPool,valPool)
%A value is an outlier if it is more than 'decades' orders of magnitude below
%the median kcat of its own EC, taken over the pool of values. The EC widens
%from the full code to its 3-, 2- and 1-level prefix while fewer than minRef
%values are available.
decades = 4;
minRef  = 5;
usePool = valPool > 0 & ~contains(ecPool,'-');
ecPool  = ecPool(usePool);
logPool = log10(valPool(usePool));
refLog  = NaN(size(ecs));
if isempty(ecs) || isempty(logPool)
    outlier = false(size(ecs));
    return
end
for level = 4:-1:1
    poolKey = ecPrefix(ecPool,level);
    [grp,keys] = findgroups(poolKey);
    nPerKey    = accumarray(grp,1);
    medPerKey  = splitapply(@median,logPool,grp);
    [found,loc] = ismember(ecPrefix(ecs,level),keys);
    use = isnan(refLog) & found;
    use(use) = nPerKey(loc(use)) >= minRef;
    refLog(use) = medPerKey(loc(use));
end
outlier = vals > 0 & ~isnan(refLog) & log10(vals) < refLog - decades;
end

function prefix = ecPrefix(ecs,level)
%The first 'level' dot-separated parts of each EC code
prefix = regexp(ecs,['^[^.]*(?:\.[^.]*){0,' num2str(level-1) '}'],'match','once');
prefix(cellfun(@isempty,prefix)) = {''};
prefix = reshape(prefix,size(ecs));
end
