function protData = loadProtData(replPerCond, varargin)
% loadProtData  Load absolute proteomics data and average over replicates.
%
% Loads absolute proteomics data (in mg/gDCW) and returns mean values
% across replicates for each condition in the data file. By default it also
% filters the data by various criteria, to remove uncertain data (see input
% parameters).
%
% Parameters
% ----------
% replPerCond : double
%     vector with number of replicates for each condition in the dataset.
%     Example: [3, 2] if first condition has triplicates and second
%     condition has duplicates. The sum must equal the number of replicate
%     columns of the file (all columns after the UniProt ID), otherwise an
%     error is raised.
%
% Name-Value Arguments
% --------------------
% protDataFile : char
%     path to file with proteomics data, where protein levels are in
%     mg/gDCW (default reads data/proteomics.tsv as specified in
%     modelAdapter). Alternatively, protDataFile can be a protData structure
%     that was previously made by loadProtData.
% filterData : logical
%     whether abundances should be filtered. If false, minVal, maxRSD,
%     maxMissing, cutLowest and addStdevs are not considered, and
%     zerosAsMeasured applies in both cases. Default true.
% modelAdapter : ModelAdapter
%     a loaded model adapter (default: the current default model adapter).
% minVal : double
%     minimum mean protein measurement per condition; proteins with a lower
%     mean abundance are discarded (set to NaN) for that condition (default
%     0, i.e. no lower threshold).
% maxRSD : double
%     maximum relative standard deviation (RSD, the standard deviation
%     across replicates divided by their mean) allowed per condition;
%     proteins exceeding it are discarded (set to NaN) for that condition
%     (default 1, i.e. RSD up to 100%).
% maxMissing : double
%     minimum fraction of replicates within a condition that must have a
%     measured value; proteins with fewer present replicates than this
%     fraction requires are discarded (set to NaN) for that condition
%     (default 2/3, i.e. at least two of three replicates must have a
%     value). If conditions have different number of replicates (as
%     indicated in replPerCond), maxMissing can also be a vector of the
%     same length as replPerCond, with individual values per condition.
% cutLowest : double
%     percentage of lowest mean values per condition to be discarded (not
%     considering NaN values) (default 5).
% addStdevs : double
%     how many standard deviations should be added to the mean value of
%     each protein measurement across replicates, broadening the confidence
%     interval (default 1).
% zerosAsMeasured : logical
%     whether an abundance of 0 is a measurement (default false). If false,
%     a 0 is read as not measured (NaN) before any filter, also in a protDataFile
%     that is a protData structure: the mean of a protein over its replicates
%     leaves out the zeros, and a protein with only zeros is removed. If true,
%     a 0 counts as a measured replicate and is kept in the mean; a protein
%     whose replicates are all 0 is returned as 0, is not removed by minVal or
%     cutLowest, and limits the usage of its enzyme to 0 in constrainEnzConcs.
%     Proteomics reports a protein below the detection limit as 0 or as
%     missing. With true, a 0 is read as a protein that is absent, for
%     example in a deletion strain. maxRSD and maxMissing still apply to a
%     protein with some zero replicates.
%
% Returns
% -------
% protData : struct
%     structure with (filtered) proteome data.
%
% Notes
% -----
% The protData structure has the following fields:
%
% - uniprotIDs : cell array with Uniprot IDs matching protData.abundances.
% - abundances : matrix of proteomics data, where each column contains mean abundances per condition.
%
% A protein that is not measured in a condition is NaN in abundances and
% leaves the usage of its enzyme drawing from the protein pool (see
% fillEnzConcs and constrainEnzConcs).
%
% Examples
% --------
%     % optional arguments may be given positionally or as name-value pairs:
%     protData = loadProtData(replPerCond, protDataFile, filterData, modelAdapter, minVal, maxRSD, maxMissing, cutLowest, addStdevs, zerosAsMeasured);
%     protData = loadProtData(replPerCond, 'filterData', false);
%     protData = loadProtData(replPerCond, 'zerosAsMeasured', true);
%
% See also
% --------
% loadFluxData

p = parseGECKOargs(varargin, { ...
    'protDataFile', []; ...
    'filterData',   true; ...
    'modelAdapter', []; ...
    'minVal',       0; ...
    'maxRSD',       1; ...
    'maxMissing',   2/3; ...
    'cutLowest',    5; ...
    'addStdevs',    1; ...
    'zerosAsMeasured', false});
protDataFile = p.protDataFile;
filterData   = p.filterData;
modelAdapter = p.modelAdapter;
minVal       = p.minVal;
maxRSD       = p.maxRSD;
maxMissing   = p.maxMissing;
cutLowest    = p.cutLowest;
addStdevs    = p.addStdevs;
zerosAsMeasured = p.zerosAsMeasured;
if isempty(zerosAsMeasured)
    zerosAsMeasured = false;
end

if isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if isempty(protDataFile)
    protDataFile = fullfile(params.path,'data','proteomics.tsv');
end

format = '%s';
for i=1:sum(replPerCond)
    format = [format ' %f'];
end
if ~isstruct(protDataFile)
    checkFileExistence(protDataFile);
    fID         = fopen(protDataFile);
    nHeaderCols = numel(strsplit(fgetl(fID),'\t','CollapseDelimiters',false));
    if nHeaderCols ~= sum(replPerCond)+1
        fclose(fID);
        error('loadProtData:replicateMismatch', ...
              ['The header of %s has %d columns but replPerCond sums to %d, ' ...
               'so %d columns are expected (the UniProt ID and one per replicate).'], ...
              protDataFile, nHeaderCols, sum(replPerCond), sum(replPerCond)+1)
    end
    frewind(fID);
    protDataRaw = textscan(fID,format,'Delimiter','\t','HeaderLines',1,'TreatAsEmpty',{'NA','na','NaN','#VALUE!'});
    uniprotIDs  = protDataRaw{1};
    abundances  = cell2mat(protDataRaw(2:end));
    fclose(fID);
else
    uniprotIDs  = protDataFile.uniprotIDs;
    abundances  = protDataFile.abundances;
end

% An abundance of 0 is not measured unless zerosAsMeasured is true
if ~zerosAsMeasured
    abundances(abundances == 0) = NaN;
end

%Remove entriew without ID
remData = cellfun(@isempty,uniprotIDs);
uniprotIDs(remData,:) = [];
abundances(remData,:) = [];
m                     = size(abundances,1);
filtAbund             = nan(m,numel(replPerCond));

if filterData
    for i=1:numel(replPerCond)
        condAbund    = abundances(:,1:replPerCond(i));
        if i<numel(replPerCond)
            abundances   = abundances(:,replPerCond(i)+1:end);
        end
        % First filter maxMissing
        if size(condAbund,2) > 1
            if numel(maxMissing)>1
                maxMisRepl = maxMissing(i);
            else
                maxMisRepl = maxMissing;
            end
            if zerosAsMeasured
                present = condAbund>=0;
            else
                present = condAbund>0;
            end
            remData = sum(present,2)<maxMisRepl*size(condAbund,2);
            condAbund(remData,:) = nan;
        end
        % Filter by RSD
        remData = (std(condAbund,0,2,'omitnan')./mean(condAbund,2,'omitnan'))>maxRSD;
        condAbund(remData,:) = nan;
        % Add stdevs
        condAbund = mean(condAbund,2,'omitnan')+(addStdevs*std(condAbund,0,2,'omitnan'));
        % Filter by minVal
        remData = mean(condAbund,2,'omitnan')<minVal;
        if zerosAsMeasured
            remData = remData & condAbund~=0;
        end
        condAbund(remData) = nan;
        % Remove the lowest cutLowest percent
        if zerosAsMeasured
            cutData = find(~isnan(condAbund) & condAbund~=0);
        else
            cutData = find(~isnan(condAbund));
        end
        [~,sortData] = sort(condAbund(cutData));
        lowCutoff = floor(numel(cutData)*cutLowest/100);
        condAbund(cutData(sortData(1:lowCutoff))) = nan;
        % Combine conditions
        filtAbund(:,i) = condAbund;
    end
else
    for i=1:numel(replPerCond)
        condAbund    = abundances(:,1:replPerCond(i));
        if i<numel(replPerCond)
            abundances = abundances(:,replPerCond(i)+1:end);
        end
        filtAbund(:,i) = mean(condAbund,2,'omitnan');
    end
end
notAllNan = logical(sum(~isnan(filtAbund),2));
protData.abundances = filtAbund(notAllNan,:);
protData.uniprotIDs = uniprotIDs(notAllNan);
end

