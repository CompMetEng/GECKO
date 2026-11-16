%run this test case with the command
%results = runtests('adapterTomlTests.m')

function tests = adapterTomlTests
    tests = functiontests(localfunctions);
end

function testReadAdapterToml_tc0061(testCase)
    % Every TOML feature that model adapter files use.
    file = writeTemp({ ...
        '# comment', ...
        'name = "a \"quoted\" # not a comment"  # comment', ...
        'literal = ''C:\path''', ...
        'int = 559_292', ...
        'float = 2.5e-1', ...
        'yes = true', ...
        'list = ["a", "b",', ...
        '        "c"]  # over two lines', ...
        'nums = [0.4, 0.2, 0.1]', ...
        'empty = []', ...
        '', ...
        '[kegg]', ...
        'id = "sce"', ...
        '', ...
        '[evotune.source_groups.okp]', ...
        'sources = ["OKP"]', ...
        'sigma0_log = 0.25', ...
        '', ...
        '[[stage]]', ...
        'action = "select"', ...
        '[[stage]]', ...
        'action = "fill"', ...
        'sources = ["isozymes"]'});
    d = readAdapterToml(file);
    verifyEqual(testCase, d.name, 'a "quoted" # not a comment')
    verifyEqual(testCase, d.literal, 'C:\path')
    verifyEqual(testCase, d.int, 559292)
    verifyEqual(testCase, d.float, 0.25)
    verifyEqual(testCase, d.yes, true)
    verifyEqual(testCase, d.list, {'a','b','c'})
    verifyEqual(testCase, d.nums, [0.4 0.2 0.1])
    verifyEqual(testCase, d.empty, {})
    verifyEqual(testCase, d.kegg.id, 'sce')
    verifyEqual(testCase, d.evotune.source_groups.okp.sources, {'OKP'})
    verifyEqual(testCase, d.evotune.source_groups.okp.sigma0_log, 0.25)
    verifyEqual(testCase, numel(d.stage), 2)
    verifyEqual(testCase, d.stage(2).action, 'fill')
    verifyEqual(testCase, d.stage(2).sources, {'isozymes'})
end

function testReadAdapterTomlErrors_tc0062(testCase)
    % Unsupported TOML stops with the line number.
    verifyError(testCase, @() readAdapterToml(writeTemp({'x = { a = 1 }'})), 'GECKO:readAdapterToml')
    verifyError(testCase, @() readAdapterToml(writeTemp({'x = 2026-10-15'})), 'GECKO:readAdapterToml')
    verifyError(testCase, @() readAdapterToml(writeTemp({'x = unquoted'})), 'GECKO:readAdapterToml')
    verifyError(testCase, @() readAdapterToml(writeTemp({'x = 1', 'x = 2'})), 'GECKO:readAdapterToml')
    try
        readAdapterToml(writeTemp({'a = 1', 'b = """x"""'}));
    catch err
        verifySubstring(testCase, err.message, 'line 2')
    end
end

function testWriteAdapterToml_tc0063(testCase)
    % What writeAdapterToml writes, readAdapterToml reads back unchanged.
    s.text = 'with "quotes" and \ backslash';
    s.number = 0.1;
    s.integer = 150;
    s.flag = false;
    s.list = {'dlkcat','brenda'};
    s.vector = [0.4 0.2];
    s.table.inner.value = 1e-8;
    file = [tempname '.toml'];
    writeAdapterToml(s, file, 'header');
    verifyEqual(testCase, readAdapterToml(file), s)
end

function testLoadAdapterParamsDefaults_tc0064(testCase)
    % A minimal file gets geckopy's defaults, under GECKO's field names.
    folder = tempname; mkdir(folder);
    writeTemp({'conv_gem = "models/x.xml"', 'org_name = "testus testus"', 'p_tot = 0.4'}, ...
        fullfile(folder, 'model_adapter.toml'));
    params = loadAdapterParams(folder);
    verifyEqual(testCase, params.path, folder)
    verifyEqual(testCase, params.convGEM, fullfile(folder, 'models/x.xml'))
    verifyEqual(testCase, params.Ptot, 0.4)
    verifyEqual(testCase, params.sigma, 0.5)
    verifyEqual(testCase, params.kegg.ID, '')
    verifyEqual(testCase, params.uniprot.ID, '')
    verifyFalse(testCase, isfield(params, 'kcatAggregation'))
    verifyFalse(testCase, isfield(params, 'gR_exp'))
    verifyFalse(testCase, isfield(params.okp, 'method'))
    verifyEqual(testCase, params.evotune.kcatSources, {'custom','brenda','ec_tree','combined'})
    verifyEqual(testCase, params.evotune.sigma0logSource, [0.1; 0.2; 0.3; 0.4])
    verifyEqual(testCase, params.evotune.kcatGroups, params.evotune.kcatSources)
    verifyEqual(testCase, params.evotune.maxGenerations, 150)
    verifyEqual(testCase, params.okp.handleLongSequences, 'truncate')
    verifyEmpty(testCase, params.complex.taxonomicID)
    % A misspelled key is an error, as in geckopy.
    writeTemp({'conv_gem = "x.xml"', 'org_name = "x"', 'sigmaa = 1'}, fullfile(folder, 'model_adapter.toml'));
    verifyError(testCase, @() loadAdapterParams(folder), 'GECKO:loadAdapterParams')
end

function testLoadAdapterParamsSourceGroups_tc0065(testCase)
    % Trust tiers with several sources, a prefix pattern and a group width.
    folder = tempname; mkdir(folder);
    writeTemp({'conv_gem = "x.xml"', 'org_name = "x"', ...
        '[evotune.source_groups.okp]', 'sources = ["okp_*"]', 'sigma0_log = 0.25', ...
        '[evotune.source_groups.database]', 'sources = ["brenda", "okp_brenda"]', 'sigma0_log = 0.2', ...
        '[matlab]', 'make_model_anaerobic = "anaerobicModel_GECKO"', ...
        '[matlab.params]', 'c_UptakeExp = 0.64'}, fullfile(folder, 'model_adapter.toml'));
    params = loadAdapterParams(fullfile(folder, 'model_adapter.toml'));
    verifyEqual(testCase, params.evotune.kcatSources, {'okp_*','brenda','okp_brenda'})
    verifyEqual(testCase, params.evotune.kcatGroups, {'okp','database','database'})
    verifyEqual(testCase, params.evotune.sigma0logSource, [0.25; 0.2; 0.2])
    verifyEqual(testCase, params.matlab.makeModelAnaerobic, 'anaerobicModel_GECKO')
    verifyEqual(testCase, params.c_UptakeExp, 0.64)
end

function testConvertAdapterToToml_tc0066(testCase)
    % Converting an adapter class and loading the TOML gives the same
    % parameters back.
    geckoPath = findGECKOroot;
    legacyFolder = tempname; mkdir(legacyFolder);
    writeTemp(legacyAdapterClass(), fullfile(legacyFolder, 'LegacyFixtureAdapter.m'));
    adapters = {fullfile(geckoPath,'test','unit_tests','ecTestGEM','TestGEMAdapter.m'), ...
                fullfile(legacyFolder,'LegacyFixtureAdapter.m')};
    for k = 1:numel(adapters)
        original = ModelAdapterManager.getAdapter(adapters{k});
        folder = tempname; mkdir(folder);
        file = fullfile(folder, 'model_adapter.toml');
        report = convertAdapterToToml(original, file);
        verifyTrue(testCase, any(strcmp(report.methods, 'getSpontaneousReactions')))
        converted = TomlAdapter(file);
        expected = original.params;
        % Fields without a key in model_adapter.toml are not carried over.
        if k == 2
            verifyEqual(testCase, report.notCarried, ...
                {'kcatAggregation', 'evotune.kcatSources', 'evotune.sigma0logSource'})
            expected = rmfield(expected, 'kcatAggregation');
            expected.evotune = rmfield(expected.evotune, {'kcatSources', 'sigma0logSource'});
        else
            verifyEmpty(testCase, report.notCarried)
        end
        actual = converted.params;
        % The project folder is where the TOML file is; the model path is
        % relative to it.
        [~, convName] = fileparts(expected.convGEM);
        [~, convName2] = fileparts(actual.convGEM);
        verifyEqual(testCase, convName2, convName)
        expected = rmfield(expected, {'path','convGEM'});
        verifyEqual(testCase, keepFields(actual, expected), expected)
    end
    verifyEqual(testCase, report.unmapped, {'c_UptakeExp'})
end

function lines = legacyAdapterClass()
% An adapter class as written for GECKO 3 and 4: fields that
% model_adapter.toml has no key for (kcatAggregation, evotune trust tiers)
% and a parameter that geckopy does not know.
lines = { ...
    'classdef LegacyFixtureAdapter < ModelAdapter', ...
    '    methods', ...
    '        function obj = LegacyFixtureAdapter()', ...
    '            obj.params.path = fileparts(mfilename(''fullpath''));', ...
    '            obj.params.convGEM = fullfile(obj.params.path,''models'',''human-GEM.yml'');', ...
    '            obj.params.sigma = 0.1;', ...
    '            obj.params.Ptot = 0.505717472;', ...
    '            obj.params.f = 0.412;', ...
    '            obj.params.gR_exp = 0.020663429;', ...
    '            obj.params.c_UptakeExp = 0.641339301;', ...
    '            obj.params.org_name = ''homo sapiens'';', ...
    '            obj.params.kcatAggregation = ''median'';', ...
    '            obj.params.complex.taxonomicID = 9606;', ...
    '            obj.params.kegg.ID = ''hsa'';', ...
    '            obj.params.kegg.geneID = ''Ensembl'';', ...
    '            obj.params.uniprot.type = ''proteome'';', ...
    '            obj.params.uniprot.ID = ''UP000005640'';', ...
    '            obj.params.uniprot.geneIDfield = ''gene_primary'';', ...
    '            obj.params.uniprot.reviewed = true;', ...
    '            obj.params.c_source = ''MAR09034'';', ...
    '            obj.params.bioRxn = ''MAR13082'';', ...
    '            obj.params.enzyme_comp = ''Cytosol'';', ...
    '            obj.params.evotune.sigma0logDefault = 0.5;', ...
    '            obj.params.evotune.kcatSources = {''dlkcat'',''brenda'',''custom''};', ...
    '            obj.params.evotune.sigma0logSource = [0.4; 0.2; 0.1];', ...
    '            obj.params.evotune.maxGrowthWeight = 2.0;', ...
    '            obj.params.evotune.priorPenaltyWeight = 0.03;', ...
    '            obj.params.evotune.tieIsozymes = true;', ...
    '            obj.params.evotune.rmseThreshold = 0.2;', ...
    '            obj.params.evotune.maxGenerations = 212;', ...
    '        end', ...
    '        function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)', ...
    '            spont = contains(model.rxnNames,''spontaneous'');', ...
    '            spontRxnNames = model.rxnNames(spont);', ...
    '        end', ...
    '    end', ...
    'end'};
end

function testTomlAdapterSpontaneous_tc0067(testCase)
    folder = tempname; mkdir(fullfile(folder,'data'));
    writeTemp({'rxns	spontaneous', 'R1	0', 'R2	1', 'R3	1'}, fullfile(folder,'data','spont.tsv'));
    writeTemp({'conv_gem = "x.xml"', 'org_name = "x"', ...
        '[spontaneous]', 'name_contains = "spontaneous"', 'file = "data/spont.tsv"'}, ...
        fullfile(folder, 'model_adapter.toml'));
    adapter = ModelAdapterManager.getAdapter(folder);
    verifyClass(testCase, adapter, 'TomlAdapter')
    model.rxns = {'R1';'R2';'R4';'R5'};
    model.rxnNames = {'a';'b';'spontaneous c';'d'};
    [spont, names] = adapter.getSpontaneousReactions(model);
    verifyEqual(testCase, spont, [false; true; true; false])
    verifyEqual(testCase, names, {'R2';'R4'})
end

function testCreateGECKOprojectToml_tc0068(testCase)
    % A new project gets model_adapter.toml, which loads as it is.
    here = pwd;
    restore = onCleanup(@() cd(here));
    root = tempname; mkdir(root);
    createGECKOproject('name', 'ecTestProject', 'path', root);
    folder = fullfile(root, 'ecTestProject');
    verifyTrue(testCase, isfile(fullfile(folder, 'model_adapter.toml')))
    verifyEqual(testCase, strtrim(fileread(fullfile(folder, '.gitignore'))), 'data/okpApiKey.txt')
    adapter = ModelAdapterManager.getAdapter(folder);
    verifyEqual(testCase, adapter.params.org_name, 'genus species')
    verifyEqual(testCase, adapter.params.evotune.kcatSources, {'custom','brenda','ec_tree','combined'})
    verifyEqual(testCase, adapter.params.evotune.sigma0logSource, [0.1; 0.2; 0.3; 0.4])
    verifyFalse(testCase, isfield(adapter.params, 'gR_exp'))
    verifyEqual(testCase, adapter.params.spontaneous.name_contains, 'spontaneous')
end

function testKcatsFromOKPDefaultMethod_tc0737(testCase)
    % Without a method, kcatsFromOKP uses the first of [kcat_sources.okp] methods.
    adapter = ModelAdapterManager.getAdapter(fullfile(findGECKOroot,'test','unit_tests','ecTestGEM', 'TestGEMAdapter.m'));
    scratchDir = tempname;
    mkdir(fullfile(scratchDir, 'data'));
    cleanupDir = onCleanup(@() rmdir(scratchDir, 's'));
    adapter.params.path = scratchDir;
    ecModel.ec.rxns = {'R1'};
    try
        kcatsFromOKP(ecModel, 'useStored', true, 'save', false, 'modelAdapter', adapter);
        err = [];
    catch err
    end
    verifyEqual(testCase, err.identifier, 'kcatsFromOKP:noStored')
    verifySubstring(testCase, err.message, 'KinForm-H')
    adapter.params.kcat_sources.okp.methods = {'UniKP', 'CataPro'};
    try
        kcatsFromOKP(ecModel, 'useStored', true, 'save', false, 'modelAdapter', adapter);
        err = [];
    catch err
    end
    verifyEqual(testCase, err.identifier, 'kcatsFromOKP:noStored')
    verifySubstring(testCase, err.message, 'UniKP')
end

function testOkpMethodsAndFiles_tc0738(testCase)
    % CatPred and TurNup are not accepted in [kcat_sources.okp] methods, and
    % an okp_<x> file that is no known method's file fails naming the file.
    verifyError(testCase, @() kcatAssignmentConfig(struct('kcat_sources', ...
        struct('okp', struct('methods', {{'CatPred'}})))), 'kcatAssignment:config')
    known = 'CataPro, CatRange, DLKcat, EITLEM, KinForm-H, KinForm-L, OmniESI, OmniESI-O2DENet, UniKP';
    try
        kcatAssignmentConfig(struct('kcat_sources', struct('okp', struct('methods', {{'Zed', 'UniKP', 'Zed', 'Abc'}}))));
        err = [];
    catch err
    end
    verifyEqual(testCase, err.message, ['Unknown OKP methods in [kcat_sources.okp] methods: Zed, Abc; methods: ' known '.'])
    verifyError(testCase, @() kcatAssignmentConfig(struct('kcat_sources', ...
        struct('okp', struct('methods', {{'TurNup'}})))), 'kcatAssignment:config')
    folder = tempname; mkdir(folder);
    cleanupDir = onCleanup(@() rmdir(folder, 's'));
    writeTemp({'conv_gem = "x.xml"', 'org_name = "x"', '[[kcat_assignment.stage]]', ...
        'sources = ["okp_nosuchmethod", "okp_kinform_h"]'}, fullfile(folder, 'model_adapter.toml'));
    adapter = ModelAdapterManager.getAdapter(folder);
    report = evalc('r = updateKcatSources(struct(), ''modelAdapter'', adapter);'); %#ok<NASGU>
    row = strcmp(r.source, 'okp_nosuchmethod');
    verifyEqual(testCase, r.status{row}, 'failed')
    verifyEqual(testCase, r.note{row}, ['okp_nosuchmethod is not the file of a known OKP method; the files are ' ...
        'okp_catapro, okp_catrange, okp_dlkcat, okp_eitlem, okp_kinform_h, okp_kinform_l, okp_omniesi, ' ...
        'okp_omniesi_o2denet, okp_unikp.'])
end

function testOkpMethodsAtLoad_tc0744(testCase)
    % The methods are checked when the file is read.
    base = {'conv_gem = "x.xml"', 'org_name = "x"', '[kcat_sources.okp]'};
    known = 'CataPro, CatRange, DLKcat, EITLEM, KinForm-H, KinForm-L, OmniESI, OmniESI-O2DENet, UniKP';
    cases = { ...
        {'methods = "KinForm-H"'}, 'kcat_sources.okp.methods must be a list of strings'; ...
        {'methods = [1, 2]'}, 'kcat_sources.okp.methods must be a list of strings'; ...
        {'methods = ["Zed", "UniKP", "Zed", "CatPred"]'}, ...
            ['Unknown OKP methods in [kcat_sources.okp] methods: Zed, CatPred; methods: ' known '.']};
    for k = 1:size(cases, 1)
        folder = tempname; mkdir(folder);
        file = fullfile(folder, 'model_adapter.toml');
        writeTemp([base, cases{k,1}], file);
        try
            loadAdapterParams(file);
            err = [];
        catch err
        end
        verifyNotEmpty(testCase, err)
        if k < 3
            verifyEqual(testCase, err.message, [file ': ' cases{k,2}])
        else
            verifyEqual(testCase, err.message, cases{k,2})
        end
    end
end

function s = keepFields(s, template)
% s reduced to the fields (recursively) that template has.
fields = fieldnames(s);
for k = 1:numel(fields)
    if ~isfield(template, fields{k})
        s = rmfield(s, fields{k});
    elseif isstruct(s.(fields{k})) && isstruct(template.(fields{k}))
        s.(fields{k}) = keepFields(s.(fields{k}), template.(fields{k}));
    end
end
s = orderfields(s, template);
end

function file = writeTemp(lines, file)
if nargin < 2
    file = [tempname '.toml'];
end
fid = fopen(file, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
