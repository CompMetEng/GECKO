%run this test case with the command
%results = runtests('brendaMeasurementsXzTests.m')

function tests = brendaMeasurementsXzTests
    tests = functiontests(localfunctions);
end

function testLoadBrendaMeasurementsReadsXz_tc0403(testCase)
    [status, ~] = system('xz --version');
    assumeEqual(testCase, status, 0)
    folder = tempname();
    mkdir(folder);
    cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>
    writeMeasurements(fullfile(folder, 'measurements.tsv'), {'1.1.1.1', 5})
    % only the plain file
    m = loadBrendaMeasurements('folder', folder);
    verifyEqual(testCase, m.kcat, 5)
    % only the compressed file
    verifyEqual(testCase, system(sprintf('xz -k "%s"', fullfile(folder, 'measurements.tsv'))), 0)
    delete(fullfile(folder, 'measurements.tsv'))
    m = loadBrendaMeasurements('folder', folder);
    verifyEqual(testCase, m.kcat, 5)
    verifyEqual(testCase, m.release, 'BRENDA release 2026.1')
    % both files: the compressed one is read
    writeMeasurements(fullfile(folder, 'measurements.tsv'), {'1.1.1.1', 7})
    m = loadBrendaMeasurements('folder', folder);
    verifyEqual(testCase, m.kcat, 5)
    % neither
    delete(fullfile(folder, 'measurements.tsv'))
    delete(fullfile(folder, 'measurements.tsv.xz'))
    verifyError(testCase, @() loadBrendaMeasurements('folder', folder), 'loadBrendaMeasurements:missing')
end

function writeMeasurements(file, rows)
    fid = fopen(file, 'w');
    fprintf(fid, '# BRENDA release 2026.1\n');
    fprintf(fid, 'ec\tsubstrate\tsubstrate_inchikey\torganism\ttaxid\tuniprot\tseq_id\tkcat\ttemperature\tph\tpmid\n');
    for i = 1:size(rows, 1)
        fprintf(fid, '%s\tS\tK\tO\t1\tP1\t1\t%g\t25\t7\t1\n', rows{i, 1}, rows{i, 2});
    end
    fclose(fid);
end
