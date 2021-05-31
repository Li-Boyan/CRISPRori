function saveLineageData(inputFile, outputFile, varargin)
load DefaultAnalysisParams.mat
% Customized parameters
for i = 1:2:numel(varargin)
    try eval(strcat(varargin{i},'=',varargin{i+1},';'))
    catch
        eval(strcat(convertCharsToStrings(varargin{i}), "='", ...
                    convertCharsToStrings(varargin{i+1}), "'",';'));
    end
end
try
    load(fullfile(inputPath, inputFile), 'lineage')
catch ME
    if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
        error([inputFile, 'could not be read. Continue to next one.'])
    end
end
lineage.timeCourse(...
    (lineage.timeCourse.Time==0) | ...
    (lineage.timeCourse.Time < tstart) |...
    (lineage.timeCourse.Time > tend), :) = [];
subplot(311)
plot(lineage.timeCourse.Time, lineage.timeCourse.Total_proteins / 3e6);
subplot(312)
plot(lineage.timeCourse.Time, lineage.timeCourse.DnaAatp_free);
subplot(313)
plot(lineage.timeCourse.Time, lineage.timeCourse.FtsZ);
if saveLineage, writetable(lineage.timeCourse, fullfile(outputPath, outputFile)); end
% 6, 7