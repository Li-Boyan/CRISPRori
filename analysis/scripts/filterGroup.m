function targetFileNames = filterGroup(varargin)
%FILTERGROUP return the simulation file names with given filtering conditions.
%
% Default parameters
load DefaultAnalysisParams.mat 
load ../../scripts/DefaultParams.mat
% Customized parameters
for i = 1:2:numel(varargin)
    try eval(strcat(varargin{i},'=',varargin{i+1},';'))
    catch
        eval(strcat(convertCharsToStrings(varargin{i}), "='", ...
                    convertCharsToStrings(varargin{i+1}), "'",';'));
    end
end
% Listing names of all files to read
if exist(inputPath, 'dir')
    if recursiveReading, files = dir(fullfile(inputPath,  '**', matFileName));
    else, files = dir(fullfile(inputPath, matFileName));
    end
    fileNames = {files.name}';
else, error('Input Path does not exist!')
end
targetFileNames = cell(1,3e2);
i = 1;
for fileidx = 1:length(fileNames)
    %
    % Load data
    try
        load(fullfile(inputPath, fileNames{fileidx}), 'lineage')
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            warning([fileNames{fileidx}, 'could not be read. Continue to next one.'])
        end
        continue
    end
    if lineage.params.parameters{end-3} == koff ...
            && strcmp(lineage.params.condition{2}, dCas9_target) ...
            && lineage.params.condition{3} == plasmidCopyNum ...
            && lineage.params.condition{1} == c_dCas9 ...
            && lineage.params.parameters{end-8} == kRIDA ...
            && lineage.params.parameters{end-7} == kdatA ...
            && lineage.params.parameters{end-15} == ve ...
            && lineage.params.parameters{9} == nDnaA
        
        targetFileNames{i} = fileNames{fileidx};
        i = i + 1;
    end
end
targetFileNames(cellfun(@isempty, targetFileNames)) = [];
