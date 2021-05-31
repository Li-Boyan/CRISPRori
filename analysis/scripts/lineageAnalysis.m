function lineageAnalysis(varargin)
%LINEAGEANALYSIS is a batch processing program for simulation results from
% IntegratedSimulation. The processing results of this program include:
%       (1) Metadata that specifies the parameters used in each group;
%       (2) Instantaneous growth rate;
%       (3) Division cycle and cell volume;
%       (4) Division adder;
%       (5) Nucleoid-to-cytoplasm ratio;
%       (6) GFP proteomic fraction;
%       (7) Initiation mass;
%       (8) Cell volume at certain time point;
%       (9) Instantaneous GFP production rate;
%
% Parameters for analysis are as follows:

%       group: number of groups in total for memory pre-allocation; Default
%              is 50;

%       inputPath: Result file path. Default is ../../results.

%       matFileName: Pattern of file name. Default is *.mat.

%       outputPath: Path to put the output analysis results. Default is
%                  ../results.
%       
%       outputVarNames: File names for 9 output results. Default is :
%                         {
%                          'InstGrowthData',
%                          'DivisionCycleData',
%                          'DivisionLenData',
%                          'AdderData',
%                          'NCRatio',
%                          'GFP',
%                          'InitiationMassData',
%                          'CellSizeData',
%                          'GFP_rate'                 
%                         }
%
%       tstart: Start time point to record. Default is 180 minutes.
%
%       tend: Last time point to record. Default is 600 minutes.

% Default parameters
load DefaultAnalysisParams.mat
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
paramGroup = cell(1, group);
nGroup = 1;
lineageData = cell(9, length(fileNames));
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
    
    % Remove empty lines
    lineage.timeCourse((lineage.timeCourse.Time==0) | (lineage.timeCourse.Time < tstart) | (lineage.timeCourse.Time > tend), :) = [];
    lineage.di((lineage.di.Time==0) | (lineage.di.Time < tstart) | (lineage.di.Time > tend), :) = [];
    lineage.ri((lineage.ri.Time==0) | (lineage.ri.Time < tstart) | (lineage.ri.Time > tend), :) = [];

    if isempty(lineage.timeCourse) || lineage.timeCourse.Time(end) < 480, continue; end    
    
    lineage.params = rmfield(lineage.params, 'initCondition');
    
    % Specify the group of the input file
    if fileidx == 1
        paramGroup{1} = lineage.params;
        currGroup = 1;
        nGroup = nGroup + 1;
    else
        isExistParams = false;
        for j = 1:numel(paramGroup)
            if isequal(lineage.params, paramGroup{j})
                isExistParams = true;
                currGroup = j;
            end
        end
        if ~isExistParams
            paramGroup{nGroup} = lineage.params;
            nGroup = nGroup + 1;
            currGroup = nGroup; 
        end
    end
    
    % Instantaneous growth
    cellvol = (lineage.timeCourse.Total_proteins - lineage.timeCourse.GFP)/3e6;
    dvdt = diff(cellvol) ./ diff(lineage.timeCourse.Time);
    v = mean([cellvol(1:(end-1)), cellvol(2:end)], 2);
    time = mean([lineage.timeCourse.Time(1:(end-1)), lineage.timeCourse.Time(2:end)], 2);
    dnaContent = lineage.timeCourse.GeneNum ./ lineage.timeCourse.Total_proteins * 3e6/4536;
    dnaContent = mean([dnaContent(1:(end-1)), dnaContent(2:end)], 2);
    %   Remove division timepoint
    v(dvdt < 0) = [];
    dnaContent(dvdt < 0) = [];
    time(dvdt < 0) = [];
    dvdt(dvdt < 0) = [];
    group = repmat(currGroup, length(v), 1);
    srcfile = repmat(fileidx, length(v), 1);
    lineageData{1, fileidx} = table(v, dvdt, dnaContent, time, group, srcfile);
    
    % Division cycle and cell volume
    dc = diff(lineage.di.Time);
    Ld = lineage.di.Total_proteins/3e6;
    group = repmat(currGroup, length(dc), 1);
    srcfile = repmat(fileidx, length(dc), 1);
    lineageData{2, fileidx} = table(dc, group, srcfile);
    group = repmat(currGroup, length(Ld), 1);
    srcfile = repmat(fileidx, length(Ld), 1);    
    lineageData{3, fileidx} = table(Ld, group, srcfile);  
    
    % Adder
    newBornTimeIdx = find([false; (cellvol(1:(end-1)) > cellvol(2:end))]);
    divTimeIdx = newBornTimeIdx(2:end) - 1;
    if ~isempty(newBornTimeIdx)
        newBornTimeIdx(end) = [];        
        Lb = cellvol(newBornTimeIdx);
        divLen = cellvol(divTimeIdx);
        adder = divLen - Lb;
        group = repmat(currGroup, length(adder), 1);
        srcfile = repmat(fileidx, length(adder), 1);
        dc = lineage.timeCourse.Time(divTimeIdx) - lineage.timeCourse.Time(newBornTimeIdx);
    else 
        Lb = nan; adder = nan; dc=nan; group = currGroup; srcfile = fileidx;
    end
    lineageData{4, fileidx} = table(Lb, adder, dc, group, srcfile);
    
    % Nucleiod-cytoplasm ratio
    ncr = lineage.timeCourse.GenomeNum(selected_tp) ./ ...
          lineage.timeCourse.Total_proteins(selected_tp) * 2e6;
    nGenome = lineage.timeCourse.GenomeNum(selected_tp);
    group = currGroup;
    srcfile = fileidx;
    lineageData{5, fileidx} = table(ncr, nGenome, group, srcfile);    
    
    % GFP proteomic fraction
    gfp = lineage.timeCourse.GFP(selected_tp2) / lineage.timeCourse.Total_proteins(selected_tp2);
    group = currGroup;
    srcfile = fileidx;
    lineageData{6, fileidx} = table(gfp, group, srcfile);    
    
    % Initiation mass
    mi = lineage.ri.Total_proteins / 3e6 ./ lineage.ri.OriCNum;
    time = lineage.ri.Time;
    group = repmat(currGroup, length(mi), 1);
    srcfile = repmat(fileidx, length(mi), 1);
    lineageData{7, fileidx} = table(mi, time, group, srcfile);
    
    % Cell volume at certain timepoint
    v_tp = lineage.timeCourse.Total_proteins(selected_tp3) / 3e6;
    rv_tp = lineage.timeCourse.Ribosome(selected_tp3);
    group = currGroup;
    srcfile = fileidx;
    lineageData{8, fileidx} = table(v_tp, rv_tp, group, srcfile);    
    
    % Instantaneous GFP production rate
    cellvol = (lineage.timeCourse.Total_proteins - lineage.timeCourse.GFP)/3e6;
    gfp_prod = diff(lineage.timeCourse.GFP) ./ diff(lineage.timeCourse.Time);
    v = mean([cellvol(1:(end-1)), cellvol(2:end)], 2);
    dnaContent = lineage.timeCourse.GeneNum ./ lineage.timeCourse.Total_proteins * 3e6/4536;
    dnaContent = mean([dnaContent(1:(end-1)), dnaContent(2:end)], 2);
    %   Remove division timepoint
    v(gfp_prod < 0) = [];
    dnaContent(gfp_prod < 0) = [];
    gfp_prod(gfp_prod < 0) = [];
    gfp_prod = gfp_prod ./ v;
    group = repmat(currGroup, length(v), 1);
    srcfile = repmat(fileidx, length(v), 1);
    lineageData{9, fileidx} = table(dnaContent, gfp_prod, group, srcfile);
    
end
paramGroup(cellfun(@isempty, paramGroup)) = [];

% Saving results to files
if ~exist(outputPath, 'dir'), mkdir(outputPath); end
for i = 1:length(outputVarNames)
    writetable(vertcat(lineageData{i,:}), fullfile(outputPath, [outputVarNames{i}, '.csv']));
end

% Saving metadata about the groups
if length(paramGroup) > 1
parameterNames = {'nu_Rb', 'nu_R', 'nu_DnaA', 'nu_FtsZ', 'nu_total', 'Kt', 'alphaf', 'KDnaA', ...
                 'nDnaA', 'taum_Rb', 'taum_R', 'taum_DnaA', 'taum_FtsZ', 'taum_total',   ...
                 'sigma_Rb', 'sigma_R', 'sigma_DnaA', 'sigma_FtsZ', 'sigma_total', 've', ...
                 'a_inact', 'Ke', 'katp', 'kadp', 'ratp2adp', 'taup_FtsZ', 'kRIDA', ...
                 'kdatA', 'khbox', 'klbox', 'kon', 'koff', 'nu_GFP', 'taum_GFP', 'sigma_GFP'};
conditionNames = {'c_dCas9', 'dCas9_target', 'plasmidCopyNum', 'model'};
boxNames = {'R1', 'R5M', 'TAU2', 'I1', 'I2', 'C3', 'C2', 'I3', 'C1', 'R4'};
paramCmpComb = nchoosek(1:length(paramGroup), 2);
paramDiffIdx = cell(1, size(paramCmpComb, 1));
conditionDiffIdx = cell(1, size(paramCmpComb, 1));

% Check each parameter pair
for i = 1:size(paramCmpComb, 1)
    paramComb = paramCmpComb(i,:);
    lineageParam1 = cell2mat(paramGroup{paramComb(1)}.parameters);
    lineageParam2 = cell2mat(paramGroup{paramComb(2)}.parameters);
    lineageCond1  = paramGroup{paramComb(1)}.condition;
    lineageCond2  = paramGroup{paramComb(2)}.condition;
    lineageCond1  = [lineageCond1{1}, find(strcmp(boxNames, lineageCond1{2})), lineageCond1{3}, 2-(lineageCond1{4}=='OSHD')];
    lineageCond2  = [lineageCond2{1}, find(strcmp(boxNames, lineageCond2{2})), lineageCond2{3}, 2-(lineageCond2{4}=='OSHD')];
    paramDiffIdx{i} = find(lineageParam1 - lineageParam2);
    conditionDiffIdx{i} = find(lineageCond1 - lineageCond2);
end
paramDiffIdx = unique(horzcat(paramDiffIdx{:}));
conditionDiffIdx = unique(horzcat(conditionDiffIdx{:}));
paramArray = zeros(length(paramGroup), length(paramDiffIdx) + length(conditionDiffIdx));
for i = 1:length(paramGroup)
    lineageParam = paramGroup{i}.parameters;
    paramArray(i, 1:length(paramDiffIdx)) = cell2mat(lineageParam(paramDiffIdx));
    lineageCond  = paramGroup{i}.condition;
    lineageCond  = [lineageCond{1}, find(strcmp(boxNames, lineageCond{2})), lineageCond{3}, 2-(lineageCond{4}=='OSHD')];
    paramArray(i, (length(paramDiffIdx)+1):end) = lineageCond(conditionDiffIdx);
end
metadata = array2table([(1:length(paramGroup))', paramArray], ...
                       'VariableNames', [{'Group'}, parameterNames(paramDiffIdx), conditionNames(conditionDiffIdx)]);
writetable(metadata, fullfile(outputPath, [metaFileName, '.csv']))
end
