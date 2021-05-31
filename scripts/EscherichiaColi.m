classdef EscherichiaColi < handle

%ESCHERICHIACOLI creates an integrative simulation platform of Escherichia coli
% which jointly took consideration of protein expression, growth and cell cycle events.
%
properties
%
% A table array recording evolution of all substances in the cell.
% The title line is the name string of all substances and the
% following lines are quantities of substances at different
% timepoints, where number of records is user-defined
substanceData
%
% Current timepoint after starting simulation
timePoint
%
% A cell array of which elements are single genomes. See class
% BACTERIALGENOME
genomes
%
% Sealing all reactions into a system which can perform evolution
% itself based on chosen method.
% See class LANGEVIN.        
reactionSyst
%
% Threshold value of intracellular number of FtsZ to trigger cell 
% division.
FtsZ_thres
%
% Target box of dCas9: high or low-affinity box (1=HIGH, 2=LOW)
dCasTarget
%
% Method of modelling
modellingMethod
end

properties(Dependent)

% Total gene copy number
totGeneCopyNum

% OriC number free for dCas9 to bind in each genome
availableOriC

% Number of oriCs that have been occupied by dCas9 in each genome.
dCasBoundOriC

% Total number of replication forks.
repForkNum

% Total number of OriC
nOriC

% Number of available datA loci on each genome
datA

% Total number of all high-affinity boxes, all low-affinity boxes,
% all free high-affinity boxes and all free low-affinity boxes,
% whose data type is CONTAINERS.MAP
boxNum
end

properties(Constant)  

% Genomic location of different genes or loci.
LocRb    = 0.7233; % Gene rpsA
LocR     = 0.1127; % Gene rpoC
LocDnaA  = 0.0185;
LocFtsZ  = 0.3540;   
LocDARS1 = 0.6584;
LocDARS2 = 0.4130;
LocDatA  = 0.2010;

% Protein density, in the unit of 1/um^3
proDensity = 3e6;

% DnaA boxes
DnaABoxes = {'R1','R5M','TAU1','I1','I2','R4','C1','I3','C2','C3'};
end

methods
function obj = EscherichiaColi(initialCondition, genomeInitCondition, ...
                               recordLineNum, params, FtsZ_thres, ...
                               c_dCas9, targetBox, plasmidNum, modellingMethod)                                   
    %ESCHERICHIACOLI creates an object of class ESCHERICHIACOLI.

    obj.timePoint = 0;
    obj.FtsZ_thres = FtsZ_thres;
    obj.modellingMethod = modellingMethod;

    % Define the initial state of the genome (see class
    % BACTERIALGENOME)
    [sqstt, cperiod, init] = genomeInitCondition{:};
    obj.genomes = {bacterialGenome(sqstt, cperiod, init, obj.modellingMethod)};

    % Define all the parameters used in kinetic equations (see below).
    [nu_Rb, nu_R, nu_DnaA, nu_FtsZ, nu_total, Kt, alphaf, KDnaA, ...
     nDnaA, taum_Rb, taum_R, taum_DnaA, taum_FtsZ, taum_total,   ...
     sigma_Rb, sigma_R, sigma_DnaA, sigma_FtsZ, sigma_total, ve, ...
     a_inact, Ke, katp, kadp, ratp2adp, taup_FtsZ, kRIDA, ...
     kdatA, khbox, klbox,kon,koff,nu_GFP,taum_GFP,sigma_GFP] = params{:};

    % Define which box for dCas9 to target.            
    if strcmp(obj.modellingMethod, 'DRF')

        switch targetBox
        case 'high', obj.dCasTarget = 1;
        case 'low',  obj.dCasTarget = 2;
        otherwise, error('Unknown box type! Select either high or low');
        end

    elseif strcmp(obj.modellingMethod, 'OSHD')
        obj.dCasTarget = find(ismember(obj.DnaABoxes,targetBox));
        if isempty(obj.dCasTarget), error('Unknown DnaA Box!'); end
    else
        error('Unknown modelling method! Select either DRF or OSHD.')
    end

    % Define all reactions and seal them into LANGEVIN class.
    % Reactions in the order:
    % 01. rpsA mRNA transcription;
    % 02. rpoC mRNA transcription;
    % 03. dnaA mRNA transcription;
    % 04. ftsZ mRNA transcription;
    % 05. Total transcription;
    % 06. rpsA mRNA degradation;
    % 07. rpoC mRNA degradation;
    % 08. dnaA mRNA degradation;
    % 09. ftsZ mRNA degradation;
    % 10. Total degradation;
    % 11. RpsA translation;
    % 12. RpoC translation;
    % 13. DnaA in ATP-bound form production;
    % 14. DnaA in ADP-bound form production;
    % 15. FtsZ translation;
    % 16. Total translation;
    % 17. FtsZ degradation;
    % 18. Replicative inactivation of DnaA;
    % 19. datA titration of DnaA-ATP;
    % 20. DnaA-ATP binding to high-affinity boxes on oriC;
    % 21. DnaA-ADP binding to low-affinity boxes on oriC;
    % 23. dCas9 binding to oriC;
    % 24. dCas9 dissociation from oriC;
    % 25. gfp mRNA transcription;
    % 26. GFP translation.

    obj.reactionSyst = Langevin({ ...
    Reaction({'Ribosome_mRNA',...
              'Total_mRNA'   },   [1,1],    {'RNAP'         },  @(x) nu_Rb    * obj.geneCopyNum(obj.LocRb) * x / (x + obj.totGeneCopyNum * Kt / alphaf)), ...
    Reaction({'RNAP_mRNA',...
              'Total_mRNA'   },   [1,1],    {'RNAP'         },  @(x) nu_R     * obj.geneCopyNum(obj.LocR) * x / (x + obj.totGeneCopyNum * Kt / alphaf)), ...
    Reaction({'DnaA_mRNA',...
              'Total_mRNA'   },   [1,1],    {'RNAP',...
                                             'DnaAatp_free' },  @(x) nu_DnaA  * obj.geneCopyNum(obj.LocDnaA) * x(1) / (x(1) + obj.totGeneCopyNum * Kt / alphaf)  / (KDnaA^nDnaA + x(2)^nDnaA)), ...
    Reaction({'FtsZ_mRNA',...
              'Total_mRNA'   },   [1,1],    {'RNAP'         },  @(x) nu_FtsZ  * obj.geneCopyNum(obj.LocFtsZ) * x / (x + obj.totGeneCopyNum * Kt / alphaf)), ...
    Reaction({'Total_mRNA'   },   1,        {'RNAP'         },  @(x) nu_total * obj.totGeneCopyNum * x / (x + obj.totGeneCopyNum * Kt / alphaf)), ...
    Reaction({'Ribosome_mRNA',...
              'Total_mRNA'   },  [-1,-1],   {'Ribosome_mRNA'},  @(x) x / taum_Rb), ...
    Reaction({'RNAP_mRNA' ,...
              'Total_mRNA'   },  [-1,-1],   {'RNAP_mRNA'    },  @(x) x / taum_R), ...
    Reaction({'DnaA_mRNA',...
              'Total_mRNA'   },  [-1,-1],   {'DnaA_mRNA'    },  @(x) x / taum_DnaA), ...
    Reaction({'FtsZ_mRNA',...
              'Total_mRNA'   },  [-1,-1],   {'FtsZ_mRNA'    },  @(x) x / taum_FtsZ), ...
    Reaction({'Total_mRNA'   },  -1,        {'Total_mRNA'   },  @(x) x / taum_total), ...                     
    Reaction({'Ribosome',...
             'Total_proteins'},  [1,1],     {'Ribosome_mRNA',...
                                             'Ribosome',...
                                             'Total_mRNA',...
                                             'Total_proteins'}, @(x) sigma_Rb   * ve * x(1) * (x(2) - x(4) * a_inact) / ((x(2) - x(4) * a_inact) + Ke * x(3))), ...
    Reaction({'RNAP',...
             'Total_proteins'},  [1,1],     {'RNAP_mRNA',...
                                          'Ribosome',...
                                          'Total_mRNA',...
                                          'Total_proteins'}, @(x) sigma_R    * ve * x(1) * (x(2) - x(4) * a_inact) / ((x(2) - x(4) * a_inact) + Ke * x(3))), ...
    Reaction({'DnaA',...
              'DnaAatp',...
              'DnaAatp_free',...
             'Total_proteins'}, [1,1,1,1],  {'DnaA_mRNA',...
                                          'Ribosome',...
                                          'Total_mRNA',...
                                         'Total_proteins'},  @(x) sigma_DnaA * ve * x(1) * (x(2) - x(4) * a_inact) / ((x(2) - x(4) * a_inact) + Ke * x(3)) * katp / (kadp + katp * ratp2adp)), ...
    Reaction({'DnaA',...
              'DnaAadp',...
             'Total_proteins'},  [1,1,1],   {'DnaA_mRNA',...
                                          'Ribosome',...
                                          'Total_mRNA',...
                                         'Total_proteins'},  @(x) sigma_DnaA * ve * x(1) * (x(2) - x(4) * a_inact) / ((x(2) - x(4) * a_inact) + Ke * x(3)) * kadp / (kadp + katp*ratp2adp)), ...
    Reaction({'FtsZ',...
             'Total_proteins'},  [1,1],     {'FtsZ_mRNA',...
                                          'Ribosome',...
                                          'Total_mRNA',...
                                         'Total_proteins'},  @(x) sigma_FtsZ * ve * x(1) * (x(2) - x(4) * a_inact) / ((x(2) - x(4) * a_inact) + Ke * x(3))), ...
    Reaction({'Total_proteins'},  1,        {'Ribosome',...
                                          'Total_mRNA',...
                                         'Total_proteins'},  @(x) sigma_total * ve * x(2)  * (x(1) - x(3) * a_inact) / ((x(1) - x(3) * a_inact) + Ke * x(2))), ...
    Reaction({'FtsZ',...
             'Total_proteins'},  [-1,-1],   {'FtsZ'      },  @(x) x / taup_FtsZ), ...
    Reaction({'DnaAatp', ...
              'DnaAadp'      }, [-1,1],     {'DnaAatp',...
                                           'DnaAatp_free'},  @(x) kRIDA * obj.repForkNum * (x(1)-x(2)) * (x(1) > x(2))), ...
    Reaction({'DnaAatp_free' },  -1,        {'DnaAatp_free',...
                                         'Total_proteins'},  @(x) kdatA * sum(obj.datA) * x(1) * 1.661 / (x(2) / obj.proDensity)),  ...
    Reaction({'DnaAatp_free' },  -1,        {'DnaAatp_free',...
                                         'Total_proteins'},  @(x) khbox * obj.boxNum('hboxfree') * x(1) * 1.661 / (x(2) / obj.proDensity)),  ...
    Reaction({'DnaAatp_free' },  -1,        {'DnaAatp_free',...
                                         'Total_proteins'},  @(x) klbox * obj.boxNum('lboxfree') * x(1) * 1.661 / (x(2) / obj.proDensity)),  ...
    Reaction({               },   [],       {            },  @(x) kon * c_dCas9 * sum(obj.availableOriC) * (obj.timePoint >= 180)),    ...
    Reaction({               },   [],       {            },  @(x) koff * sum(obj.dCasBoundOriC)), ...
    Reaction({'GFP_mRNA',...
              'Total_mRNA'   },   [1,1],    {'RNAP', ...
                                         'Total_proteins'},  @(x) nu_GFP * floor(plasmidNum * x(2) / obj.proDensity) * x(1) / (x(1) + obj.totGeneCopyNum * Kt / alphaf)), ...
    Reaction({'GFP_mRNA' ,...
              'Total_mRNA'   },  [-1,-1],   {'GFP_mRNA'  },  @(x) x/taum_GFP),...
    Reaction({'GFP'   ,...
             'Total_proteins'},  [1,1],     {'GFP_mRNA',...
                                          'Ribosome',...
                                          'Total_mRNA',...
                                         'Total_proteins'}, @(x) sigma_GFP * ve * x(1) * (x(2) - x(4) * a_inact) / ((x(2) - x(4) * a_inact) + Ke * x(3)))}, initialCondition);

    % Record the initial state.
    obj.recordData(1,0,recordLineNum);

end             

function recordData(obj,line,label,recordLineNum)            
    %RECORDDATA Records information about the lineage into a table array.

    if nargin < 3, recordLineNum = []; end

    appendData = [obj.timePoint, ...
                  obj.totGeneCopyNum, ...
                  length(obj.genomes), ...
                  obj.reactionSyst.currSubstances('Ribosome') / obj.reactionSyst.currSubstances('Total_proteins'), ...
                  sum(obj.dCasBoundOriC), ...
                  sum(obj.nOriC), ...
                  sum(obj.repForkNum),...
                  sum(obj.datA),...
                  obj.boxNum('hboxfree'),...
                  obj.boxNum('lboxfree'),...
                  label, ...
                  obj.totGeneCopyNum/obj.reactionSyst.currSubstances('Total_proteins')*obj.proDensity];
    appendVar  = {'Time', ...
                  'GeneNum', ...
                  'GenomeNum', ...
                  'RibosomeContent' ...
                  'dCasBound', ...
                  'OriCNum',...
                  'ForkNum',...
                  'datANum',...
                  'hbox',...
                  'lbox',...
                  'label',...
                  'DNAContent'};

    % Pre-define a table for data record   
    currLine = [array2table(cell2mat(obj.reactionSyst.currSubstances.values), 'VariableNames', obj.reactionSyst.currSubstances.keys), ...
                     array2table(appendData, 'VariableNames', appendVar)];

    if line == 1
        obj.substanceData = [currLine; ...
                             array2table(zeros(recordLineNum, width(currLine)),'VariableNames',currLine.Properties.VariableNames)];
    else, obj.substanceData(line,:) = currLine;
    end

end

function geneNum = get.totGeneCopyNum(obj)
    %TOTGENECOPYNUM Determines the total gene copy number.

    geneNum = 0;

    % Calculte gene copy number on every genome.
    for i = 1:length(obj.genomes)
        geneNum = geneNum + obj.genomes{i}.totGeneCopyNum;
    end

end

function geneNum = geneCopyNum(obj, locOfGene)
    %GENECOPYNUM calculates the total copy number of a gene with given location on the bacterial genome.

    geneNum = 0;

    for i = 1:length(obj.genomes)
        geneNum = geneNum + obj.genomes{i}.geneCopyNum(locOfGene); 
    end

end        

function y = get.boxNum(obj)
    %BOXNUM Returns all free high-affinity and low-affinity boxes on all genomes in the cell.

    hboxtot_all  = 0;
    lboxtot_all  = 0;
    hboxfree_all = 0;
    lboxfree_all = 0;

    for i = 1:length(obj.genomes)
        [hboxTot, lboxTot, hboxFree, lboxFree] = obj.genomes{i}.boxNum();
        hboxfree_all = hboxfree_all + sum(sum(hboxFree));
        lboxfree_all = lboxfree_all + sum(sum(lboxFree));
        hboxtot_all  = hboxtot_all + sum(sum(hboxTot));
        lboxtot_all  = lboxtot_all + sum(sum(lboxTot));
    end

    y = containers.Map({'hboxtot','lboxtot','hboxfree','lboxfree'}, ...
                       [hboxtot_all,lboxtot_all,hboxfree_all,lboxfree_all]);
end

function bn = specBoxNum(obj,k)
    %SPECBOXNUM Returns the total number of free kth boxes on all genomes in the cell.

    bn = 0;

    for i = 1:length(obj.genomes)
        bn = bn + sum(obj.genomes{i}.numOfSpecBox(k));
    end

end

function generateLineage(obj,dt,ttol,showGenomes,filename,exittime)
    %GENERATELINEAGE is the main simulation function.

    rcdIntv = ttol / (height(obj.substanceData) - 1);    
    % Recorded time points (uniformly distributed)    
    rcdCount = 2;   
    frameNo = 1;
    delayTime = 1e-1;

    % Maximize figure window
    if showGenomes, figure('units','normalized','outerposition',[0 0 1 1]); end

    runtime = 0;

    while obj.timePoint <= ttol                   
        tic

        % Update the substances quantity
        [~, reactionUpdate] = obj.reactionSyst.evolution(dt, false);
        hboxON = reactionUpdate(20);
        lboxON = reactionUpdate(21);
        dCasON = reactionUpdate(22);
        dCasOFF = reactionUpdate(23);
        datAON = reactionUpdate(19) - reactionUpdate(18);
        obj.timePoint = obj.timePoint + dt;

        % Pre-define the (possibly) newly-generated and replaced genomes
        div_genome = zeros(1,10);
        div_count = 1;
        new_genome = cell(1,20);

        % Number of bound high-affinity and low-affinity boxes
        hboxBound = zeros(1,length(obj.genomes));
        lboxBound = zeros(1,length(obj.genomes));
        dCasBound = zeros(1,length(obj.genomes));
        dCasDisso = zeros(1,length(obj.genomes));
        datABound = zeros(1,length(obj.genomes));

        for i = 1:length(obj.genomes)
            [~,~,hboxfree,lboxfree] = obj.genomes{i}.boxNum();
            if obj.boxNum('hboxfree') > 0
                hboxBound(i) = hboxON * sum(sum(hboxfree)) / obj.boxNum('hboxfree');
            end
            if obj.boxNum('lboxfree') > 0
                lboxBound(i) = lboxON * sum(sum(lboxfree)) / obj.boxNum('lboxfree');
            end
        end

        if sum(obj.availableOriC) > 0
            dCasBound = dCasON * obj.availableOriC ./ sum(obj.availableOriC);
        end

        if sum(obj.dCasBoundOriC) > 0
            dCasDisso = dCasOFF * obj.dCasBoundOriC / sum(obj.dCasBoundOriC);
        end

        if sum(obj.datA) > 0
            datABound = datAON * obj.datA / sum(obj.datA);
        end

        for i = 1:length(obj.genomes)
            % Allocate the bound boxes to different genomes.
            hboxBoundHere = hboxBound(i); 
            lboxBoundHere = lboxBound(i);
            dCasBoundHere = dCasBound(i);
            dCasDissoHere = dCasDisso(i);
            datABoundHere = datABound(i);
            obj.genomes{i}.dCasBinding(dCasBoundHere,obj.dCasTarget)
            obj.genomes{i}.dCasDissociation(dCasDissoHere,obj.dCasTarget)
            obj.genomes{i}.datABinding(datABoundHere)
            [nri_1,~,releasedDnaA_1] = obj.genomes{i}.oriCBinding(1,hboxBoundHere);
            [nri_2,~,releasedDnaA_2] = obj.genomes{i}.oriCBinding(2,lboxBoundHere);

            if nri_1 + nri_2 > 0, obj.recordData(rcdCount,1); rcdCount = rcdCount + 1; end

            obj.reactionSyst.currSubstances('DnaAatp_free') = obj.reactionSyst.currSubstances('DnaAatp_free') ...
                                                   + (releasedDnaA_1 + releasedDnaA_2);
            % Debug
            % if nri_1 + nri_2 > 0, disp([releasedDnaA_1, releasedDnaA_2]); end
            
            % Replication termination                    
            datA_before = obj.genomes{i}.geneCopyNum(obj.LocDatA);

            if obj.genomes{i}.replisomeAdvancing(dt)
                datA_after = obj.genomes{i}.geneCopyNum(obj.LocDatA); 
                obj.genomes{i}.datA = obj.genomes{i}.datA + ...
                    (datA_after - datA_before) * obj.genomes{i}.DnaAPerDatA * (datA_after >= datA_before);                        
                % Define two new genomes
                currGenome_1 = copy(obj.genomes{i});
                currGenome_2 = copy(obj.genomes{i});
                currGenome_1.nucleoidDivision(1)
                currGenome_2.nucleoidDivision(2)
                div_genome(div_count) = i;
                new_genome([2*div_count-1,2*div_count]) = {currGenome_1, currGenome_2};
                div_count = div_count + 1;
            end                    

        end

        % Delete the old genome and add two new ones
        obj.genomes(div_genome~=0) = [];
        obj.genomes = [obj.genomes, new_genome(1:(2*(div_count-1)))];
        obj.recordData(rcdCount,2);
        iscd = obj.cellDivision();  

        if iscd, rcdCount = rcdCount + 1; end

        % Record the data and plot the genome
        if rem(obj.timePoint, rcdIntv) < dt
            % fprintf('Time: %d; AvailableOri: %d', obj.timePoint, obj.availableOriC);
            obj.recordData(rcdCount,0)
            rcdCount = rcdCount + 1;

            if showGenomes
                % Genome(s)
                subplot(2,1,1), obj.plotGenomes(5), drawnow
                % Cell Volume
                subplot(6,1,4)
                curves_1 = plotTable(obj.substanceData, 'Time', {'Total_proteins'}, @(x) x/obj.proDensity);
                maxh = 0;

                for c = 1:length(curves_1)
                    currLine = curves_1{c};
                    currYData = get(currLine,'ydata');
                    maxh = maxh * (maxh > max(currYData)) + max(currYData) * (maxh <= max(currYData));
                    currLine.LineWidth = 2;
                    currLine.Color = 'k';
                end

                maxh = ceil((maxh)/10) * 10;
                hold on, plot([obj.timePoint, obj.timePoint], [0, maxh], 'LineWidth', 1, 'Color', 'k'), hold off
                axis([0,ttol,0,maxh])
                set(gca,'xtick',[])
                ylabel({'Cell volume';'(\mum^3)'},'fontsize',12)

                % Total Gene Copy Number
                subplot(6,1,5)
                curves_1 = plotTable(obj.substanceData, 'Time', {'GeneNum'}, @(x) x/obj.genomes{1}.Gnum);
                maxh = 0;

                for c = 1:length(curves_1)
                    currLine = curves_1{c};
                    currYData = get(currLine,'ydata');
                    maxh = maxh * (maxh > max(currYData)) + max(currYData) * (maxh <= max(currYData));
                    currLine.LineWidth = 2;
                    currLine.Color = 'k';
                end

                maxh = ceil((maxh));
                hold on, plot([obj.timePoint, obj.timePoint], [0, maxh], 'LineWidth', 1, 'Color', 'k'), hold off
                axis([0,ttol,0,maxh])
                set(gca,'xtick',[])
                ylabel({'Total Gene';'Copy Number'},'fontsize',12)

                % DNA concentration
                subplot(6,1,6)
                curves_1 = plotTable(obj.substanceData, 'Time', {'DNAContent'}, @(x) x/obj.genomes{1}.Gnum);
                maxh = 0;

                for c = 1:length(curves_1)
                    currLine = curves_1{c};
                    currYData = get(currLine,'ydata');
                    maxh = maxh * (maxh > max(currYData)) + max(currYData) * (maxh <= max(currYData));
                    currLine.LineWidth = 2;
                    currLine.Color = 'k';
                end

                maxh = ceil((maxh)/0.5) * 0.5;
                hold on, plot([obj.timePoint, obj.timePoint], [0, maxh], 'LineWidth', 1, 'Color', 'k'), hold off
                axis([0,ttol,0,maxh])     
                ylabel({'DNA';'Content'},'fontsize',12)
                xlabel('Time (min)','fontsize',12)

                % Record into .gif file
                currFrame = getframe(gcf);
                imind = frame2im(currFrame);
                [imind,map] = rgb2ind(imind, 256);

                if frameNo == 1
                imwrite(imind,map,filename,'gif','Loopcount',inf,'DelayTime',delayTime,'DisposalMethod','restoreBG')
                else
                imwrite(imind,map,filename,'gif','WriteMode','append','DelayTime',delayTime,'DisposalMethod','restoreBG')
                end

                frameNo = frameNo + 1;
                % Save all frames in the pdf format
                % currFig = gcf;
                % currFig.PaperPositionMode='auto';
                % pdfName = ['snapshot', num2str(rcdCount-1), '.pdf'];
                % print(pdfName,'-dpdf','-fillpage')
            end

        end        

        runtime = runtime + toc;

        if runtime > exittime, break; end

        % Update waitbar
        % wtstr = strcat(sprintf("Simulation completes: %d", ceil(obj.timePoint/ttol*100)), "%");
        % waitbar(obj.timePoint/ttol, wtbar, wtstr) 
    end

    % close(wtbar)            
end

function iscd = cellDivision(obj)
    %CEllDIVISION Judges if cell division occurs.

    if obj.reactionSyst.currSubstances('FtsZ') >= obj.FtsZ_thres && length(obj.genomes) > 1
        iscd = true;
        % Allocate genomes
        obj.genomes = obj.genomes(1:floor(length(obj.genomes)/2));

        % Allocate substances based on binomial distribution
        for i = 1:length(obj.reactionSyst.currSubstances.keys)
            all_substance = obj.reactionSyst.currSubstances.keys;
            substance = all_substance{i};

            if ~strcmp(substance,'hbox') && ~strcmp(substance,'lbox')
                obj.reactionSyst.currSubstances(substance) = binornd(floor(obj.reactionSyst.currSubstances(substance)),0.5);
            end

        end

    else, iscd = false;
    end

end

function y = get.availableOriC(obj)
    %AVAILABLEORIC Returns the number of target boxes on each genome in the cell.

    y = zeros(1,length(obj.genomes));

    for i = 1:length(obj.genomes)
        currGenome = obj.genomes{i};

        switch obj.modellingMethod
            case 'DRF'
                y(i) = sum(currGenome.oriCstate(1,:) == 1 ...
                            & currGenome.oriC(obj.dCasTarget,:) > 0 ...
                            & currGenome.oriC(obj.dCasTarget+2,:) == 0);
            case 'OSHD'
                y(i) = sum(currGenome.numOfSpecBox(obj.dCasTarget));
        end

    end

end

function y = get.nOriC(obj)
    %NORIC Returns total number of oriCs in the cell.

    y = zeros(1,length(obj.genomes));

    for i = 1:length(obj.genomes)
        currGenome = obj.genomes{i};
        y(i) = size(currGenome.oriC,2);
    end

end

function y = get.dCasBoundOriC(obj)
    %DCASBOUNDORIC Returns number of oriCs bound by dCas9 in the cell.

    y = zeros(1,length(obj.genomes));

    for i = 1:length(obj.genomes)
        currGenome = obj.genomes{i};

        switch obj.modellingMethod
            case 'DRF'
                y(i) = sum((currGenome.oriC(obj.dCasTarget+2,:)~=0));
            case 'OSHD'
                y(i) = sum((currGenome.oriC(ceil(obj.dCasTarget/5)+2,:)~=6));
        end

    end

end

function y = get.repForkNum(obj)
    %REPFORKNUM returns total number of replication forks.

    y = 0;

    for i = 1:length(obj.genomes)
        y = y + size(obj.genomes{i}.repFork,2);
    end

end

function y = get.datA(obj)
    %DATA returns all the empty datA boxes in the cell.

    y = zeros(1,length(obj.genomes));

    for i = 1:length(obj.genomes)
        y(i) = obj.genomes{i}.datA;
    end

end

function plotGenomes(obj,xIntv)
    %PLOTGENOMES Plots and aligns all genomes in one figure
    %   PLOTGENOMES(xIntv) plot all genomes with horizontal 
    %   equal to xIntv between adjacent ones.

    nGenomes = length(obj.genomes);
    maxOriNum = 0;
    sc = 0.05;

    for i = 1:nGenomes                
        if size(obj.genomes{i}.oriC,2) > maxOriNum
            maxOriNum = size(obj.genomes{i}.oriC,2);
        end
    end

    if rem(nGenomes,2) == 0
        r0_right = (xIntv/2) : xIntv : (-xIntv/2 + xIntv*nGenomes/2);
        r0 = [-fliplr(r0_right), r0_right] * (1 - sc*maxOriNum);
    else
        r0_right = xIntv : xIntv : (xIntv*floor(nGenomes/2));
        r0 = [-fliplr(r0_right), 0, r0_right] * (1 - sc*maxOriNum);
    end

    for i = 1:nGenomes
        obj.genomes{i}.plotGenome(1-rem(i,2),-sc*maxOriNum,0,r0(i)), hold on
    end

    hold off
end

end
    
end
