classdef bacterialGenome < matlab.mixin.Copyable
    
%BACTERIALGENOME creates aggregation of all information of genomes in a bacterial cell
% and exert functions on DNA replication-associated events.

properties

% Information about all the replication forks on the genome,
% represented by an 2*n double array, where n is the number of
% current replication forks:
% the first line is the location of the forks (>0);
% the second line is the serial number of the replication forks in
% the binary tree.
repFork

% Information about all the replication origins on the genome,
% represented by an 2*m double array, where m is the number of
% replication origins:
% the first line and the second line are number of empty
% high-affinity and low-affinity DnaA-boxes, respectively.
% The third and fourth line refers to dCas9 binding state, 
% 0: bound and 1:unbound, to high-affinity and low-affinity DnaA
% boxes, respectively.
oriC

% Sequestration period (min)
sqstTime

% C-period (min)
C

% Total number of available datA locus
datA

% Model used: DRF=1, OSHD=2
model
end

properties(Constant)

% Total number of genes on E. coli K-12 M1655 genome.
Gnum     = 4538;  

% DnaA titrated by datA
DnaAPerDatA = 200;

% Number of high-affinity boxes per OriC
hboxPerOri = 3;

% Number of low-affinity boxes per OriC
lboxPerOri = 8;
end

properties(Dependent)

% Replication forks in the last layer of the binary tree:
% the first line is the locations;
% the second line is the serial number;
% the third line is the number of corresponding replication
% origins, which can be either 1 or 2.        
lastLayerFork

% State of all replication origins:
% the first line is the sequestration state: 1=true, 0=false;
% the second line is the serial number of corresponding
% replication forks.        
oriCstate

% Total gene copy number
totGeneCopyNum

% Used in OSHD model, to calculate replication initiation
% probability in each oriC
riProb
end

methods
function obj = bacterialGenome(sqstt, cperiod, init, model)
    %BACTERIALGENOME is the contructor function for BACTERIALGENOME
    %
    % Set the initial values for all properties:
    % Sequestration period (Time of the region near oriC to be bound
    % by SeqA)
    obj.sqstTime = sqstt;
    % C period
    obj.C = cperiod;
    obj.repFork = init;
    obj.datA = 200;

    if strcmp(model,'DRF') % DnaA random filling model
        obj.model = 1;
        obj.oriC = repmat([obj.hboxPerOri;obj.lboxPerOri;0;0], 1, size(obj.oriCstate,2));
    elseif strcmp(model,'OSHD') % Oligomerization Started by High-affinity DnaA
        obj.model = 2;
        obj.oriC = repmat([5; 5; 6; 6], 1, size(obj.oriCstate,2));
    end

end

function p = get.riProb(obj)
    %RIPROB Returns the initiation probability of all replication origins, based on the OSHD model.

    occupyFreq = [5;5] - obj.oriC(1:2,:);
    box_score = [6.5, 1.0, 0.5, 0.5, 0.5; ...
                 0.4, 0.4, 0.4, 0.4, 0.4];  
    score = zeros(1,size(occupyFreq,2));
    exp_idx = 3;

    for i = 1:size(occupyFreq,2)
        score(i) = sum(box_score(1,1:occupyFreq(1,i))) ...
                 + sum(box_score(2,1:occupyFreq(2,i)));
    end

    p = (exp(score*exp_idx)-1) / (exp(exp_idx*10)-1);
    fullscore = sum(sum(box_score));
    fullp = (exp(fullscore*exp_idx)-1) / (exp(exp_idx*10)-1); 
    p(abs(p-fullp)>1e-10) = p(abs(p-fullp)>1e-10) / 50;
    p(abs(p-fullp)<=1e-10) = ones(1,sum(abs(p-fullp)<=1e-10));
end

function [isri, releasedDnaA] = replicationInitiation(obj)
    %REPLICATIONINITIATION Judges if a replication initiation event occurs.

    if obj.model == 1
        judgeBox = find(obj.oriC(1,:) <= 0 & obj.oriC(2,:)<=0);
    elseif obj.model == 2
        judgeBox = find(obj.riProb >= rand(1,size(obj.oriC,2)));
    end

    if ~isempty(judgeBox)
        isri = sum(judgeBox); 
        occupyFreq = [5;5] - obj.oriC(1:2,:);
        releasedDnaA = sum(sum(occupyFreq(:, judgeBox)));
        % Debug
        % disp(obj.oriC)                
        % disp(obj.riProb)
        % Determine the serial number of the new replication fork

        if ~isempty(obj.repFork)
            newForkNode = obj.oriCstate(2,judgeBox) * 2 + obj.judgeDaughter(judgeBox);
        else, newForkNode = 1;
        end

        % Add new replication fork to obj.repFork
        obj.repFork = [obj.repFork, [zeros(size(newForkNode)); newForkNode]];

        % Sort the replication forks in the order of binary tree
        % number.
        obj.repFork = (sortrows(obj.repFork',2))';

        % Add new replication origins to the system and remove old
        % ones.
        newOriC = zeros(4,size(obj.oriCstate,2));
        newOriC(:,~ismember(obj.oriCstate(2,:),newForkNode)) = obj.oriC(:,~ismember(1:size(obj.oriC,2),judgeBox));

        if obj.model == 1
            newOriC(:, ismember(obj.oriCstate(2,:),newForkNode)) = repmat([obj.hboxPerOri; obj.lboxPerOri; 0; 0],1,sum(ismember(obj.oriCstate(2,:),newForkNode)));
        elseif obj.model == 2
            newOriC(:, ismember(obj.oriCstate(2,:),newForkNode)) = repmat([5; 5; 6; 6],1,sum(ismember(obj.oriCstate(2,:),newForkNode)));
        end

        obj.oriC = newOriC;
    else
        isri = 0;
        releasedDnaA = 0;
    end      

end

function daughterNum = judgeDaughter(obj, daughterLoc)
    %JUDGEDAUGHTER Determines an oriC is located on which daughter node of corresponding replication fork.
    %   Left daughter : Return 0,
    %   Right daughter: Return 1.

    repForkNum = obj.oriCstate(2, daughterLoc);
    daughterNum = zeros(size(daughterLoc));

    for i = 1:length(daughterLoc)
        if daughterLoc(i) > 1 && obj.oriCstate(2,daughterLoc(i)-1) == repForkNum(i)
            % The left sibling exists.
            daughterNum(i) = 1;
        elseif daughterLoc(i) < size(obj.repFork,2) && obj.oriCstate(2,daughterLoc(i)+1) == repForkNum(i)
            % The right sibling exists.
            daughterNum(i) = 0;
        else
            % Neither the left nor the right sibling exists
            if isempty(find(obj.repFork(2,:) == repForkNum(i) * 2, 1))
                daughterNum(i) = 0;
            else, daughterNum(i) = 1;
            end
        end
    end

end

function ister = replisomeAdvancing(obj,dt)
    %REPLISOMEADVANCING updates the location of replication forks and judges if termincation occurs.
    %   Replisome is assumed to be continuously advancing in a 
    %   constant speed without any halt. Once a replication fork 
    %   reaches location of 1, replication termination occurs and 
    %   the current genome divides into two.

    if ~isempty(obj.repFork)
        obj.repFork(1,:) = obj.repFork(1,:) + dt / obj.C;
        ister = ~isempty(find(obj.repFork(1,:) >=1, 1));
    else, ister = false;
    end

end

function nucleoidDivision(obj,n)
    %NUCLEOIDDIVISION Divides one genome into two when termination occurs.
    % Because this function directly alters the properties of
    % genome, to generate two genomes, make deep copy of this
    % object first before going on nucleoid separation.
    % Specifically, n = 1 and n = 2 refers to the daughter binary
    % tree with the mother node number equal to 2 and 3 in original
    % object, respectively.           

    % DEBUG % Show genome info when division
    % disp("GENOME: ")
    % disp("oriC:")
    % disp(obj.oriC)
    % disp("repFork:")
    % disp(obj.repFork)
    % disp("oriCstate:")
    % disp(obj.oriCstate)    

    maxSerialNum = max(obj.repFork(2,:));
    maxLayer = log2(maxSerialNum);
    maxLayer = ceil(maxLayer) + (rem(maxLayer,1) == 0);

    % Map old binary serial number to new one
    div_tab = obj.divBinaryTree(maxLayer,n);

    % if maxium layer is 1 before termination, no replication fork
    % exists after then.
    if isempty(div_tab)
        obj.repFork = []; 
        obj.oriC = obj.oriC(:,n);
        return
    end

    % Generate new oriCs and replication fork binary tree.
    original_oriCstate = obj.oriCstate;
    original_oriC = obj.oriC;
    obj.oriC(:,~ismember(obj.oriCstate(2,:),div_tab(1,:))) = [];
    obj.repFork(:,~ismember(obj.repFork(2,:),div_tab(1,:))) = [];
    obj.repFork(2,:) = div_tab(2,ismember(div_tab(1,:),obj.repFork(2,:)));

    if isempty(obj.oriC)
        obj.oriC = original_oriC(:,original_oriCstate(2,:)==1);
    end

    if size(obj.oriC,2) ~= size(obj.oriCstate,2), error('WRONG!'); end

end

function [nri,unbound,dissoDnaA] = oriCBinding(obj,x,n)
    %ORICBINDING Decides how n DnaAatp molecules is bound to high-affinity or low-affinity DnaA-boxes on different available replication origins.
    %   Parameters:
    %       x: 
    %           x = 1: high-affinity box binding
    %           x = 2: low-affinity box binding
    %       n: number of DnaAatp molecules to bind to boxes
    %       y: dCas9-occupied box (1: high; 2: low)

    % Notably, n may not be a integer. Under that circumstances, an
    % integer is randomly chosen based on n's location between two
    % nearest integers.
    m = n;
    if rem(n,1) ~= 0, n = floor(n) + (rem(n,1) > rand); end

    unbound = m-n;
    nri = 0;
    dissoDnaA = 0;

    if n == 0
        [isri, releasedDnaA] = obj.replicationInitiation();
        nri = nri + isri; 
        dissoDnaA = dissoDnaA + releasedDnaA;
        return; 
    end                

    for i = 1:n
        [~, ~, hboxFree, lboxFree] = obj.boxNum();
        % Determine the available oriCs
        box_count = {hboxFree, lboxFree};
        freeBox = box_count{x};

        if sum(freeBox) == 0
            unbound = n;    
            [isri, releasedDnaA] = obj.replicationInitiation();
            nri = nri + isri; 
            dissoDnaA = dissoDnaA + releasedDnaA;
            return; 
        end

        % The probability of binding to different DnaA-boxes with
        % the same affinity is the same.
        bindingProb = cumsum(sum(freeBox,1)./sum(sum(freeBox,1)));
        boundOriC = find(bindingProb >= rand, 1);

        if obj.model == 1
            freeBox(boundOriC) = freeBox(boundOriC) - 1;
            obj.oriC(x,:) = freeBox;                                   
        elseif obj.model == 2
            bindingOligomerProb = cumsum((freeBox(:,boundOriC))' ./ sum(freeBox(:,boundOriC)));
            boundOligomer = find(bindingOligomerProb >= rand, 1);
            obj.oriC(boundOligomer, boundOriC) = obj.oriC(boundOligomer, boundOriC) - 1;                    
            if obj.oriC(boundOligomer, boundOriC) < 0
                error('Wrong!')
            end

        end               

        % Judge if the binding can lead to replication initiation,
        % which can cause release of DnaAatp and empty DnaA boxes
        % and thus change state of the system.    
        [isri, releasedDnaA] = obj.replicationInitiation();
        nri = nri + isri; 
        dissoDnaA = dissoDnaA + releasedDnaA;
    end

end

function [hboxTot, lboxTot, hboxFree, lboxFree] = boxNum(obj)
    %BOXNUM Returns the total number of high and low-affinity boxes and number of free boxes, which refers to those can still be bound. 
    % The methods of calculation can be different between two 
    % models, due to the difference in defining oriC data.
    % Specifically, those boxes which are located in sequestered
    % oriC or occupied by CRISPR-dCas9 should be excluded from
    % total to obtain the number of free boxes.

    if obj.model == 1
        hboxTot  = obj.oriC(1,:) - obj.oriC(3,:);
        lboxTot  = obj.oriC(2,:) - obj.oriC(4,:);
        hboxFree = obj.oriC(1,:) .* obj.oriCstate(1,:);
        lboxFree = obj.oriC(2,:) .* obj.oriCstate(1,:);
    elseif obj.model == 2
        hboxTot =  [(obj.oriC(1,:)==5);(obj.oriC(2,:)==5)];
        lboxTot =  [sum(obj.oriC([1,3],:)) - 6; sum(obj.oriC([2,4],:)) - 6];
        hboxFree = [(obj.oriC(1,:)==5 .* (obj.oriC(3,:)>1) .*obj.oriCstate(1,:)); (obj.oriC(2,:)==5 .* (obj.oriC(4,:)>1) .* obj.oriCstate(1,:))];
        lboxFree = [(obj.oriC(3,:) + obj.oriC(1,:) - 6) .* obj.oriCstate(1,:); ...
                    (obj.oriC(4,:) + obj.oriC(2,:) - 6) .* obj.oriCstate(1,:)];
    end

end

function bn = numOfSpecBox(obj, k)
    %NUMOFSPECBOX returns the number of kth box in oriC.

    oligoNo = floor(k/5) + (rem(k,5) ~= 0);
    serialNoInOligo = rem(k,5);
    serialNoInOligo = serialNoInOligo + 5 * (serialNoInOligo==0);

    switch serialNoInOligo
        case 1, bn = (obj.oriC(oligoNo,:) == 5) .* (obj.oriC(oligoNo+2,:) > serialNoInOligo) .* obj.oriCstate(1,:);
        otherwise, bn = (obj.oriC(oligoNo,:) > 5-serialNoInOligo) ...
                      .*(obj.oriC(oligoNo+2,:) > serialNoInOligo) ...
                      .* obj.oriCstate(1,:);
    end                          

end  

function datABinding(obj,n)
    %DATABINDING update DnaA binding to datA loci.
    %   n refers to number of this binding event.

    if rem(n,1) ~= 0, n = floor(n) + (rem(n,1) + (n<0) > rand); end

    obj.datA = obj.datA - n; 
    obj.datA = obj.datA * (obj.datA >= 0);
end

function dCasBinding(obj,n,tgt)
    %DCASBINDING updates dCas9 binding to DnaA boxes.
    % This function is also different between two models.
    % For DnaA random filling model, tgt can be 1 or 2, which
    % refer to high and low-affinity boxes, respectively.
    % For Oligomerization started by high-affinity DnaA model, tgt
    % can be any integer from 1 to 10, which refers to R1, R5M, 
    % TAU1, I1, I2, R4, C1, I3, C2 and C3, in order.

    if rem(n,1) ~= 0, n = floor(n) + (rem(n,1) > rand); end            
    if n == 0, return; end

    if obj.model == 1
        isTgtBoxEmpty = zeros(1,size(obj.oriC,2));
        boxPerOri = [obj.hboxPerOri, obj.lboxPerOri];
        tgtBoxPerOri = boxPerOri(tgt);

        for i = 1:size(obj.oriC,2)
            occupiedBox = tgtBoxPerOri - obj.oriC(tgt,i);
            if occupiedBox == 0, isTgtBoxEmpty(i) = 1;
            else
                isTgtBoxEmpty(i) = 1 - nchoosek(tgtBoxPerOri-1, occupiedBox-1)/nchoosek(tgtBoxPerOri, occupiedBox);
            end
        end

        isTgtBoxEmpty = (rand(1,size(obj.oriC,2)) < isTgtBoxEmpty) .* (obj.oriC(tgt+2,:)==0);                
    elseif obj.model == 2
        isTgtBoxEmpty = logical(obj.numOfSpecBox(tgt));
    end

    if n >= sum(isTgtBoxEmpty)
        switch obj.model
            case 1, obj.oriC(tgt+2,isTgtBoxEmpty) = 1;
            case 2, obj.oriC(ceil(tgt/5)+2,isTgtBoxEmpty) = rem(tgt,5) + 5*(rem(tgt,5)==0);
        end
    else                
        nEmptyBoxes = sum(isTgtBoxEmpty);
        randomPerm = randperm(nEmptyBoxes);
        emptyLoc = find(isTgtBoxEmpty);

        switch obj.model
            case 1, obj.oriC(tgt+2,emptyLoc(randomPerm(1:n))) = 1;
            case 2, obj.oriC(ceil(tgt/5)+2,emptyLoc(randomPerm(1:n))) = rem(tgt,5) + 5*(rem(tgt,5)==0);
        end

    end

end

function dCasDissociation(obj,n,tgt)
    %DCASDISSOCIATION updates dCas9 dissociating from DNA.
    if rem(n,1) ~= 0, n = floor(n) + (rem(n,1) > rand); end
    if n == 0, return; end

    switch obj.model
        case 1, nBoundOriC = sum(obj.oriC(tgt+2,:)~=0);
        case 2, nBoundOriC = sum(obj.oriC(ceil(tgt/5)+2,:)<6);
    end

    if n >= nBoundOriC
        switch obj.model
            case 1, obj.oriC(tgt,:) = zeros(1,size(obj.oriC,2));
            case 2, obj.oriC(ceil(tgt/5)+2,:) =6 * ones(1,size(obj.oriC,2));
        end
    else
        randomPerm = randperm(nBoundOriC);

        switch obj.model
            case 1
                boundLoc = find(obj.oriC(tgt+2,:)~=0);
                obj.oriC(tgt,boundLoc(randomPerm(1:n))) = 0;
            case 2
                boundLoc = find(obj.oriC(ceil(tgt/5)+2,:)<6);
                obj.oriC(ceil(tgt/5)+2,boundLoc(randomPerm(1:n))) = 6;
        end

    end
end


function value = get.lastLayerFork(obj)
    %LASTLAYERFORK Querys all replication forks on the last layer of binary tree,
    % which means the corresponding node has 0 or 1 daughter node.

    value = zeros(3,20);
    n = 1;

    for i = 1:size(obj.repFork,2)
        leftDaughter  = find(obj.repFork(2,:)==2*obj.repFork(2,i),1);
        rightDaughter = find(obj.repFork(2,:)==2*obj.repFork(2,i)+1,1);   

        if isempty(leftDaughter) || isempty(rightDaughter)
        if isempty(leftDaughter) && isempty(rightDaughter)                      
            value(:,n) = [obj.repFork(:,i); 2];
        else, value(:,n) = [obj.repFork(:,i); 1];                            
        end
        n = n + 1;
        end

    end

    % Remove blank terms
    value = value(:,value(2,:)~=0);
end

function oriC = get.oriCstate(obj)
    %ORICSTATE Determines the sequestration state of all oriCs according to
    % the locations of last-layer replication forks.

    if isempty(obj.repFork), oriC = [1;0];
    else
    % Pre-define the oriC state information.
    oriC = nan(2, 2*size(obj.lastLayerFork,2)); 
    n = 1;

    for i = 1:size(obj.lastLayerFork,2)
        % OriC is unavailable if sequested
        if obj.lastLayerFork(1,i) <= obj.sqstTime / obj.C
            oriC(:,n) = [0;obj.lastLayerFork(2,i)];
        else, oriC(:,n) = [1;obj.lastLayerFork(2,i)];
        end                

        if obj.lastLayerFork(3,i) == 2
            n = n + 1; 
            oriC(:,n) = oriC(:,n-1); 
        end

        n = n + 1;
    end

    % Remove blank terms
    oriC = oriC(:,~isnan(oriC(1,:)));
    end
end

function totGene = get.totGeneCopyNum(obj)
    %TOTGENECOPYNUM Returns the total gene copy number of the genome

    if isempty(obj.repFork)
        totGene = 1;
    else

    totGene = 0;

    for i = 1:size(obj.repFork,2)
        if obj.repFork(2,i) == 1
            totGene = totGene + 1 - obj.repFork(1,i);
        else
            node = obj.repFork(2,i);
            parentNode = floor(node/2);
            parentLoc  = obj.repFork(1,obj.repFork(2,:)==parentNode);
            totGene = totGene + parentLoc - obj.repFork(1,i);
        end
    end

    totGene = totGene + sum(obj.lastLayerFork(1,:) .* obj.lastLayerFork(3,:));
    end

totGene = totGene * obj.Gnum;
end

function geneNum = geneCopyNum(obj, geneLoc)
    %GENECOPYNUM Returns the copy number of active gene with given gene
    % location on the genome.

    if isempty(obj.repFork)
        geneNum = 1;
    else

    geneNum = sum(obj.lastLayerFork(3,(obj.lastLayerFork(1,:) > geneLoc) ...
                  & (obj.lastLayerFork(1,:) > obj.sqstTime/obj.C)));

    for i = 1:size(obj.repFork,2)
        forkLoc = obj.repFork(1,i);

        if forkLoc < geneLoc
            if obj.repFork(2,i) == 1
                geneNum = geneNum + 1;
            else
                node = obj.repFork(2,i);
                parentNode = floor(node/2);
                parentLoc  = obj.repFork(1,obj.repFork(2,:)==parentNode);
                geneNum = geneNum + (parentLoc > geneLoc & forkLoc < parentLoc);
            end
        end

    end            
    end

end

function plotGenome(obj,x0,y0,theta0,r0)
    %PLOTGENOME Draws quantitative sketch map of the cell genome.

    % Separate part of replication fork
    raiseStage = .25;  % Replication fork seperation fraction
    h0 = .1;   
    dx = 1e-3; % Point interval 

    % Define the line specificity
    linespec1 = {2,'-','k','none',1}; % Genome
    linespec2 = {1,'none','r','.',20}; % Free oriC
    linespec3 = {1,'none','k','.',20}; % Ter
    linespec4 = {1,'none',[.5,.5,.5],'.',20}; % dCas9-bound oriC
    linespec5 = {1,'none','b','.',20}; % Sequestered oriC

    if isempty(obj.repFork)
        obj.line2circlePlot(-1:dx:1,zeros(size(-1:dx:1)),linespec1,x0,y0,theta0,r0), hold on
        dCasOccupyMarker = 6 * (obj.model==2);
        if obj.oriC(3)~=dCasOccupyMarker || obj.oriC(4)~=dCasOccupyMarker
            obj.line2circlePlot(0,0,linespec4,x0,y0,theta0,r0);
        elseif obj.oriCstate(1) == 0
            obj.line2circlePlot(0,0,linespec5,x0,y0,theta0,r0);
        else
            obj.line2circlePlot(0,0,linespec2,x0,y0,theta0,r0);
        end
    else

    for i = 1:size(obj.repFork,2)
        % For each replication fork, three positions need
        % determination: positions of itself loc1, its left 
        % daughter loc21 and right daughter loc 22.
        loc1 = obj.repFork(1,i);
        layer = floor(log2(obj.repFork(2,i))) + 1;    % Layer on the binary tree
        nodeNum = obj.repFork(2,i) - 2^(layer-1) + 1; % which one in the layer
        % Determine the replication fork size
        h = h0 / (2^(layer-1));
        % Middle location of the replication fork in y direction
        hloc = -h0 * 2 + 2 * h + (nodeNum - 1) * h * 4;
        % Determine the position of daughter nodes
        daughterNode = [2 * obj.repFork(2,i), 2 * obj.repFork(2,i) + 1];
        loc21 = obj.repFork(1,obj.repFork(2,:)==daughterNode(1));
        loc22 = obj.repFork(1,obj.repFork(2,:)==daughterNode(2));

        if ~isempty([loc21,loc22])
            minIntv = min([loc1-loc21,loc1-loc22]);
        else, minIntv = 2;
        end

        raiseIntv = min([loc1 * 2 * raiseStage, minIntv]);
        leftRaise  = -loc1 : dx : (-loc1 + raiseIntv);
        rightDescend = (loc1 - raiseIntv) : dx : loc1;

        if i == 1
            obj.line2circlePlot(-1:dx:-loc1,zeros(1,length(-1:dx:-loc1)), linespec1,x0,y0,theta0,r0)
            obj.line2circlePlot(loc1:dx:1,zeros(1,length(loc1:dx:1)), linespec1,x0,y0,theta0,r0)
        end

        obj.line2circlePlot(leftRaise, obj.sigmoid(leftRaise,h,hloc), linespec1,x0,y0,theta0,r0)
        obj.line2circlePlot(leftRaise, obj.sigmoid(leftRaise,-h,hloc), linespec1,x0,y0,theta0,r0)
        obj.line2circlePlot(rightDescend, obj.sigmoid(rightDescend,-h,h+hloc), linespec1,x0,y0,theta0,r0)
        obj.line2circlePlot(rightDescend, obj.sigmoid(rightDescend,h,-h+hloc), linespec1,x0,y0,theta0,r0)

        if ~isempty(loc21)
            obj.line2circlePlot(leftRaise(end):dx:-loc21, (hloc-h)*ones(1,length(leftRaise(end):dx:-loc21)), linespec1,x0,y0,theta0,r0)
            obj.line2circlePlot(loc21:dx:rightDescend(1), (hloc-h)*ones(1,length(loc21:dx:rightDescend(1))), linespec1,x0,y0,theta0,r0)
        else
            obj.line2circlePlot(leftRaise(end):dx:rightDescend(1), (hloc-h)*ones(1,length(leftRaise(end):dx:rightDescend(1))), linespec1,x0,y0,theta0,r0)
            forkNo = obj.repFork(2,i);
            daughterLoc = find(obj.oriCstate(2,:)==forkNo);
            daughterNode = obj.judgeDaughter(daughterLoc);
            corrOriC = obj.oriC(:,daughterLoc(daughterNode==0));
            dCasOccupyMarker = 6 * (obj.model==2);

            if corrOriC(3) ~= dCasOccupyMarker || corrOriC(4) ~= dCasOccupyMarker
                obj.line2circlePlot(0,hloc-h, linespec4 ,x0,y0,theta0,r0)
            elseif obj.oriCstate(1,daughterLoc(daughterNode==0)) == 0
                obj.line2circlePlot(0,hloc-h, linespec5 ,x0,y0,theta0,r0)
            else
                obj.line2circlePlot(0,hloc-h, linespec2 ,x0,y0,theta0,r0)
            end

        end

        if ~isempty(loc22)
            obj.line2circlePlot(leftRaise(end):dx:-loc22, (hloc+h)*ones(1,length(leftRaise(end):dx:-loc22)), linespec1 ,x0,y0,theta0,r0)
            obj.line2circlePlot(loc22:dx:rightDescend(1), (hloc+h)*ones(1,length(loc22:dx:rightDescend(1))), linespec1,x0,y0,theta0,r0)
        else
            obj.line2circlePlot(leftRaise(end):dx:rightDescend(1), (hloc+h)*ones(1,length(leftRaise(end):dx:rightDescend(1))), linespec1,x0,y0,theta0,r0)
            forkNo = obj.repFork(2,i);
            daughterLoc = find(obj.oriCstate(2,:)==forkNo);
            daughterNode = obj.judgeDaughter(daughterLoc);
            corrOriC = obj.oriC(:,daughterLoc(daughterNode==1));
            dCasOccupyMarker = 6 * (obj.model==2);

            if corrOriC(3) ~= dCasOccupyMarker || corrOriC(4) ~= dCasOccupyMarker
                obj.line2circlePlot(0,hloc+h, linespec4 ,x0,y0,theta0,r0)
            elseif obj.oriCstate(1,daughterLoc(daughterNode==1)) == 0
                obj.line2circlePlot(0,hloc+h, linespec5 ,x0,y0,theta0,r0)
            else
                obj.line2circlePlot(0,hloc+h, linespec2 ,x0,y0,theta0,r0)
            end

        end      

    end

    end

    obj.line2circlePlot(1,0,linespec3,x0,y0,theta0,r0)
    hold off, axis off
end    

end

methods(Static)
function div_table = divBinaryTree(layer,n)
    %DIVBINARYTREE Divides a binary tree into two.
    %   DIVBINARYTREE(layer, n) Returns an array whose first line 
    %   is the serial number of nodes in old tree and second line 
    %   is the serial number of new tree.
    %
    %   Parameters:
    %       layer: number of layers in old tree
    %       n:
    %           n = 1: Left daughter tree
    %           n = 2: Right daughter tree

    N = 2^layer - 1; div_table = zeros(2,(N-1)/2); allNum = 2:N;

    if n == 1
        div_table(1,:) = allNum(allNum-2.^(floor(log2(allNum))) - 2.^(floor(log2(allNum))-1) < 0);
    else
        div_table(1,:) = allNum(allNum-2.^(floor(log2(allNum))) - 2.^(floor(log2(allNum))-1) >= 0);
    end

    div_table(2,:) = div_table(1,:) - 2.^(floor(log2(div_table(1,:)))-1) * n;            
end  

function y = sigmoid(x,h,b)
    %SIGMOID Defines the shape of DNA separation at the replication fork.

    xRg = 5; 
    xmax = max(x); 
    xmin = min(x);
    y = b + h ./ (1 + exp(- (x-xmax) / (xmax - xmin) * 2 * xRg - xRg));    
end

function line2circlePlot(x,y,linespec,x0,y0,theta0,r0)
    %LINE2CIRCLEPLOT Maps linear genome to circular genome.
    %   x0 is for adjusting angle of genome:
    %       x0 is even: oriC at left
    %       x0 is odd: oriC at right
    %   y0 is for adjusting the shape of the genome

    r = y + y0 + 1; 
    theta = 0 + (x + x0 + 1) * pi; 
    [theta_new, r_new] = cart2pol(r .* cos(theta) + r0 * cos(theta0), r .* sin(theta) + r0 * sin(theta0));
    p = polarplot(theta_new,r_new);
    [p.LineWidth, p.LineStyle, p.Color, p.Marker, p.MarkerSize] = linespec{:}; 
    hold on
end        

end
        
end
