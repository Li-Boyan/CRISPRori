function IntegratedSimulation(varargin)
%INTEGRATEDSIMULATION Simulates cell growth, gene expression and cell cycle processing in Escherichia coli, with customized parameters and conditions.
%     IntegratedSimulation() runs the simulation with default conditions, where
%     concentration of dCas9 is 0, total simulation time is 1,440 minutes, 
%     default model is DDOM and no exogenous plasmids is introduced. All other 
%     parameters are using default values. The output will be a structure array
%     which contains:
%             1) All parameters and conditions used;
%             2) Time-course lineage information;
%             3) Lineage information at cell division;
%             4) Lineage information at DNA replication initiation;
%     The last three are all using tables as datatype. The output is saved to a
%     .mat file named according to the conditions used under an automatically 
%     created folder called SimulationResults.
% 
%     IntegratedSimulation(param1, value1, param2, value2, ..., paramn, valuen)
%     enables simulation with customized parameters and conditions. NOTE that
%     the parameter space is not thoroughly tested and hence some unappropriate
%     customized parameters can lead to error.
% 
%     Basic cell parameters can be specified as follows:
% 
%             nDnaA:     Hill coefficient of autorepression function in DnaA transcription. Default is 0.1.
%             nu_Rb:     Maximum transcription frequency of ribosomal proteins. Default is 30 /min.
%             nu_R:      Maximum transcription frequency of RNAp proteins. Default is 30 /min.
%             nu_DnaA:   Maximum transcription frequency of DnaA. Default is 2 /min.
%             nu_FtsZ:   Maximum transcription frequency of FtsZ. Default is 5 /min.
%             nu_total:  Maximum transcription frequency of general proteins. Default is 1 /min.
%             Kt:        Number of RNAp half saturating a promoter. Default is 0.4.
%             alphaf:    Fraction of free RNAp with respect to total. Default is 0.22.
%             taum_DnaA: Degradation time of DnaA mRNA. Default is 13.6 minutes.
%             taum_Rb:   Degradation time of ribosomal mRNA. Default is 17.5 minutes.
%             taum_total:Degradation time of general mRNA. Default is 9.8 minutes.
%             taum_FtsZ: Degradation time of FtsZ mRNA. Default is 16.3 minutes.
%             taum_R:    Degradation time of RNAP mRNA. Default is 9.8 minutes.
%             ve:        Speed of the ribosome. Default is 2.88 Kbp/min.
%             sigma_DnaA:Translation rate of DnaA. Default is 3 /min.
%             sigma_total, sigma_FtsZ, sigma_R, sigma_Rb: 
%                        Translation rate of general mRNA, FtsZ, RNAp and ribosomes,
%                        respectively. Defaults are all 10 /min.
%             Ke:        Number of ribosomes half saturating a mRNA. Default is 7.2.
%             a_inact:   Fraction of inactive ribosomes. Default is 0.0020.
%             ratp2adp:  Ratio of intracellular ATP to ADP. Default is 3.2.
%             katp:      Equilibrium constant of ATP + DnaA <=> DnaA-ATP. Default is 0.05.
%             kadp:      Equilibrium constant of ATP + DnaA <=> DnaA-ATP. Default is 0.0125.
%             taup_FtsZ: Degradation time of FtsZ. Default is 105 minutes.
%             C:         C-period. Default is 40 minutes.
%             sqstt:     Sequestration period. Default is 10 minutes.
%             khbox:     Binding rate coefficient of DnaA-ATP to high-affinity boxes. Default is 0.06.
%             klbox:     Binding rate coefficient of DnaA-ATP to low-affinity boxes. Default is 0.006.
%             kdatA:     Binding rate coefficient of DnaA-ATP to boxes in datA locus. Default is 0.006.
%             kRIDA:     Rate cofficient of conversion of DnaA-ATP to DnaA-ADP through
%                        replicative inactivation of DnaA-ATP. Default is 0.2.
%             kon:       Association rate coefficient of dCas9 to DNA. Default is 0.005.
%             koff:      Dissociation rate coefficient of dCas9 from DNA. Default is 0.005.
%             nu_GFP:    Transcription frequency of plasmid-encoded GFP. Default is 10 /min.
%             taum_GFP:  Degradation time of GFP mRNA. Default is 9.8 minutes.
%             sigma_GFP: Translation rate of GFP. Default is 10 /min.
% 
%     The initial conditions can be specified as follows:
% 
%             Rb_init:   Number of ribosomes. Default is 45,000 /um^3.
%             R_init:    Number of RNAp. Default is 40,000 /um^3.
%             DnaA_init: Number of DnaA. Default is 1,000 /um^3.
%             FtsZ_init: Number of FtsZ. Default is 10,000 /um^3.
%             total_proteins_init:
%                        Number of total proteins. Default is 7,500,000 /um^3.
%             total_mRNA_init:
%                        Number of total mRNA. Default is 8,400 /um^3.
%             initRepFork:
%                        Initially existing replication forks. Default is an empty array.
% 
%     The simulation conditions can be specified as follows:
% 
%             c_dCas9:   Concentration of dCas9. Default is 0.
%             dCas9_target:
%                        Target box of dCas9. Default is R1.
%             ttol:      Total time of simulation. Default is 1,440 minutes.
%             plasmidCopyNum:
%                        Concentration of plasmids. Default is 0 /um^3. 
% 
%     Simulation parameters can be specified as follows:
% 
%             nRcd:      Total recorded timepoints. Default is 10,000.
%             dt:        Update interval. Default is 0.01 min.
%             gif:       Boolean variable to decide if a gif is generated or not. Default is false. 
%             outputPath:  Path for output. Default is "SimulationResults".
%
%
load DefaultParams.mat % Default parameters and conditions.
% Customize the parameters
for i = 1:2:numel(varargin)
    try eval(strcat(varargin{i},'=',varargin{i+1},';'))
    catch
        eval(strcat(convertCharsToStrings(varargin{i}), "='", ...
                    convertCharsToStrings(varargin{i+1}), "'",';'));
    end
end
Parameters = {nu_Rb, nu_R, nu_DnaA, nu_FtsZ, nu_total, Kt, alphaf, KDnaA, ...
              nDnaA, taum_Rb, taum_R, taum_DnaA, taum_FtsZ, taum_total,   ...
              sigma_Rb, sigma_R, sigma_DnaA, sigma_FtsZ, sigma_total, ve, ...
              a_inact, Ke, katp, kadp, ratp2adp, taup_FtsZ, kRIDA, ...
              kdatA, khbox, klbox, kon, koff,nu_GFP,taum_GFP,sigma_GFP};   
genomeInitCondition = {sqstt, C, initRepFork};
proteins_0 = [Rb_init, R_init, DnaA_init, DnaA_init*12/13.6, DnaA_init*6/13.6, ...
              DnaA_init/13.6, FtsZ_init, total_proteins_init, 0, 0];
mRNA_0 = [total_mRNAs_init * [Rb_init, R_init, DnaA_init, FtsZ_init] / total_proteins_init, total_mRNAs_init];
substrates = {'Ribosome_mRNA', 'RNAP_mRNA', 'DnaA_mRNA', 'FtsZ_mRNA', 'Total_mRNA', ...
              'Ribosome', 'RNAP', 'DnaA', 'DnaAatp', 'DnaAatp_free', 'DnaAadp', ...
              'FtsZ', 'Total_proteins', 'GFP_mRNA','GFP'};
initCondition = containers.Map(substrates, [mRNA_0, proteins_0]);

% Simulation
rng('shuffle');
Ecoli = EscherichiaColi(initCondition, genomeInitCondition, nRcd, Parameters, FtsZ_init*ftsZthres, c_dCas9, dCas9_target, plasmidCopyNum, model);
Ecoli.generateLineage(dt, ttol, gif, fullfile(outputPath, 'test'), exittime);

% Results output
if ~exist(outputPath,'dir'), mkdir(outputPath); end
randNum = ceil(26*rand(1,4));
alphabet = 'A':'Z';
newFileName = ['lineage-', datestr(now,'mmddHHMMSSFFF'), alphabet(randNum), '.mat'];

% Save results
switch outputType
case 'data'
varName = Ecoli.substanceData.Properties.VariableNames;
lineage.timeCourse = array2table(Ecoli.substanceData{Ecoli.substanceData.label==0,:}, 'VariableNames', varName);
lineage.ri = array2table(Ecoli.substanceData{Ecoli.substanceData.label==1,:}, 'VariableNames', varName);
lineage.di = array2table(Ecoli.substanceData{Ecoli.substanceData.label==2,:}, 'VariableNames', varName);              
lineage.params.parameters = Parameters;
lineage.params.initCondition = initCondition;
lineage.params.genomeInitCondition = genomeInitCondition;
lineage.params.condition = {c_dCas9, dCas9_target, plasmidCopyNum, model};
save(fullfile(outputPath, newFileName), 'lineage')
case 'object'
save(fullfile(outputPath, newFileName), 'Ecoli')
end