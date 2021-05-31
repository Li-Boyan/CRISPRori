# CRISPRori
Numerical simulation programs for modeling DNA replication, gene expression, and cell division of bacteria. 

In our paper, we report a novel discovery of DNA's role in governing cell growth homeostasis and resource allocation. For the details of this model or experimental validations, Please refer to our paper:
**Scaling between DNA and cell size governs bacterial growth homeostasis and resource allocation**

This code repository includes the main programs for the numerical simulation and data analysis that has been tested on MATLAB R2020a. To start the simulation, download the whole repository to your local directory, and run the following command in the terminal:

    matlab -nosplash -nodesktop -r "IntegratedSimulation"
    
Primary paramters include concentration of dCas9, target of dCas9, total time, and plasmid copy number. For example, if you want to perform the simulation with dCas9 concentration of 50 targeting R4 box, total time of 600 minutes, GFP-encoded plasmid number of 10, and returning lineage GIF with 100 frames, the command should be:

    matlab -nosplash -nodesktop -r "IntegratedSimulation c_dCas9 50 dCas9_target R4 ttol 600 plasmidCopyNum 10 gif True nRcd 100"
    
There are a number of other parameters that can be altered. For details, please run **help IntegratedSimulation** in MATLAB. The simulation results will be in the default directory "../results" in the \*.mat format. We also provide a batch processing program to automatically deal with those result files. To run the processing program, switch to the *analysis* directory, and run the following command directly:

    matlab -nosplash -nodesktop -r "lineageAnalysis"
    
If you want to change the default parameters, please run **help lineageAnalysis** in MATLAB. The program will return following results in csv files:

    (1) Metadata that specifies the parameters used in each group;
    (2) Instantaneous growth rate;
    (3) Division cycle and cell volume;
    (4) Division adder;
    (5) Nucleoid-to-cytoplasm ratio;
    (6) GFP proteomic fraction;
    (7) Initiation mass;
    (8) Cell volume at certain time point;
    (9) Instantaneous GFP production rate.
    
The results directory gives some sample result files and analysis/results are the analysis results of these simulations.
