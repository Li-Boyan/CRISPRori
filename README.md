# CRISPRori
Numerical simulation programs for modeling DNA replication, gene expression, and cell division of *Escherichia coli* jointly. In our paper, we report a novel discovery of DNA's role in governing cell growth homeostasis and resource allocation.

For the details of this model or experimental validations, Please refer to our paper:
**Scaling between DNA and cell size governs bacterial growth homeostasis and resource allocation**

This code repository includes the main programs for the numerical simulation and data analysis that has been test on MATLAB R2020a. To start the simulation, download the whole repository to your local directory, and run the following command in the terminal:

    matlab -nosplash -nodesktop -r IntegratedSimulation
    
Primary paramters include concentration of dCas9, target of dCas9, total time, plasmid copy number. For example, if you want to perform the simulation with dCas9 concentration of 50 targeting R4 box, total time of 600 minutes, GFP-encoded plasmid number of 10, and returning lineage GIF, the command should be:

    matlab -nosplash -nodesktop -r IntegratedSimulation c_dCas9 50 dCas9_target R4 ttol 600 plasmidCopyNum 10 gif True
    
There are a number of other parameters that can be altered. For details, please run *help IntegratedSimulation* in MATLAB. The simulation results will be in the default directory "../results" in .mat format. We also provided a batch processing program to automatically deal with those result files. 
