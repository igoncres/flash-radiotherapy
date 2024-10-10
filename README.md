# flash-radiotherapy
This folder contains the files to reproduce the results presented in the article:
I. González-Crespo, F. Gómez, Ó. López-Pouso and J. Pardo-Montero 2024 "An in-silico study of conventional and FLASH radiotherapy iso-effectiveness: potential impact of radiolytic oxygen depletion on tumor growth curves and tumor control probability", Physics in Medicine and Biology.

----------------------------------------------------------------------------------

Folders:
    Data_figures: contains the .mat files to reproduce the figures of the article
    MATLAB_functions: contains the code to reproduce the results of the article
    
    
----------------------------------------------------------------------------------


    Data_figures

The name of the file indicates the corresponding figure in the article.

** TumorSample **
This folder contains the parameters to generate the 100 simulated tumors with heterogeneous oxygenation in the Omega sample.
    
    
----------------------------------------------------------------------------------
    
    MATLAB_functions

** mainFLASH **
Generates a simulated tumor with a heterogeneous distribution of capillaries and calculates the initial oxygen distribution as well as its evolution during FLASH-RT.

** createGeometry **
Generates the geometry of a simulated tumor: a squared domain with a capillary distribution represented by inner circles.
(Edited from Rodríguez-Barbeito et al.)

** decsg_fix **
This function fixes a bug in the MATLAB function decsg.
(Taken from Rodríguez-Barbeito et al.)

** pO2DiffusionSolver **
Defines the equation for the oxygenation problem in CONV-RT.
(Edited from Rodríguez-Barbeito et al.)

** flashDiffusionSolver **
Defines the equation for the oxygenation problem in FLASH-RT.

** calculateSF **
Evaluate the surviving fraction (SF) in CONV-RT using the LQ model with OERs.

** calculateFlashSF **
Evaluate the surviving fraction (SF) in FLASH-RT using the LQ model with OERs and the methodology proposed by Taylor et al.

** getSF **
Uses the above two functions to obtain the SF with CONV-RT and FLASH-RT.

** ponderation **
Function to evaluate the mean value of any variable from its values in the nodes of the mesh.
(Edited from Rodríguez-Barbeito et al.)

** getVolume **
Code to solve the model of tumor evolution.

** annealing_optimization **
Implementation of the simulated annealing algorithm to fit the volume curves to experimental data.

** calculateTCP_hom **
Calculates tumor control probability (TCP) in homogeneously oxygenated tumors.

** calculateTCP_het **
Calculates the tumor control probability (TCP) in heterogeneously oxygenated tumors.

** getD50_D90 **
Obtains the D50 and D90 from TCP curves.

----------------------------------------------------------------------------------
