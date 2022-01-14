codes in this folder are written for the analyses of hM4Di single electrode data for the Rocchi et al paper.
requires the preprocessed .mat files. 

How to use:
run the Main.m routine for each animal separately, at the beginning of the routine, it will ask you to go to the animal folder and select the traces you want for analyses. it will save the results in the current directory in a folder called results\'(depends on the input of the code, please check)'
For the spike-LFP relationships, run the spiking_Properties.m routine for each animal separately, it will ask you to go to the animal folder and select the traces you want for analyses. it will save the results in the current directory in a folder called results\MUA_statistics
After finishing, with both routines run the results_check.m, it will ask you for the folder that results from the Main routine are saved and asks you to select them, then it visualizes the results with statistical assessment
then it will ask you for the folder that results from the spiking Properties routine are saved and asks you to select them, then it visualizes the results with statistical assessment

--------------------------------------------------------------
Shahryar Noei