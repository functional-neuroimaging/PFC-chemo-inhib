codes in this folder are written for the analyses of hM4Di multi-electrode data for the Rocchi et al paper.
requires the preprocessed .mat files. 

How to use:
There are 3 routines for these analyses:
1- the main.m routine calculates the spectrum and coherency
2- the Phase_relationships.m routine calculates the PLV
3- the Gamma_envelope.m routine calculates the coherency of gamma envelope based on Nir paper

The codes need to be performed "only once", AT THE BEGINNING OF EACH ROUTINE THERE ARE MANY PARAMETERS TO BE DEFINED, PLEASE READ THEM CAREFULLY.
-run the Main routine once. it will save the results in the current directory in a folder called results\'(name_of_each_animal)'.
-run the Phase_relationships routine once. it will save the results in the current directory in a folder called Results_angle\.
-run the Gamma_envelop routine once. it will save the results in the current directory in a folder called Results_gamma_envelope\.


After finishing, you can use the following codes for the statistics and results visualization:
AT THE BEGINNING OF EACH CODE THERE ARE MANY PARAMETERS TO BE DEFINED, PLEASE READ THEM CAREFULLY.
-stats_main.m: visualizes the statistics for the coherency and spectrum, requires results from the Main.m
-Phase_results.m: visualizes the statistics for the PLV, requires results from the Phase_relationships.m
-gamma_envelope_results.m visualizes the statistics for the gamma_envelope coherency, requires results from the Gamma_envelope.m
-visualize_abs_coherency.m visualizes the raw coherency, requires the results from the Main.m
-visualize_abs_gamma_envelope.m visualizes the raw coherency of gamma_envelope, requires the results from the Gamma_envelope.m


in addition, there is a file named spike_synchrony_power_relation.m, it is just click and run. it is a small simulation to demonstrate how increased spike coupling can lead to increase of power
--------------------------------------------------------------
Shahryar Noei