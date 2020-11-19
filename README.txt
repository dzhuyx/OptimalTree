This folder contains codes used for simulations studies and data analysis for the manuscript "Obtaining Optimal Cutoff Values for Tree Classifiers Using Multiple Biomarkers."

For simulation studies under correctly specified models, relevant files are:
1. "simulation_correct.R", which is the code for running one replication of one scenario of simulation under a correctly specified model;
2. "simulation_correct_trueROC.R", which calculates true optimal ROC curves for each scenario through Monte-Carlo;
3. "simulation_correct_runbatch.R", which is provided as an illustration of the whole simulation study and can be used to run the whole simulation study under correctly specified models, but a parallel computing equivalence on a computing cluster is probably preferred;
4. "simulation_correct_summary.R", which summarizes simulation results obtained by running "simulation_correct.R" and "simulation_correct_runbatch.R", and formats them into tables shown in the manuscript;
5. "simulation_wboot.R", which is for running one replication of simulation in Section 7.1 with bootstrapping for variance estimation;
6. "simulation_wboot_runbatch.R", which is provided as an illustration of the whole simulation study in Section 7.1 and can be used to run the whole simulation study, but a parallel computing equivalence on a computing cluster is probably preferred;
7. "simulation_wboot_summary.R", which summarizes simulation results obtained from "simulation_wboot.R" and "simulation_wboot_runbatch.R".

For simulation studies under misspecified models, relevant files are:
1. "simulation_misspecified.R", which is the code for running one replication of one scenario of simulation under a correctly specified model;
2. "simulation_misspecified_runbatch.R", which is provided as an illustration of the whole simulation study and can be used to run the whole simulation study under correctly specified models, but a parallel computing equivalence on a computing cluster is probably preferred;
3. "simulation_misspecified_summary.R", which summarizes simulation results obtained by running "simulation_misspecified.R" and "simulation_misspecified_runbatch.R", and formats them into tables shown in the manuscript.

File "simulation_correct_trueROC.R" should be run before files "simulation_correct_summary.R" and "simulation_misspecified_summary.R" can be run.

For real data analysis, the relevant file is "DREAM_analysis.R".

The remaining two files, "tau-a.cpp" and "fun-tau_n.cpp", are implementations of the proposed estimation method in the manuscript.

