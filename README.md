This code calibrates and evaluates nested Species Distribution Models (SDMs).

The scripts folder contains one sub-folder per step of modelling, to be run in order. We describe here their content.
- 1_data_preparation : Download and clean species and environment data.
- 2_M1 : Calibrate the coarse-grain model, at 10 km spatial resolution and pluriannual temporal resolution.
- 3_M2: Calibrate the fine-grain model, at 1 km spatial resolution and annual temporal resolution.
- 4_combine : Combine these two models, with the geometric mean or the permanence of ratios.
- 5_stats_and_maps : Generate the figures and tables for the report.
