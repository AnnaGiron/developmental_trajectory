# developmental_trajectory
## Preprint:
https://charleywu.github.io/downloads/giron2022developmental.pdf

## Datasets:    
- `data/behavioralData.csv`: behavioral data of all 281 participants   
- `data/modelFit.csv`: parameter estimates of all models for all participants   
- `data/modelFit_OriginalID`: parameter estimates of all models for all participants with original IDs from the Schulz et al and Meder et al (needed for recovery analyses)
- `data/paramsYoungestAgegroup.csv`: cross-validated parameter estimates of the youngest age group (used as starting points for the optimization algorithms)

## Scripts:   
- `dataProcessing.R`: import and pre-process behavioral data and parameter estimates (added for reference, since we already include the generated outputs instead of the inputs)
- `statisticalTests.R` contains wrapper functions for all statistical tests used for the analyses
-  `behavior_tests.R` and `behavior_plots.R`: analyze and plot participant's behavior in the multi-armed bandit task    
   (generates Figure 2)
-  `crossvalidation.R`: optimize parameters of the GP-UCB model and the lesioned models. All models being fit here are defined in `models.R`.
-  `learningCurves.R`: simulate learning curves
-  `PXP.ipynb` compute and save the protected exceedance probability (*pxp*) for all models and age groups. Functions to compute the *pxp* are defined in `bms.py` and files containing the negative log likelihoods that are imported in the notebook are created in `modelResults_tests.R`.
-  `modelResults_tests.R` and `modelResults_plots.R`: analyze and plot model results.  Therefore, simulated learning curves from `learningCurves.R` and *pxp*s from `PXP.ipynb` are imported.    
   (generates Figure 3 and S2)
-  `reliabilityChecks_tests.R` and `reliabilityChecks_plots.R`: compare participant's performance in the multi-armed bandit task and model results across experiments    
   (generates Figure S1)
-  `simulateModels.R` and `simulateModels_plots.R`: simulate the GP-UCB model with different parameter combinations and plot expected rewards    
   (Generates Figure S5)
-  `hillClimbingAlgorithm.R`, `hillClimbingAlgorithm_tests.R` and `hillClimbingAlgorithm_plots.R`: run optimization algorithms in the parameter space calculated in `simulateModels.R`   
   (generates Figure 4 and S6)
   
## Model and Parameter Recovery
- `Recovery/Model_Recovery.Rmd` import, analyze and plot model recovery results generated with `Recovery/Model_Recovery_Cluster.R`, `Recovery/Model_Recovery_Cluster_Meder.R` and `Recovery/Model_Recovery_Cluster_Schulz.R` and saved in `Recovery/modelRecovery/`. Models used for recovery are defined in `fit_parallel_cluster.R`.    
  (generates Figure S3 )
- `Recovery/Parameter_Recovery_check.R` import, analyze and plot parameter recovery results generated with `Recovery/Parameter_Recovery_Cluster.R` and saved in `Recovery/parameterRecovery/`. Models used for recovery are defined in `Recovery/fit_parallel_cluster.R`.    
  (generates Figure S4)
