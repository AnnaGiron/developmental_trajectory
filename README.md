# developmental_trajectory

## Datasets:    
- `data/behavioralData.csv`: behavioral data of all 281 participants   
- `data/modelFit.csv`: parameter estimates of all models for all participants   
- `data/modelFit_OriginalID`: parameter estimates of all models for all participants with original IDs from the Schulz et al and Meder et al (needed for recovery analyses)
- `data/paramsYoungestAgegroup.csv`: cross-validated parameter estimates of the youngest age group (used as starting points for the optimization algorithms)

## Scripts:   
- `dataProcessing.R`: import and pre-process behavioral data and parameter estimates
-  `behavior_tests.R` and `behavior_plots.R`: analyze and plot participant's behavior in the multi-armed bandit task
-  `crossvalidation.R`: optimize parameters
-  `learningCurves.R`: simulate learning curves
-  `modelResults_tests.R` and `modelResults_plots.R`: analyze and plot model results
-  `reliabilityChecks_tests.R` and `reliabilityChecks_plots.R`: compare participant's performance in the multi-armed bandit task and model results across experiments
-  `simulateModels.R` and `simulateModels_plots.R`: simulate the GP-UCB model with different parameter combinations and plot expected rewards
-  `hillClimbingAlgorithm.R`, `hillClimbingAlgorithm_tests.R` and `hillClimbingAlgorithm_plots.R`: run optimization algorithms in the parameter space calculated in `simulateModels.R`
