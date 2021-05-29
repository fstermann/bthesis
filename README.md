# Bachelors Thesis in Statistics @LMU Munich
## Einfluss von Abhängigkeitsstrukturen auf Clustering Performance nach Dimensionsreduktion
\
\
\
This repository contains code used in my bachelors thesis, as well as plots and results of the analysis.




    ├── data                        <- The data folder is generated with the function create_folders()
    │   │                              and relevant for the data generation
    │   │
    │   ├── datasets                <- Datasets with distributionparameters, 
    │   │                              correlation matrices, covariance matrices, 
    │   │                              computation time and seed
    │   ├── correlation_matrices    <- Correlation matrices
    │   └── correlation_bounds      <- Upper and lower correlation bounds
    │       └── quantiles           <- Calculation of the quantile function
    │                                  that is used for the upper and lower 
    │                                  bounds
    │
    ├── ari                         <- Results from the ARI analysis
    ├── graphics                    <- Produced plots
    │		
    ├── execution_file.R            <- Main code that calls functions from 
    │                                  the files below
    │           
    ├── data_generation.R           <- Functions for the data generation
    ├── dim_red.R                   <- Functions for dimensionality reduction
    ├── clustering.R                <- Functions for cluster analysis
    ├── correlations.R              <- Functions to visualize the correlations
    ├── auxillary.R                 <- Auxillary functions
    │
    └── session_info.txt            <- Info about the R version and loaded packages
