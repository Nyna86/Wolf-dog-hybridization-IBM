# Simulating the efficacy of wolf–dog hybridization management with individual-based modelings

This repository contains the R code to run the individual-based model (IBM) presented in Santostasi et al. 2024 () which simulated wolf life cycle including individual behavior, pack dynamics and hybridization dynamics.
The model is described in the following 2 papers: 

- Santostasi, N. L., Bauduin, S., Grente, O., Gimenez, O., & Ciucci, P. (2024). Simulating the efficacy of wolf–dog hybridization management with individual-based modeling. Conservation Biology.

- Bauduin, S., Grente, O., Santostasi, N. L., Ciucci, P., Duchamp, C., & Gimenez, O. (2020). An individual-based model to explore the impacts of lesser-known social dynamics on wolf populations. Ecological Modelling

The file initParam.R builds the initial wolf population needed to run the IBM simulations and set up the model parameters. 
The file submodels.R includes all sub-models used in the wolf IBM and detailed in the Methods section of the paper. 
The file run.R runs the wolf IBM and extracts some results: 
1) it calls the submodel.R file to load all sub-models,
2) it calls the initParam.R file to create the initial population and load the model parameters,
3) a loop organizes the different sub-models, runs the simulation and records the outputs and
4) some results are extracted from saved model outputs and figures are produced. Detailed comments are included in each file.

Nina L. Santostasi October 2021
