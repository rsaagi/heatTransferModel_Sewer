# Heat transfer model for sewer pipe - Malmö case-study
This repo contains the code used to model heat transfer for sewer pipes. 
The model is used for a case-study in Malmö.
The model here generates the results mentioned in the Environmental Modelling & Software paper submitted (under revision)

This model is mainly meant for reviewers at EMS to access the code, if required. Any other readers are suggested to wait until the paper is published to download the final version of the code.

<strong>Steps to run the model</strong>
1. mexall - This will mex the c-files that will be used in the simulation
2. run_sewermodel_malmo_mechanistic - This will run the mechanistic model
3. run_sewermodel_malmo_conceptual - This will run the conceptual model
4. malmo_runall - This will run both the mechanistic and conceptual model
