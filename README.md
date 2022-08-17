# PosterCNMAC_2022
 Numerical Experiments for CNMAC 2022
 # Folders and Files

* `ArtificialAVTA.m`: Routine used to generate the results of poster

- Functions: auxiliary functions used by ArtificialAVTA.m
    * `AVTA_?.m`: Solve the irredundance problem for each algorithm: Triangle Algorithm (TA) `ta_anti_warm.m`, Spherical Triangle (SPH)`Spherical_TA11.m`, Frank Wolfe (FW) `GreedyTriangleAlgorithmAVTA.m` and Away Step Frank Wolfe (ASFW) `AwayStepFrankWolfeAlgorithmAVTA.m`.
    * `Random_pts.m`: Generates A randomly according to a uniform distribution on the unit ball of Rm
    * `Random_cvs.m`: Generates randomly redundances points for A
    

- Figures: Figures included in Poster

To reproduce the experiments as they are in the poster, add the folder 'Functions' to MATLAB path and use the routine ArtigicialAVTA.m. 


