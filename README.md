# MFM-IGau
This repository provides R codes for MFM-IGau, DMFM-IGau, MFM-Ga and DMFM-IGau. These methods are dealt with in the following paper.

Iwashige. F., Hashimoto. S. (2025). Bayesian mixture modeling using a mixture of finite mixtures with normalized inverse Gaussian weights. 	arXiv:2501.18854

This repository has the following function files. Please refer to the comments for details on the model, inputs and outputs.

* ```IG_functions_MFM.R``` and ```IG_functions_DMFM.R```: for cluster analysis and density estimation via MFM-IGau and DMFM-IGau
* ```IG_functions_MFM_network.R``` and  ```IG_functions_DMFM_network.R```: for community detection via MFM-IGau and DMFM-IGau
* ```IG_functions_general.R```: This file is basic fail for all files involved MFM-IGau and DMFM-IGau. Please load this file before above files are loaded.
* ```summary.R```:  This file is used to summarize MCMC output.

The same applies to the function files for MFM-Ga and DMFM-Ga (```Ga_functions_MFM.R```, ```Ga_functions_DMFM.R```, ```Ga_functions_MFM_network.R```,  ```Ga_functions_DMFM_network.R``` and ```Ga_functions_general.R```)

This repository has the following files for one shot example.

* ```cluster_analysis_git.R```: cluster analysis via MFM-IGau, DMFM-IGau, MFM-Ga and DMFM-Ga.
* ```commu_git.R```: community detection via MFM-IGau, DMFM-IGau, MFM-Ga and DMFM-Ga.
* ```density_estimate_git.R```:density estimation for the galaxy data via MFM-IGau and MFM-Ga.
