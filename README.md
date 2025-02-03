# MFM-Inv-Ga
This repository provides R codes for MFM-Inv-Ga and MFM-Ga. These methods are dealt with in the following paper.

Iwashige. F., Hashimoto. S. (2025). Bayesian mixture modeling using a mixture of finite mixtures with normalized inverse Gaussian weights. 	arXiv:2501.18854

This repository has the following function files. 

* ```IG_functions.R```: for cluster analysis and density estimation with MFM-Inv-Ga  
* ```Ga_functions.R```: for cluster analysis and density estimation with MFM-Ga  
* ```IG_functions_network.R```: for community detection with MFM-Inv-Ga  
* ```Ga_functions_network.R```: for community detection with MFM-Ga
* ```density_output.R```: function file for MCMC summary of density estimation 

This repositoru has the following files for one shot example.

* ```cluster_analysis.R```: cluster analysis with MFM-Inv-Ga and MFM-Ga 
* ```density_estimate_galaxy```: density estimation for galaxy data with MFM-Inv-Ga and MFM-Ga
* ```community_detection.R```: community detection with MFM-Inv-Ga and MFM-Ga

Note that call the appropriate function file according to your purpose. For example, to perform cluster analysis, call ```IG_functions.R``` and ```Ga_functions.R``` in ```cluster_analysis.R```.
