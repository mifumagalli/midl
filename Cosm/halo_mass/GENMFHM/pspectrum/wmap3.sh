#!/bin/bash
TOPGRID=0.03    # Top grid spacing (Mpc)
TOPSIZE=1000    # Top grid size (N where simulation is NxNxN)
ISEED=1915574050 # Random seed

# Parameters taken from WMAP3 data and can be found at:
# http://lambda.gsfc.nasa.gov/product/map/dr2/params/lcdm_wmap.cfm
# Parameters from the PDG can be found at:
# http://pdg.lbl.gov/2006/reviews/hubblerpp.pdf
OMEGAB=0.041    # Baryon fraction (WMAP3)
OMEGAC=0.196    # CDM fraction (WMAP3)
OMEGAV=0.763    # Dark Energy fraction (WMAP3)
OMEGAN=0.000    # None
H0=73.5         # Hubble parameter (WMAP3)
SPECTRAL=0.951  # spectral index n (PDG)
SIGMA8=-0.742   # Sigma 8 - normalization (WMAP3)
