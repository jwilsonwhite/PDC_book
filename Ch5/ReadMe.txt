ReadMe for Chapter 5

example_dndm_figs.m creates Fig. 5.1

urchin_sizedist_plot.m creates Fig. 5.2, using data in Urchin_data_tegner.csv. Also calls urchin_IPM.m

urchin_IPM.m creates one of the curves that was layered on to Fig. 5.2, and creates Fig. 5.7. It calls makeSimpVec.m and kernmatSimp.m.

makeSimpVec.m creates a vector of size 1xn that, when vector multiplied by vector N, returns the integral of N with interval dy by Simpson's rule. Inherited from Marissa Baskett.

kernmatSimp.m creates an IPM kernel with a given integration mesh and parameter set. Adapted from code originally written by Easterling et al. (2001, Ecology) and adapted by White et al. (2016, Ecological Applications). Calls mkkern.m

mkkern.m creates an IPM kernel given the parameters assigned in kernmatSimp.m. Adapted from code originally written by Easterling et al. (2001, Ecology) and adapted by White et al. (2016, Ecological Applications). Calls mkkern.m

urchin_growth_mortality_plot.m creates Fig. 5.3a

red_urchin.m creates Fig. 5.3b,c, Fig. 5.4, and Fig. 5.5

dungeness_cohorts creates Fig. 5.6. The shading was added later in Adobe Illustrator. This code also creates a sixth panel with size-dependent mortality that was not included in the book.




