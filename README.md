# mirahelsinki

This page will contain R code used in statistical analyses and building the figures of MIRA Helsinki -study. You may contact the main author of the code via email: _topi.hovinen@helsinki.fi_

# mira_statistics_figures.R

This is the main R-file used. It contains most of the statistical analysis and the code to create figures/ggplots. Please note, that some of the files required for running the code smoothly are missing here. Please contact the author or otherwise specified instance below if you wish to use these files:

"HUSLAB Data final.txt" AND "MIRA Metabolomics final.txt": Due to small sample size and risk of identifying participants of the study, we choose not to publish the original data files online.

"thl_kasvukayrat.txt": The original data of Finnish growth curves were provided for the scientific purposes of this study, and as we do not own the rights for this data, please contact THL (Terveyden ja Hyvinvoinnin Laitos) in order to gain permission for this data.

# p.values.adjusted.RData

As running the permutation tests may be time consuming, here are the p-values calculated during our study, if one wishes to run other parts of code that requires the p-values (e.g. plotting significancy stars with permira.test.gg)

# WHO_MUAC_z-scores_(...).txt

WHO tables for Mid-Upper Arm Circumference (MUAC) z-scores and extensions by Mramba et al. used in this study for reference.
