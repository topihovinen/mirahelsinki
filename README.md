# mirahelsinki

This page will contain R code used in statistical analyses and building the figures of MIRA Helsinki -study. You may contact the main author of the code via email: _topi.hovinen@helsinki.fi_

## mira_statistics_figures.R - Main code

This is the main R-file used. It contains the statistical analysis of blood and urine samples and the code to create figures/ggplots. Please note, that some of the files required for running the code smoothly are missing here. Please contact the author or otherwise specified instance below if you wish to use these files:

*HUSLAB Data final.txt* AND *MIRA Metabolomics final.txt*: Due to small sample size and risk of identifying participants, we choose not to publish the original data files online.

*thl_kasvukayrat.txt*: The original data of Finnish growth curves were provided for the scientific purposes of this study, and as we do not own the rights for this data, please contact THL (Terveyden ja Hyvinvoinnin Laitos) in order to gain permission for this data.

### permira.test.R

The source code for functions ``permira.test`` and ``permira.test.gg``. The former is the actual permutation test algorithm and the latter is a dirty ad-hoc plug-in for ggplot to plot p-values from ``permira.test`` in ggplot with the ``geom_signif`` -function.

### p.values.adjusted.RData

As running the permutation tests may be time consuming, here are the Benjamini-Hochberg adjusted p-values calculated during our study, if one wishes to run other parts of code that requires the p-values (e.g. plotting significancy stars with ``permira.test.gg``)

### WHO_MUAC_z-scores_(...).txt

WHO tables for Mid-Upper Arm Circumference (MUAC) z-scores and extensions by Mramba et al. used in this study for reference.

## mira_nutrients_analysis.R - Code for nutrient intake analysis

This contains the statistical analysis of nutrient intake based on food frequency questionnaires and food diaries. The code requires ``permira.test`` -function that is found in another file here.

*muuttujat_analyysiin.txt*: This would be the data input file of nutrient intake. Due to small sample size and risk of identifying participants, we choose not to publish the original data files online.

## mira_pathway_analysis.R - Code for pathway analysis

This file includes the code used in pathway analysis. All the required files to run this code are included in this project folder

### pathwaydata.txt

This is the untargeted metabolomics data used in pathway analysis. Includes e.g. fold changes and p-values of each ion recognized in untargeted metabolomics. Comparison is done between omnivores and vegans.

### hmdbPtw_v3.0.mo

A text file including Human Metabolome DataBase pathways v.3.0, used in Pathway analysis.
