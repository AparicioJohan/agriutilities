---
title: '`agriutilities`: Utilities for Data Analysis in Agriculture' 
tags:
  - R package
  - Agriculture
  - Plant Breeding
  - Field trials
  - Experimental design
  - Single trial analysis
  - MET-analysis
authors:
  - name: Johan S. Aparicio
    orcid:
    affiliation: 1
  - name: Alexia N. Bornhorst
    orcid: 0009-0009-2512-5524
    affiliation: 1
affiliations:
 - name: The Alliance of Bioversity International and CIAT, Palmira, Colombia
   index: 1
date: 28 March 2023
bibliography: paper.bib
---

# Summary

`agriutilities` is an R package designed to make the analysis of field trials easier and more accessible for everyone working in plant breeding. It has a simple and intuitive interface to carry out Single Trial Analysis (STA) and Multi-Environmental Trial Analysis (MET), as well as, to estimate genetic parameters such as heritability, BLUEs, BLUPs, variance components and stability indices, among others. With this package you can perform analysis for Completely Randomized Design [@kuehl2001], Randomized Complete Block Design [@kuehl2001], Partially Replicated Designs [@Cullis2006], Alpha-Lattice Design [@kuehl2001], and Row-Column Designs.

Row-Column designs are more desirable for accounting for field spatial variations. Spatial Analysis is the most recommended when the row and column information is available. “Many factors combine to generate micro-environments that differ from plot to plot, strongly influencing yield and other traits. It is necessary to correct them when estimating treatment and/or genotypic effects” [@rodriguez2018]. 

Whether you are a beginner or an experienced user, `agriutilities` will help you to carry out complex analyses quickly and easily with confidence. With built-in functions for fitting Linear Mixed Models (LMM) and S3 methods for plotting results, `agriutilities` is the ideal choice for anyone who wants to save time and focus on interpreting their results. 

# Statement of need

In agriculture, breeding has been used to continuously improve genotypes that are used to produce the crops that feed humanity. An essential part of the breeding programs is the field trials, where experiments are performed to evaluate different characteristics and properties of crops. In plant breeding, it is necessary to identify which genotypes have the desired characteristics, for example, nutritional characteristics or tolerance to drought, so it is very important to develop crop varieties that are healthy, resistant, and have a higher yield. When we carry out field trials, they must have appropriate management so that the experimental error is minimal between plots, and thus the best genotypes are easily detected. Therefore, several questions can appear around how to model the different experimental designs, how to fit a spatial model, how to analyze multiple trials and what to do with multiple traits, along with others. 

In breeding programs, there are complex scenarios that the data analyst faces when it comes to analyze the data collected. Some of those scenarios occur when there are multiple trials in different locations and in different years. The measured traits are not found in all the trials; not all genotypes are found in the locations; there are different experimental designs between trials, among others. For this reason, the package `agriutilities` was designed to facilitate the analysis of field trials.

# Usage

When we have the data, the first thing that we must do is to identify the experimental design with the function check_design_met(), so we can check the quality of the data and identify the design. This works as a quality check before we fit any model. This function, in addition, allows the application of filters to the data, creates a summary of the traits and the experiments, creates a list with the experimental design detected, and prepares the data to be analyzed. For example, the summary of the traits by trial can provide us with information about the mean, median, variation, minimum and maximum of each trait by trial. Below is the code used to perform the data analysis pipeline.

``` r
# Data Analysis Pipeline

# (1) Installation 
install.packages("agriutilities")

# (2) Load library
library(agriutilities)

# (3) Read Data
datos <- readRDS("data/data_beans_VEF.RDS")
head(datos)

# (4) Check data
results <- check_design_met(
  data = datos,
  genotype = "Linea",
  trial = "trial",
  traits = c("YDHA", "DF", "DPM"),
  rep = "REP",
  block = "BLOCK",
  col = "col",
  row = "row"
)
results

# (5) Single Trial Analysis
out <- single_trial_analysis(results)
out

# (6) Multi-Environmental Trial Anaysis
met_results <- met_analysis(out)
met_results
```

Next, there is a function called single_trial_analysis() where the previous results are used to adjust single-trial models. The mixed modeling engine used in this function is either `lme4` or `ASRreml` [@asreml]. The software `ASReml` provides a great deal of flexibility and utility for the analysis of field trials in breeding programs. For spatial designs, `SpATS` is always used, and for other designs, `ASReml` is used as a default. `agriutilities` relies on the `statgenSTA` package for most of the single trial analysis.  

The single_trial_analysis() function removes extreme observations and generates BLUEs/BLUPs, heritability, and summaries by trial/trait. The summary of the fitted models shows the variation and heritability of each trait and the experimental design detected. Heritability describes how much variation in each trait can be attributed to genetic variation. You can also graph the coefficient of variation and heritability for each trait with the command plot, as shown in Figure 1.

![\label{fig:Fig}](summaryPlot.png)
<div align="center"> Fig 1. Summary of the fitted models. Coefficient of variation and heritability.</div>


The results of the single_trial_analysis() function are used in met_analysis() to fit multi-environmental trial models. Most field tests for plant breeding are replicated across different environments to measure the performance of breeding stocks across a range of environmental conditions to which a cultivar might be exposed [@isik2017genetic]. The function met_analysis() filters trials by heritability, fits GxE models by trait, estimates genetic correlations and overall heritability, calculates stability indices, and provides overall breeding values across locations. It integrates simpler variance covariances such as corv to more complex ones such as corgh or factor analytics.

In this package there are S3 methods that automatically generate a plot of the resulting objects check_desing_met, single_trial_analysis and met_analysis. This can generate a plot of missing data, connectivity, boxplot, and a summary of the variation and heritability. Also, there are plots of correlation, covariance, multi-traits, and spatial plots. 

# Availability 

This package can be downloaded directly from the R platform. It is also available in the [GitHub](https://github.com/AparicioJohan/agriutilities) repository. Some of the functions require the R package `ASReml` which can be obtained upon purchase from [VSN international](https://vsni.co.uk/software/asreml-r).

# References
