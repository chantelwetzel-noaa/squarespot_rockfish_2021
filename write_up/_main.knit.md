---
geometry: margin=1in
month: "December"
year: "2020"
preamble: |
output: sa4ss::techreport_pdf
bibliography: sa4ss.bib
lang: en
papersize: a4
---



<!--chapter:end:00a.Rmd-->

---
author:
  - name: Chantel R. Wetzel
    code: 1
    first: C
    middle: R
    family: Wetzel
  - name: Brian J. Langseth
    code: 1
    first: B
    middle: J
    family: Langseth
  - name: Jason M. Cope
    code: 1
    first: J
    middle: M
    family: Cope
  - name: John E. Budrick
    code: 2
    first: J
    middle: E
    family: Budrick
author_list: Wetzel, C.R., B.J. Langseth, J.M. Cope, J.E. Budrick
affiliation:
  - code: 1
    address: Northwest Fisheries Science Center, U.S. Department of Commerce, National
      Oceanic and Atmospheric Administration, National Marine Fisheries Service, 2725
      Montlake Boulevard East, Seattle, Washington 98112
  - code: 2
    address: California Department of Fish and Wildlife, 350 Harbor Boulevard, Belmont,
      California 94002
address:
  - ^1^Northwest Fisheries Science Center, U.S. Department of Commerce, National Oceanic
    and Atmospheric Administration, National Marine Fisheries Service, 2725 Montlake
    Boulevard East, Seattle, Washington 98112
  - ^2^California Department of Fish and Wildlife, 350 Harbor Boulevard, Belmont,
    California 94002
---

<!--chapter:end:00authors.Rmd-->

---
title: Squarespot Rockfish (_Sebastes hopkinsi_) along the California U.S. West Coast
  in 2020
---

<!--chapter:end:00title.Rmd-->

\pagebreak
\pagenumbering{roman}
\setcounter{page}{1}

<!--chapter:end:01a.Rmd-->


<!--chapter:end:01executive.Rmd-->

\pagebreak
\setlength{\parskip}{5mm plus1mm minus1mm}
\pagenumbering{arabic}
\setcounter{page}{1}
\renewcommand{\thefigure}{\arabic{figure}}
\renewcommand{\thetable}{\arabic{table}}

<!--chapter:end:10a.Rmd-->

# Introduction
## Basic Information


## Life History
Replace text.


## Historical and Current Fishery Information
Replace text.

## Summary of Management History and Performance
Replace text.



<!--chapter:end:11introduction.Rmd-->

# Data


<!--chapter:end:20data.Rmd-->

## Fishery-Dependent Data

<!--chapter:end:21f-.Rmd-->

## Fishery-Independent Data

<!--chapter:end:21s-.Rmd-->

### \acrlong{s-hkl}



<!--chapter:end:21s-hkl.RMd-->

### \acrlong{s-wcgbt}

The \Gls{s-wcgbt} is based on a random-grid design;
covering the coastal waters from a depth of 55-1,280 m [@bradburn_2003_2011].
This design generally uses four industry-chartered vessels per year assigned to a roughly equal number of randomly selected grid cells and divided into two 'passes' of the coast.
Two vessels fish from north to south during each pass between late May to early October.
This design therefore incorporates both vessel-to-vessel differences in catchability,
as well as variance associated with selecting a relatively small number (approximately 700) of possible cells from a very large set of possible cells spread from the Mexican to the Canadian borders.

<!--chapter:end:21s-wcgbts.Rmd-->

## Biological Data

### Natural Mortality



### Length-Weight Relationship



### Growth (Length-at-Age)



### Maturation and Fecundity



### Sex Ratio



<!--chapter:end:22biology.Rmd-->

# Assessment Model

## Summary of Previous Assessments 



### Bridging Analysis




## Model Structure and Assumptions

Squarespot rockfish area assessed using a two-sex model with sex specific life history parameters.   

### Modeling Platform and Structure

Stock Synthesis version 3.30.16 was used to estimate the parameters in the model. The R package r4ss, version 1.38.0, along with R version 4.0.1 were used to investigate and plot model fits. 

### Priors



### Data Weighting



### Estimated and Fixed Parameters



## Model Selection and Evaluation

The base assessment model for squarespot rockfish was developed to balance parsimony and realism, and the goal was to estimate a spawning output trajectory for the population of squarespot rockfish off the California coast. The model contains many assumptions to achieve parsimony and uses many different sources of data to estimate reality. A series of investigative model runs were done to achieve the final base model.



<!--chapter:end:30model.Rmd-->

## Base Model Results

 


### Parameter Estimates




### Fits to the Data




### Population Trajectory




### Reference Points



<!--chapter:end:33results.Rmd-->

## Model Diagnostics

### Convergence



### Sensitivity Analyses





### Retrospective Analysis

A five-year retrospective analysis was conducted by running the model using data only through 2015, 2016, 2017, 2018, 2019 and 2020. 

### Likelihood Profiles

Likelihood profiles were conducted for $R_0$, steepness, maximum length, and female natural mortality values separately. These likelihood profiles were conducted by fixing the parameter at specific values and estimated the remaining parameters based on the fixed parameter value.



### Unresolved Problems and Major Uncertainties


<!--chapter:end:34diagnostics.Rmd-->

# Management 

## Reference Points

## Unresolved Problems and Major Uncertainties

## Harvest Projections and Decision Tables

## Evaluation of Scientific Uncertainty

## Research and Data Needs

<!--chapter:end:40management.Rmd-->

# Acknowledgments


\newpage

<!--chapter:end:41acknowledgments.Rmd-->

\clearpage

# References
<!-- If you want to references to appear somewhere before the end, add: -->
<div id="refs"></div>
<!-- where you want it to appear -->
<!-- The following sets the appropriate indentation for the references
  but it cannot be used with bookdown and the make file because it leads
  to a bad pdf.
\noindent
\vspace{-2em}
\setlength{\parindent}{-0.2in}
\setlength{\leftskip}{0.2in}
\setlength{\parskip}{8pt}
 -->

<!--chapter:end:49bibliography.Rmd-->

\clearpage
# Tables

<!-- ======================================================= -->
<!-- ***************    Catches      *********************** --> 
<!-- ======================================================= -->




<!-- ======================================================= -->
<!-- ***************       Data      *********************** --> 
<!-- ======================================================= -->



<!-- ======================================================= -->
<!-- ***************    Biology      *********************** --> 
<!-- ======================================================= -->



<!-- ======================================================= -->
<!-- ***********   Model Parameters     ******************** --> 
<!-- ======================================================= -->


<!-- ======================================================= -->
<!-- ***********       Time Series      ******************** --> 
<!-- ======================================================= -->



<!-- ======================================================= -->
<!-- ****************     Sensitivities      *************** --> 
<!-- ======================================================= -->



<!-- ======================================================= -->
<!-- ********  Reference Points & Management *************** --> 
<!-- ======================================================= -->


<!--chapter:end:52tables.Rmd-->

\clearpage
# Figures


<!-- ====================================================================== -->  
<!-- ****************** Catches Used in the Model ************************* --> 
<!-- ====================================================================== -->  



<!-- ====================================================================== --> 
<!-- ******************* Data Used in the Model *************************** --> 
<!-- ====================================================================== --> 


![Summary of data sources used in the base model.\label{fig:data-plot}](C:/Assessments/2021/squarespot_rockfish_2021/models/2.3_hessian/plots/data_plot.png){width=100% height=100% alt="Summary of data sources used in the base model."}


<!-- ====================================================================== -->
<!-- ****************   Commercial Length Samples    ********************** --> 
<!-- ====================================================================== -->


<!-- ====================================================================== -->
<!-- **************** Recreational Length Samples    ********************** --> 
<!-- ====================================================================== -->




<!-- ====================================================================== -->
<!-- *************************     Biology     **************************** --> 
<!-- ====================================================================== -->


```

<!-- ====================================================================== -->
<!-- *********************    Selectivity            ********************** --> 
<!-- ====================================================================== -->



<!-- ====================================================================== -->
<!-- *********************   Recruitment     ****************************** --> 
<!-- ====================================================================== -->



<!-- ====================================================================== -->
<!-- ****************** Fit to the Length Data **************************** --> 
<!-- ====================================================================== -->



<!-- ====================================================================== -->
<!-- ******************      Time Series       **************************** --> 
<!-- ====================================================================== -->


![Estimated time series of spawning output.\label{fig:ssb}](C:/Assessments/2021/squarespot_rockfish_2021/models/2.3_hessian/plots/ts7_Spawning_output_with_95_asymptotic_intervals_intervals.png){width=100% height=100% alt="Estimated time series of spawning output."}


![Estimated time series of total biomass.\label{fig:tot-bio}](C:/Assessments/2021/squarespot_rockfish_2021/models/2.3_hessian/plots/ts1_Total_biomass_(mt).png){width=100% height=100% alt="Estimated time series of total biomass."}


![Estimated time series of fraction of unfished spawning output.\label{fig:depl}](C:/Assessments/2021/squarespot_rockfish_2021/models/2.3_hessian/plots/ts9_Fraction_of_unfished_with_95_asymptotic_intervals_intervals.png){width=100% height=100% alt="Estimated time series of fraction of unfished spawning output."}


![Stock-recruit curve. Point colors indicate year, with warmer colors indicating earlier years and cooler colors in showing later years.\label{fig:bh-curve}](C:/Assessments/2021/squarespot_rockfish_2021/models/2.3_hessian/plots/SR_curve.png){width=100% height=100% alt="Stock-recruit curve. Point colors indicate year, with warmer colors indicating earlier years and cooler colors in showing later years."}


<!-- ====================================================================== -->
<!-- ******************       Sensitivity     ***************************** --> 
<!-- ====================================================================== -->



<!-- ====================================================================== -->
<!-- ******************     Retrospectives    ***************************** --> 
<!-- ====================================================================== -->


<!-- ====================================================================== -->
<!-- ******************      Likelihoods      ***************************** --> 
<!-- ====================================================================== -->



<!-- ====================================================================== -->
<!-- ******************    Reference Points    **************************** --> 
<!-- ====================================================================== -->



![Estimated 1 - relative spawning ratio (SPR) by year.\label{fig:1-spr}](C:/Assessments/2021/squarespot_rockfish_2021/models/2.3_hessian/plots/SPR2_minusSPRseries.png){width=100% height=100% alt="Estimated 1 - relative spawning ratio (SPR) by year."}


![Equilibrium yield curve for the base case model. Values are based on the 2020
fishery selectivity and with steepness fixed at 0.72.\label{fig:yield}](C:/Assessments/2021/squarespot_rockfish_2021/models/2.3_hessian/plots/yield2_yield_curve_with_refpoints.png){width=100% height=100% alt="Equilibrium yield curve for the base case model. Values are based on the 2020
fishery selectivity and with steepness fixed at 0.72."}


\newpage

<!--chapter:end:53figures.Rmd-->

\clearpage
# Appendix A. Detailed Fit to Length Composition Data 



<!--chapter:end:54appendix.Rmd-->

