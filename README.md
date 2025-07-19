## Background
This repository contains the R programming scripts used for model development and implementation for the the paper entitled "A Bayesian statistical model-based approach for disaggregating small area population estimates by demographic characteristics". 
We present a flexible model-based appraoch for disaggregating small area population estimates by age, sex and social groups in order to facilitate more efficient target group-based interventions. The methodology is a multi-level Binomial regression model-based appraoch and parameter inference was based on Bayesian statisics, thus allowing easy uncertainty quantification.

## R scripts
There are two major R scripts contained in this repository 
- CMR_application_main_final.R: This includes the method application script which shows the step by step appraoch for the method implementation using datasets on Cameroon as a case study.
- simulation_study-main_final.R: The R scripts developed for the purpose of simulation study which aimed primarily to evaluate method performance over different sample sizes.

## Data
Please note we have only included the simulated data in this repository. The application data may be obtained from the National Statistical Office, Cameroon.

## Using the jollofR package
To facilitate the use of the methodology described in the paper, we have developed a user-friendly R package known as "JollofR version 0.3.0" which automatically disaggregates population estimates at both administrative uni and grid cell levels, along with the estimates of uncertainty. Please visit https://github.com/wpgp/jollofR to learn more about the "jollofR" package. 
