#Estimating dynamic treatment regimes for ordinal outcomes with household interference: Application in household smoking cessation

This repository is the official implementation of [Estimating dynamic treatment regimes for ordinal outcomes with household interference: Application in household smoking cessation](https://arxiv.org/abs/2306.12865).

## Introduction
The paper addresses a gap in precision medicine research, focusing on the estimation of dynamic treatment regimes (DTRs) with ordinal outcomes, especially in the presence of interference where one patient's treatment may affect another's outcome. We introduce the Weighted Proportional Odds Model (WPOM), a regression-based, doubly-robust approach for single-stage DTR estimation that considers ordinal outcomes and accounts for interference, particularly within households. The proposed method incorporates covariate balancing weights derived from joint propensity scores to address interference. Through simulation studies, the authors demonstrate the double robustness of WPOM with adjusted weights. The paper also extends WPOM to multi-stage DTR estimation with household interference and applies the methodology to analyze longitudinal survey data from the Population Assessment of Tobacco and Health study.

## Requirements

 - R 3.6
 - `simstudy`
 - `ordinal`
 - `mets`
 - `dplyr`
 - `datarium`
 - `caret`
 - `brant`
 
## Contents

  1. `README.txt`: implementation details of source code and contents 

  2. Three R folders: (1) Study1Singlestage, (2) Study2Multistage, correspond to Study 1 (Simulation 3.1) and Study 2 (Simulation 3.2) in the paper. (3) function contains  key funcitons in the project. 
      a). For single-stage decision settings, in Study 1a, we compare the WPOM with methods that ignore interference. In Study 1b, which is illustrated as Simulation 3.1 in the paper, we consider four different scenarios to verify the double robustness property of our method.
      b). For multi-stage decision settings, we aim to illustrate our estimation of a two-stage treatment decision problem with ordinal outcomes under household interference.
     
 

