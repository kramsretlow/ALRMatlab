# ALRMatlab

Autologistic regression modelling. MATLAB code supporting the paper "Better Autologistic Regression"

## Summary

This repository contains the MATLAB code and data files written in support of the paper:

Mark A. Wolters (2017), Better Autologistic Regression, *Frontiers in Applied Mathematics and Statistics,* 3(24), 1-20, https://doi.org/10.3389/fams.2017.00024.

Please cite this paper if you make use of this repository.

Autologistic regression is a statistical model for analysis of correlated binary responses.  For details, see the above-linked article.  [An online demo](http://www.mwolters.com/misc/2016/1114-ALRcoding.html) can also provide some background. Using the code in this repository, one can fit autologistic regression models and also explore the characteristics of such models, e.g. by drawing samples from them.  

While considerable effort has been made to ensure that the software is correct and useful, it is research-level software.  A more comprehensive and performant implementation is currently under development (using the [Julia](https://julialang.org/) language) at the time of writing.

## Brief Description of the Code 

The repository consists of two folders.  `ALRclasses` contains MATLAB OOP code implementing the autologistic models. It should be on your search path.

Create a model by instantiating an object and setting its members (parameters, covariates, responses, etc.). Parameter estimation can be done using matlab optimization functions.

The folder `PaperCode` contains scripts, functions, and data files used to do the analyses in the paper.  The scripts have filenames beginning with `Script_`. A good place to start is `Script_ExtraAnalysisForFrontiersFINAL.m`.

## License 

Copyright 2018 Mark Anthony Wolters.  This repository is under the MIT license.  See the license file.
