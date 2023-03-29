# glmhd

## Overview

This R package estimates the inflation and variance of the MLE from a high-dimensional binary regression model. The package supports two methods: 

- High-dimensional theory (which provides an exact asymptotic distribution assuming the covariates are multivariate Gaussian or sub-Gaussian).

- Resized bootstrap methods (which approximates the MLE distribution but does not make distributional assumptions of the covariates)

## Getting started

- You can find the package **vignette** in the [Github -> vignettes](https://github.com/zq00/glmhd/tree/master/vignettes) folder. 
- You can read more about the theory of high dimensional logistic MLE, and the methods used in the package [here](https://arxiv.org/abs/2001.09351).  
- You can read more about the resized bootstrap method [here](https://arxiv.org/abs/2208.08944).

## Installation

You can install the package using 

```R
install.packages("devtools")
devtools::install_github("zq00/glmhd")
```

To install with vignettes, please run

```R
devtools::install_github("zq00/glmhd", build_vignettes = TRUE)
```

Note that the vignette takes about 20 min to knit, so feel free to download the [Rmarkdown file](https://github.com/zq00/glmhd/tree/master/vignettes) and run code line by line.

## Source code

The source code is located at the [Github -> R](https://github.com/zq00/glmhd/tree/master/R) folder. 

## Feedback

If you encounter error or would like to provide feedback, please use [Github -> Issues](https://github.com/zq00/glmhd/issues) to reach us. Thank you! 
