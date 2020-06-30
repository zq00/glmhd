# glmhd

## Overview

This R package estimates the bias and variance of the MLE from a high-dimensional binary regression model. 

## Getting started

The package vignette is under the Articles tab at this [link](https://zq00.github.io/glmhd/articles/my-vignette.html). You can read more about the theory of high dimensional logistic MLE, and the methods used in the package [here](https://arxiv.org/abs/2001.09351).  

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

Note that the vignette takes about 20 min to knit, so feel free to download the [Rmarkdown file](https://github.com/zq00/glmhd/tree/master/vignettes) and run code line by line..  

## Function documentation

You can find the function documentations under the [Reference](https://zq00.github.io/glmhd/reference/index.html) tab. To get started, you can take a look at the function `adjust_glm`, which computes the adjusted coefficient and standard error estimates.

## Feedback

If you encounter errors or would like to provide feedbacks, please use [Github -> Issues](https://github.com/zq00/glmhd/issues) to reach us. Thank you! 
