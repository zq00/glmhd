# glmhd

## Overview

This R package estimates the inflation and variance of the MLE from a high-dimensional binary regression model. (Note: This package is currently under construction to fix numerical inaccuracy when gamma or beta0 is large.)

## Getting started

- You can find the package **vignette** under the Articles tab (https://zq00.github.io/glmhd/articles/my-vignette.html). 
- You can read more about the theory of high dimensional logistic MLE, and the methods used in the package [here](https://arxiv.org/abs/2001.09351).  

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

## Function documentation

You can find the function documentations under the [Reference](https://zq00.github.io/glmhd/reference/index.html) tab. To get started, you can take a look at the function `adjust_glm`, which computes the adjusted coefficient and standard error estimates.

## Source code

The source code is located at the [Github -> R](https://github.com/zq00/glmhd/tree/master/R) folder. 

## Feedback

If you encounter error or would like to provide feedback, please use [Github -> Issues](https://github.com/zq00/glmhd/issues) to reach us. Thank you! 
