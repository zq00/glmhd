# glmhd

This R package estimates the bias and variance of the MLE from a high-dimensional binary regression model. 

You can install the package using 

```R
install.packages("devtools")
devtools::install_github("zq00/glmhd")
```

To get started, you can take a look at package vignettes (it takes about 20 min to knit the whole document, feel free to download the markdown file and run code line by line). To install with vignettes, run

```R
devtools::install_github("zq00/glmhd", build_vignettes = TRUE)
```

Alternatively, you can start by checking documentation for function `adjust_glm`. 

You can read more about the theory of high dimensional logistic MLE [here](https://arxiv.org/abs/2001.09351).

If you encounter errors or would like to provide feedbacks, please use Github -> Issues to reach us. Thank you! 
