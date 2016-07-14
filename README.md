# Precision

This is the precision R package, to accompany the paper 'Detecting non-binomial sex allocation when developmental mortality operates' by Richard Wilkinson, Apostolos Kapranas and Ian Hardy.

## Installation

The easiest way to  install the package is to use devtools to install directly from github.

```
devtools::install_github('rich-d-wilkinson/precision')
```

## Data 

All of the datasets used in the paper are included in the R package. 

```
library(precision)
data(package='precision')$results[,c('Item')]
```

For example, the C. florus secondary dataset is

```
data(CflorusSecondary)
tail(CflorusSecondary)
```

To use your own dataset, specify a C x 2 matrix, with the first column containing the clutch size, and the second the number of males. It is necessary to label your columns as n and m.

```
my_data <- matrix(c(3,2,4,3,5,1,6,2,7,1,7,1), nc=2, byrow=TRUE)
colnames(my_data) <- c('n', 'm')
my_data
```

## Tutorial

See the package vignette for details on how to run the R code.