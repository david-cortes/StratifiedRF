# StratifiedRF

R package producing a Random Forest-like ensemble of C5.0 decision trees that samples variables in groups independently.

In a traditional random forest, each tree is built by sampling a number of variables at random, but for some datasets that have variables containing information about different things in groups of too different sizes (e.g. data with information about a user and a product in the same row, with 500 variables about the product but only 10 about the user), it might be desirable to force the trees to always consider at least some variables from one group.

The trees are built using the C5.0 algorithm rather than the usual CART, thus it works for classification only (i.e. doesn’t perform regression). Supports missing values in predictors and categorical variables natively, without doing one-hot encoding.

## Instalation

Package is available on CRAN, can be installed with

```install.packages(“StratifiedRF”)```

## Usage
```
data(iris)
groups <- list(c("Sepal.Length","Sepal.Width"),c("Petal.Length","Petal.Width"))
mtry <- c(1,1)
m <- stratified_rf(iris,"Species",groups,mtry,ntrees=2,multicore=FALSE)
summary(m)
predict(m,iris)
```

## Documentation

All functions are documented internally (e.g. you can try `??StratifiedRF`). You can find the full documentation as PDF in the [CRAN link](https://CRAN.R-project.org/package=StratifiedRF).

## Implementation notes

The individual trees are built using the [C50 pacakge](https://CRAN.R-project.org/package=C50). Implementation of everything outside of tree-building is in native R code, thus speed is not great.
