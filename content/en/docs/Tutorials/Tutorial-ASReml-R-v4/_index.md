---
title: "Practical: Animal models in ASReml-R v4 Accompanying Wilson et al. 2010. An Ecologist’s Guide to the Animal Model"
linkTitle: "Animal models in ASReml-R"
weight: 2
description: >
  Using ASReml-R
output: 
  html_document: 
    keep_md: yes
---

Erik Postma - e.postma@exeter.ac.uk

## Tutorial 1 - Estimating the heritability of birth weight

This tutorial will demonstrate how to run a univariate animal model using the software ASReml–R and example data files provided.

### 1) Scenario

In a population of gryphons there is strong positive selection on birth weight with heavier born individuals having, on average higher fitness. To find out whether increased birth weight will evolve in response to the selection, and if so how quickly, we want to estimate the heritability of birth weight.

### 2) Data files

Open `gryphonped.txt` and `gryphon.txt` in your text editor. The structure and contents of these files is fairly self-explanatory. 
The pedigree file `gryphonped.txt` contains three columns containing unique IDs that correspond to each animal, its father, and its mother.
Note that this is a multigenerational pedigree, with the earliest generation (for which parentage information is necessarily missing) at the beginning of the file.
For later-born individuals maternal identities are all known but paternity information is incomplete (a common situation in real world applications).
The phenotype data, as well as additional factors and covariates that we may wish to include in our model are contained in `gryphon.txt`.
Columns correspond to individual identity (`ANIMAL`), maternal identity (`MOTHER`), year of birth (`BYEAR`), sex (`SEX`, where `1` is female and `2` is male), birth weight (`BWT`), and tarsus length (`TARSUS`).
Each row of the data file contains a record for a different offspring individual. Note that all individuals included in the data file must be included as offspring in the pedigree file.

### 3) Running the model

To run an analysis line by line, copy the code below into the R console.
Alternatively, you can create a new R script. Note that the file pathway (i.e. `C://...` or `/users/...`) will have to be changed to direct R to the location of the tutorial data files on your computer.
This tutorial assumes you have already installed `ASReml-R version 4`. See the installation instructions for more details on how to do this.

First we need to load the `asreml` library:

```r
library(asreml)
```

```
### Loading required package: Matrix
```

Now we can read in the data file, telling `R` that `NA` is the symbol for missing values and that the first line of the file contains the column headers:

``` r
gryphon <- read.table("/the/path/to/gryphon.txt", 
                header=T, na.strings="NA")
```

It is a good idea to make sure that all variables are correctly assigned as numeric or factors:

``` r
gryphon$ANIMAL<-as.factor(gryphon$ANIMAL)
gryphon$MOTHER<-as.factor(gryphon$MOTHER)
gryphon$BYEAR<-as.factor(gryphon$BYEAR)
gryphon$SEX<-as.factor(gryphon$SEX)
gryphon$BWT<-as.numeric(gryphon$BWT)
gryphon$TARSUS<-as.numeric(gryphon$TARSUS)
```

Similarly we can read in the pedigree file, again telling `R` that `NA` is the symbol for missing values and that the first line of the file contains the column headers.

``` r
gryphonped <- read.table("/the/path/to/gryphonped.txt", 
              header=T, na.strings="NA")

gryphonped$ID<-as.factor(gryphonped$ID)
gryphonped$FATHER<-as.factor(gryphonped$FATHER)
gryphonped$MOTHER<-as.factor(gryphonped$MOTHER)
``` 

Now that we have imported the data and the pedigree file, we are ready to fit an animal model. To do this, we first have use our pedigree data to create (the inverse of) the relationship matrix:

``` r
ainv <- ainverse(gryphonped)
```

We are now ready to specify our first model:

``` r
model1 <- asreml(fixed = BWT ~ 1, random =~vm(ANIMAL, ainv), 
                 residual=~idv(units),
                 data=gryphon, 
                 na.action = na.method(x="omit", y="omit"))
``` 

``` 
## License check Wed Sep 16 09:47:29 2020
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:31 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -4128.454           1.0    853 09:47:31    0.0
##  2     -3284.272           1.0    853 09:47:31    0.0
##  3     -2354.992           1.0    853 09:47:31    0.0
##  4     -1710.357           1.0    853 09:47:31    0.0
##  5     -1363.555           1.0    853 09:47:31    0.0
##  6     -1263.516           1.0    853 09:47:31    0.0
##  7     -1247.854           1.0    853 09:47:31    0.0
##  8     -1247.185           1.0    853 09:47:31    0.0
##  9     -1247.183           1.0    853 09:47:31    0.0
``` 

In this model, `BWT` is the response variable and the only fixed effect is the mean (the intercept, denoted as `1`).
The only random effect we have fitted is `ANIMAL`, which will provide an estimate of $V_A$. Our random ANIMAL effect is connected to the inverse related matrix `ainv`.
`data=` specifies the name of the dataframe that contains our variables. Finally, we tell `asreml()` what to when it encounters NAs in either the dependent or predictor variables (in this case we choose to remove the records).

A note of the specification of the structure of the residuals: This simple univariate model will run fine without `residual=~idv(units)`.
However, if you are going to use `vpredict()` to calculate the heritability (see below), not specifying the residuals in this way will result in a standard error for the heritability that is incorrect.

To see the estimates for the variance components, we run:

``` r
summary(model1)$varcomp
``` 

``` 
##                  component std.error  z.ratio bound %ch
## vm(ANIMAL, ainv)  3.395398 0.6349915 5.347154     P   0
## units!units       3.828602 0.5185919 7.382687     P   0
## units!R           1.000000        NA       NA     F   0
``` 

We fitted a single random effect so have partitioned the phenotypic variance into two components.
The `vm(ANIMAL, ainv)` variance component is $V_A$ and is estimated as `3.4`.
Given that the ratio of $V_A$ to its standard error (`z.ratio`) is considerably larger than 2 (i.e. the parameter estimate is more than 2 SEs from zero) this looks likely to be highly significant.
The `units!units` component refers to the residual variance $V_R$, and `units$R` should be ignored.
If you don’t include `residual=~idv(units)` in your model specification, `units$R` will provide you with the residual variance.

## 4) Estimating heritability

We can calculate the $h^2$ of birth weight from the components above since $h^2=V_A/V_P=V_A/(V_A+V_R)$. Thus according to this model, $h^2 = 3.4 / (3.4 + 3.83) = 0.47$.

Alternatively we can use the `vpredict()` function to calculate h2h2 and its standard error:


``` r
vpredict(model1, h2.BWT ~ V1 / (V1 + V2))
``` 
``` 
##         Estimate         SE
## h2.BWT 0.4700163 0.07650881
``` 

### 5) Adding fixed effects

To add fixed effects to a univariate model simply modify the model statement. For example we might know (or suspect) that birth weight is a sexually dimorphic trait and therefore fit a model

``` r
model2 <- asreml(fixed = BWT ~ 1 + SEX, 
             random= ~vm(ANIMAL, ainv), 
             residual=~idv(units),
             data=gryphon, 
             na.action = na.method(x="omit", y="omit"))
``` 
``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:31 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -3364.126           1.0    852 09:47:31    0.0
##  2     -2702.117           1.0    852 09:47:31    0.0
##  3     -1978.916           1.0    852 09:47:31    0.0
##  4     -1487.834           1.0    852 09:47:31    0.0
##  5     -1236.350           1.0    852 09:47:31    0.0
##  6     -1172.771           1.0    852 09:47:31    0.0
##  7     -1165.270           1.0    852 09:47:31    0.0
##  8     -1165.093           1.0    852 09:47:31    0.0
##  9     -1165.093           1.0    852 09:47:31    0.0
``` 

Now we can look at the fixed effects parameters and assess their significance with a conditional Wald F-test:

``` r
summary(model2, coef = TRUE)$coef.fixed
wald.asreml(model2, ssType="conditional", denDF="numeric")
``` 

``` 
##             solution std error  z.ratio
## SEX_1       0.000000        NA       NA
## SEX_2       2.206996 0.1619974 13.62365
## (Intercept) 6.058669 0.1718244 35.26082

## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:31 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -1165.093           1.0    852 09:47:31    0.0
##  2     -1165.093           1.0    852 09:47:31    0.0
## Calculating denominator DF

## 
##             Df denDF  F.inc  F.con Margin          Pr
## (Intercept)  1   251 3491.0 3491.0        0.00000e+00
## SEX          1   831  185.6  185.6      A 2.70204e-38
``` 

The very small probability (`Pr`) in the Wald test above shows that `SEX` is a highly significant fixed effect, and from the parameter estimates we can see that the average male (sex 2) is 2.2 kg ($\pm$ 0.16 SE) heavier than the average female (sex 1).
However, when we look at the variance components in the model including `SEX` as a fixed effect, we see that they have changed slightly from the previous model:

``` r
summary(model2)$varcomp
``` 

``` 
##                  component std.error  z.ratio bound %ch
## vm(ANIMAL, ainv)  3.060441 0.5243571 5.836558     P   0
## units!units       2.938412 0.4161473 7.060991     P   0
## units!R           1.000000        NA       NA     F   0
``` 

In fact since `SEX` effects were previously contributing to the residual variance of the model, our estimate of $V_R$ (denoted `units!R` in the output) is now slightly lower than before.
This has an important consequence for estimating heritability since if we calculate $V_P$ as $V_A+V_R$ then as we include fixed effects we will soak up more residual variance driving $V_P$.
Assuming that $V_A$ is more or less unaffected by the fixed effects fitted then as $V_P$ goes down we expect our estimate of $h^2$ will go up:

``` r
vpredict(model2, h2.BWT ~ V1 / (V1 + V2))
``` 
``` 
##        Estimate         SE
## h2.BWT 0.510171 0.07432388
``` 
Here $h^2$ has increased slightly from `0.47` to `0.51`. Which is the better estimate? It depends on what your question is.
The first is an estimate of the proportion of variance in birth weight explained by additive effects, the latter is an estimate of the proportion of variance in birth weight after conditioning on sex that is explained by additive effects.

### 6) Adding random effects

This is done by simply modifying the model statement in the same way. For instance fitting

``` r
model3 <- asreml(fixed = BWT ~ 1 + SEX, 
             random= ~vm(ANIMAL, ainv) + BYEAR, 
             residual=~idv(units),
             data=gryphon, 
             na.action = na.method(x="omit", y="omit"))
``` 

``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:31 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -2742.658           1.0    852 09:47:31    0.0
##  2     -2237.268           1.0    852 09:47:31    0.0
##  3     -1690.453           1.0    852 09:47:31    0.0
##  4     -1328.910           1.0    852 09:47:31    0.0
##  5     -1154.597           1.0    852 09:47:31    0.0
##  6     -1116.992           1.0    852 09:47:31    0.0
##  7     -1113.809           1.0    852 09:47:31    0.0
##  8     -1113.772           1.0    852 09:47:31    0.0
##  9     -1113.772           1.0    852 09:47:31    0.0
``` 

results in an additional variance component of birth year:

``` r
summary(model3)$varcomp
``` 

``` 
##                  component std.error  z.ratio bound %ch
## BYEAR            0.8862604 0.2695918 3.287416     P   0
## vm(ANIMAL, ainv) 2.7068665 0.4422140 6.121169     P   0
## units!units      2.3092415 0.3451025 6.691466     P   0
## units!R          1.0000000        NA       NA     F   0
``` 

Here the variance in `BWT` explained by `BYEAR` is `0.886` and, based on the `z.ratio`, appears to be significant.
Thus we would conclude that year-to-year variation (e.g., in weather, resource abundance) contributes to $V_P$.
Note that although $V_A$ has changed somewhat, as most of what is now partitioned as a birth year effect was previously partitioned as $V_R$.
Thus what we have really done here is to partition environmental effects into those arising from year-to-year differences versus everything else, and we do not really expect much change in $h^2$ (since now $h^2=V_A/(V_A+V_{BY}+V_R)$).

However, we get a somewhat different result if we also add a random effect of `MOTHER` to test for maternal effects:

``` r
model4 <- asreml(fixed = BWT ~ 1 + SEX, 
             random= ~vm(ANIMAL, ainv) + BYEAR + MOTHER, 
             residual=~idv(units),
             data=gryphon, 
             na.action = na.method(x="omit", y="omit"))
``` 

``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:31 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -2033.178           1.0    852 09:47:31    0.0
##  2     -1723.734           1.0    852 09:47:31    0.0
##  3     -1396.354           1.0    852 09:47:31    0.0
##  4     -1193.012           1.0    852 09:47:31    0.0
##  5     -1107.946           1.0    852 09:47:31    0.0
##  6     -1095.327           1.0    852 09:47:31    0.0
##  7     -1094.816           1.0    852 09:47:31    0.0
##  8     -1094.815           1.0    852 09:47:31    0.0
``` 

Gives estimated variance components of

``` r
summary(model4)$varcomp
``` 
``` 
##                  component std.error  z.ratio bound %ch
## BYEAR            0.8820313 0.2632455 3.350604     P   0
## MOTHER           1.1184698 0.2386239 4.687167     P   0
## vm(ANIMAL, ainv) 2.2985320 0.4962496 4.631806     P   0
## units!units      1.6290034 0.3714154 4.385934     P   0
## units!R          1.0000000        NA       NA     F   0
``` 

Here partitioning of significant maternal variance has resulted in a further decrease in $V_R$ but also a decrease in $V_A$.
The latter is because maternal effects of the sort we simulated (fixed differences between mothers) will have the consequence of increasing similarity among maternal siblings.
Consequently they can look very much like additive genetic effects and if present, but unmodelled, represent a type of "common environment effect" that can - and will - cause upward bias in $V_A$ and so $h^2$.

### 7) Testing significance of random effects

A final point to note in this tutorial is that while the z.ratio (component/std.error) reported is a good indicator of likely statistical significance (>1.96?), the standard errors are approximate and are not recommended for formal hypothesis testing. A better approach is to use likelihood-ratio tests.

For example, to test the significance of maternal effects we could compare models with and without the inclusion of maternal identity as a random effect and compare the final log-likelihoods of these models.

``` r
model4$loglik
``` 

``` 
## [1] -1094.815
``` 

shows that the model including maternal identity has a log-likelihood of `-1094.815`, and

``` r
model3$loglik
``` 
``` 
## [1] -1113.772
``` 
shows that the model excluding maternal identity has a log-likelihood of `-1113.772`.

A test statistic equal to twice the absolute difference in these log-likelihoods is assumed to be distributed as Chi square with one degree of freedom.
In this case we would conclude that the maternal effects are highly significant since: `2 * (-1094.8145793 - -1113.7719147)` = `37.9146708`, and the p-value that comes with this is:

``` r
1 - pchisq(2 * (model4$loglik - model3$loglik),1)
``` 
``` 
## [1] 7.390738e-10
``` 
As P < 0.0001 we would therefore conclude that the additional of maternal identity as a random effect significantly improves the model, given an increase in log-likelihood of approximately 19.

## Tutorial 2 – A bivariate animal model

This tutorial will demonstrate how to run a multivariate animal model using the software ASReml-R and example data files provided.

### 1) Scenario

Since natural selection rarely acts on single traits, to understand how birth weight might evolve in our population of gryphons, we may also want to think about possible covariance with other traits.
If tarsus length at fledging is also under positive selection, what implications does this have for birth weight and vice versa?
If the two traits are positively genetically correlated then this will facilitate evolution of larger size (since response of one trait will induce a positively correlated response in the other).
If there is negative genetic covariance then this could act as an evolutionary constraint.

Using multivariate models allows the estimation of parameters relating to each trait alone (i.e. $V_A$, $h^2$, etc), but also yields estimates of covariance components between traits.
These include the (additive) genetic covariance $COV_A$ which is often rescaled to give the additive genetic correlation $r_A$.
However, covariance can also arise through other random effects (e.g. maternal covariance) and these sources can also be explicitly modelled in a bivariate analysis.

### 2) Data files

Pedigree and phenotypic data files are the same as those used in tutorial 1 (i.e, `gryphonped.txt` and `gryphon.txt` respectively).

### 3) Running the model

For running multivariate analyses in `ASReml-R`, the code is slightly more complex than for the univariate case.
This is because `ASReml-R` allows us to make different assumptions about the way in which traits might be related and so we need to explicitly code a model of the (co)variance structure we want to fit.
We have also specified some starting values. These are can be very approximate guestimates, but having reasonable starting values can aid convergence.
Finally, we have increased the default maximum number of iterations (`maxiter`) which can help to achieve convergence for more complicated models.

``` r
modela<-asreml(fixed = cbind(BWT, TARSUS) ~ trait,
               random=~ us(trait):vm(ANIMAL, ainv), init=c(1, 0.1, 1),
               residual= ~ id(units):us(trait, init=c(1, 0.1, 1)),
               data=gryphon,
               maxit=20)
``` 

``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:32 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -5118.122           1.0   1535 09:47:32    0.0
##  2     -4358.769           1.0   1535 09:47:32    0.0
##  3     -3540.792           1.0   1535 09:47:32    0.0
##  4     -3004.970           1.0   1535 09:47:32    0.0
##  5     -2747.831           1.0   1535 09:47:32    0.0
##  6     -2687.807           1.0   1535 09:47:32    0.0
##  7     -2680.057           1.0   1535 09:47:32    0.0
##  8     -2679.743           1.0   1535 09:47:32    0.0
##  9     -2679.741           1.0   1535 09:47:32    0.0
``` 

This has fitted a bivariate model of `BWT` and `TARSUS`, with the mean for each of the traits as a fixed effect (trait).
The additive genetic variance-covariance matrix ($\mathbf{G}$) is unstructured (us; i.e. all elements are free to vary) and the starting values for $V_A$ for `BWT`, $COV_A$ between `BWT` and `TARSUS`, and $V_A$ for `TARSUS` are set to 1, 0.1 and 1, respectively. Similarly, the residual matrix is unstructured and uses the same starting values.

Let’s have a look at the variance components, and notice that there are now six (co)variance components reported in the table:

``` r
summary(modela)$varcomp
``` 
``` 
##                                            component std.error  z.ratio bound
## trait:vm(ANIMAL, ainv)!trait_BWT:BWT        3.368498 0.6349441 5.305188     P
## trait:vm(ANIMAL, ainv)!trait_TARSUS:BWT     2.460195 1.0735865 2.291567     P
## trait:vm(ANIMAL, ainv)!trait_TARSUS:TARSUS 12.346837 3.0753775 4.014739     P
## units:trait!R                               1.000000        NA       NA     F
## units:trait!trait_BWT:BWT                   3.849836 0.5199897 7.403677     P
## units:trait!trait_TARSUS:BWT                3.312978 0.9128834 3.629136     P
## units:trait!trait_TARSUS:TARSUS            17.645584 2.6668492 6.616641     P
##                                            %ch
## trait:vm(ANIMAL, ainv)!trait_BWT:BWT       0.0
## trait:vm(ANIMAL, ainv)!trait_TARSUS:BWT    0.2
## trait:vm(ANIMAL, ainv)!trait_TARSUS:TARSUS 0.1
## units:trait!R                              0.0
## units:trait!trait_BWT:BWT                  0.0
## units:trait!trait_TARSUS:BWT               0.1
## units:trait!trait_TARSUS:TARSUS            0.1
``` 

The first three terms relate to the genetic matrix and, in order are $V_{A,BWT}$, $COV_A$, $V_{A,TARSUS}$.
Below is again a line where the `units:traitr!R` component equals 1, which again can be ignored.
The final three terms relate to the residual matrix and correspond to $V_{R,BWT}$, $COV_R$, $V_{R,TARSUS}$. Based on our quick and dirty check (is z.ratio > 1.96?) all components look to be statistically significant.

We can calculate the genetic correlation as $COV_A/\sqrt{V_{A,BWT}V_{A,TARSUS}}$.
Thus this model gives an estimate of rArA​ = 0.38. Although we can calculate this by hand, we can also use `vpredict()`, which also provides an (approximate) standard error:

``` r
vpredict(modela, rA ~ V2/sqrt(V1*V3))
``` 
``` 
##     Estimate       SE
## rA 0.3814816 0.129991
``` 
Of course we can also calculate the heritability of BWT and TARSUS from this model:

``` r
vpredict(modela, h2.bwt ~ V1/(V1+V5))
``` 

``` 
##         Estimate         SE
## h2.bwt 0.4666586 0.07672219
``` 

``` r
vpredict(modela, h2.tarsus ~ V3/(V3+V7))
``` 
``` 
##            Estimate        SE
## h2.tarsus 0.4116652 0.0930738
``` 

### 4) Adding fixed and random effects

Fixed and random effects can be added just as for the univariate case.
Given that our full model of BWT from tutorial 1 had SEX as a fixed effect as well as random effects of BIRTH YEAR and MOTHER we could specify a bivariate formulation of this using:

``` r
modelb <- asreml(fixed = cbind(BWT, TARSUS) ~ trait + trait:SEX,
               random=~ us(trait, init=c(1,0.1,1)):vm(ANIMAL, ainv) +
                 us(trait, init=c(1,0.1,1)):BYEAR +
                 us(trait, init=c(1,0.1,1)):MOTHER,
               residual= ~ id(units):us(trait, init=c(1,0.1,1)),
               data=gryphon,
               maxit=20)
``` 
``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:32 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -4672.301           1.0   1533 09:47:32    0.1
##  2     -4005.615           1.0   1533 09:47:32    0.0
##  3     -3271.483           1.0   1533 09:47:32    0.0 (1 restrained)
##  4     -2761.414           1.0   1533 09:47:32    0.0 (1 restrained)
##  5     -2481.356           1.0   1533 09:47:32    0.0
##  6     -2395.858           1.0   1533 09:47:32    0.0
##  7     -2381.050           1.0   1533 09:47:32    0.0
##  8     -2380.251           1.0   1533 09:47:32    0.0
##  9     -2380.246           1.0   1533 09:47:32    0.0
``` 

Note that we have specified a covariance structure for each random effect and an estimate of the effect of sex on both birth weight and tarsus length by interacting sex with trait in the fixed effect structure.

There will now be twelve (co)variance components reported after running the code:

``` r
summary(modelb)$varcomp
``` 
``` 
##                                             component std.error    z.ratio
## trait:BYEAR!trait_BWT:BWT                   0.9746500 0.2825643  3.4493030
## trait:BYEAR!trait_TARSUS:BWT                0.1624455 0.4185658  0.3881003
## trait:BYEAR!trait_TARSUS:TARSUS             3.7384670 1.2069893  3.0973489
## trait:MOTHER!trait_BWT:BWT                  1.1445030 0.2302062  4.9716424
## trait:MOTHER!trait_TARSUS:BWT              -1.5566953 0.4051695 -3.8420841
## trait:MOTHER!trait_TARSUS:TARSUS            4.8206658 1.3200970  3.6517512
## trait:vm(ANIMAL, ainv)!trait_BWT:BWT        1.9893623 0.4410265  4.5107542
## trait:vm(ANIMAL, ainv)!trait_TARSUS:BWT     3.3168101 0.9030192  3.6730229
## trait:vm(ANIMAL, ainv)!trait_TARSUS:TARSUS 10.2277887 2.8065673  3.6442343
## units:trait!R                               1.0000000        NA         NA
## units:trait!trait_BWT:BWT                   1.8443114 0.3443161  5.3564488
## units:trait!trait_TARSUS:BWT                4.0144042 0.7411868  5.4161844
## units:trait!trait_TARSUS:TARSUS            12.4857174 2.2890289  5.4545915
##                                            bound %ch
## trait:BYEAR!trait_BWT:BWT                      P 0.0
## trait:BYEAR!trait_TARSUS:BWT                   P 0.0
## trait:BYEAR!trait_TARSUS:TARSUS                P 0.0
## trait:MOTHER!trait_BWT:BWT                     P 0.0
## trait:MOTHER!trait_TARSUS:BWT                  P 0.0
## trait:MOTHER!trait_TARSUS:TARSUS               P 0.0
## trait:vm(ANIMAL, ainv)!trait_BWT:BWT           P 0.0
## trait:vm(ANIMAL, ainv)!trait_TARSUS:BWT        P 0.1
## trait:vm(ANIMAL, ainv)!trait_TARSUS:TARSUS     P 0.1
## units:trait!R                                  F 0.0
## units:trait!trait_BWT:BWT                      P 0.0
## units:trait!trait_TARSUS:BWT                   P 0.0
## units:trait!trait_TARSUS:TARSUS                P 0.1
``` 

### 5) Significance testing

Under the model above $r_M$ is estimated as -0.66 and the z.ratio associated with the corresponding covariance ($COV_M$) is >2 (in absolute terms).
We might therefore infer that there is evidence for a strong negative correlation between the traits with respect to the mother and that while maternal identity explains variance in both traits those mothers that tend to produce heavier offspring actually tend to produce offspring with shorter tarsus lengths.

To formally test if $COV_M$ is significantly different from zero, we can compare the log-likelihood for this model:

``` r
modelb$loglik
``` 

``` 
## [1] -2380.246
``` 

to a model in which we specify that $COV_M=0$. Since this constraint reduces the number of parameters to be estimated by one, we can use a likelihood ratio test with one degree of freedom.
To run the constrained model, we modify the G structure definition for the MOTHER random effect to diagonal (`diag`), which means we only estimate the variances (the diagonal of the matrix) but not the covariance:

``` r
modelc <- asreml(fixed = cbind(BWT, TARSUS) ~ trait + trait:SEX,
               random=~ us(trait, init=c(1,0.1,1)):vm(ANIMAL, ainv) +
                 us(trait, init=c(1,0.1,1)):BYEAR +
                 diag(trait, init=c(1,1)):MOTHER,
               residual= ~ id(units):us(trait, init=c(1,0.1,1)),
               data=gryphon,
               maxit=20)
``` 

``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:32 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -4677.820           1.0   1533 09:47:32    0.1
##  2     -4010.442           1.0   1533 09:47:32    0.0
##  3     -3275.409           1.0   1533 09:47:32    0.0
##  4     -2763.519           1.0   1533 09:47:32    0.0
##  5     -2483.732           1.0   1533 09:47:32    0.0
##  6     -2400.242           1.0   1533 09:47:32    0.0
##  7     -2386.663           1.0   1533 09:47:33    0.0
##  8     -2386.049           1.0   1533 09:47:33    0.0
##  9     -2386.045           1.0   1533 09:47:33    0.0
``` 

You can run `summary(modelc)$varcomp` to confirm this worked. We can now obtain the log-likelihood of this model and compare this to that of `modelb` using a likelihood ratio test:

``` r
modelc$loglik
``` 

``` 
## [1] -2386.045
``` 

We can see that the model log-likelihood is now -2386.05. And comparing the models using a likelihood ratio test:

``` r
2*(modelb$loglik-modelc$loglik)
``` 

``` 
## [1] 11.59831
``` 

So our chi-square test statistic is $\chi^2_1= 11.6$. The p-value that goes with this is obtained by:

``` r
1-pchisq(2*(modelb$loglik-modelc$loglik),1)
``` 

``` 
## [1] 0.000660117
``` 

We would therefore conclude that the maternal covariance is significantly different from zero.

We could apply the same procedure to show that the residual (environmental) covariance and the genetic covariance estimates are significantly greater than zero (i.e., heavier individuals tend to have longer tarsus lengths).
In contrast, we should find that the `BYEAR` covariance between the two traits is non-significant.

## Tutorial 3 – A repeated measures animal model

This tutorial will demonstrate how to run a univariate animal model for a trait with repeated observations using the software `ASReml-R` and example data files provided.

### 1) Scenario

Since gryphons are iteroparous, multiple observations of reproductive traits are available for some individuals. Here we have repeated measures of lay date (measured in days after January 1) for individual females varying in age from 2 (age of maturation) up until age 6. Not all females lay every year so the number of observations per female is variable. We want to know how repeatable the trait is, and (assuming it is repeatable) how heritable it is.

### 2) Data files

The pedigree file `gryphonped.txt` is that used in the preceding tutorials but we now use a new data file `gryphonRM.txt`. Columns correspond to individual identity (`ANIMAL`), birth year (`BYEAR`), age in years (`AGE`), year of measurement (`YEAR`) and lay date (`LAYDATE`).
Each row of the data file corresponds to a single phenotypic observation. Here data are sorted by identity and then age so that the repeated observations on individuals are readily apparent. 
However this is not a requirement for analysis - data could equally be sorted by some other variable (e.g., measurement year) or be in a random order.

``` r
gryphonRM <- read.table("/the/path/to/gryphonRM.txt", 
                        header=T, na.strings="NA")

gryphonRM$ANIMAL <- as.factor(gryphonRM$ANIMAL)
gryphonRM$BYEAR <- as.factor(gryphonRM$BYEAR)
gryphonRM$AGE <- as.factor(gryphonRM$AGE)
gryphonRM$YEAR <- as.factor(gryphonRM$YEAR)
gryphonRM$LAYDATE <- as.numeric(gryphonRM$LAYDATE)

head(gryphonRM)
``` 

``` 
##   ANIMAL BYEAR AGE YEAR LAYDATE
## 1      1   990   2  992      19
## 2      1   990   3  993      23
## 3      1   990   4  994      24
## 4      1   990   5  995      23
## 5      1   990   6  996      29
## 6      2   990   2  992      21
``` 

### 3) Estimating repeatability

With repeated measures on individuals it is often of interest, prior to fitting a genetic model, to see how repeatable a trait is.
We can estimate the repeatability of a trait as the proportion of phenotypic variance explained by individual identity.

``` r
modelv <- asreml(fixed = LAYDATE ~ 1, 
                 random =~ ANIMAL,
                 residual=~idv(units),
                 data=gryphonRM, 
                 na.action = na.method(x="omit", y="omit"))
``` 

``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:33 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -10182.83           1.0   1606 09:47:33    0.0
##  2      -8266.10           1.0   1606 09:47:33    0.0
##  3      -6145.01           1.0   1606 09:47:33    0.0
##  4      -4651.57           1.0   1606 09:47:33    0.0
##  5      -3819.31           1.0   1606 09:47:33    0.0
##  6      -3554.22           1.0   1606 09:47:33    0.0
##  7      -3501.56           1.0   1606 09:47:33    0.0
##  8      -3497.58           1.0   1606 09:47:33    0.0
##  9      -3497.54           1.0   1606 09:47:33    0.0
## 10      -3497.54           1.0   1606 09:47:33    0.0
``` 

Note that since we want to estimate the amount of variance explained by individual identity (rather than by additive effects), we fit `ANIMAL` as a normal random effect and we don’t associate it with the pedigree.

This model partitions the phenotypic variance in `LAYDATE` as follows:

``` r
summary(modelv)$varcomp
``` 

``` 
##             component std.error   z.ratio bound %ch
## ANIMAL       11.08634 1.1794319  9.399728     P   0
## units!units  21.29643 0.8896196 23.938798     P   0
## units!R       1.00000        NA        NA     F   0
``` 

Between-individual variance is given by the `ANIMAL` component, while the residual component (`units!units`) represents within-individual variance.
Here then the repeatability of the trait can be determined by hand as 0.34 (i.e., as `11.086/(11.086 + 21.296)`).

Mean lay date might change with age, so we could ask what the repeatability of lay date is after conditioning on age. This would be done by adding `AGE` into the model as a fixed effect.

``` r
modelw <- asreml(fixed = LAYDATE ~ AGE, 
                 random =~ ANIMAL,  
                 residual=~idv(units),
                 data=gryphonRM, 
                 na.action = na.method(x="omit", y="omit"))
``` 

``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:33 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -8402.968           1.0   1602 09:47:33    0.0
##  2     -6912.361           1.0   1602 09:47:33    0.0
##  3     -5274.379           1.0   1602 09:47:33    0.0
##  4     -4143.634           1.0   1602 09:47:33    0.0
##  5     -3541.895           1.0   1602 09:47:33    0.0
##  6     -3372.909           1.0   1602 09:47:33    0.0
##  7     -3347.670           1.0   1602 09:47:33    0.0
##  8     -3346.655           1.0   1602 09:47:33    0.0
##  9     -3346.652           1.0   1602 09:47:33    0.0
``` 

``` r
summary(modelw)$varcomp
``` 

``` 
##             component std.error  z.ratio bound %ch
## ANIMAL       12.28982  1.156115 10.63027     P   0
## units!units  16.37989  0.686619 23.85586     P   0
## units!R       1.00000        NA       NA     F   0
``` 

The repeatability of lay date, after accounting for age effects, is now estimated as 0.43 (i.e., as `12.29/(12.29 + 16.38)`). 
So, just as we saw when estimating $h^2$ in Tutorial 1, the inclusion of fixed effects will alter the estimated effect size if we determine total phenotypic variance as the sum of the variance components.
Thus, proper interpretation is vital.

Here age is modelled as a 5-level factor (specified using the function `as.factor()` at the beginning of the analysis).
We could equally have fitted it as a continuous variable, in which case, given potential for a late life decline, we would probably also include a quadratic term.

### 4) Partitioning additive and permanent environment effects

Generally we expect that the repeatability will set the upper limit for heritability since, while additive genetic effects will cause among-individual variation, so will other types of effects.
Non-additive contributions to fixed among-individual differences are normally referred to as permanent environment effects.
If a trait has repeated measures then it is necessary to model permanent environment effects in an animal model to prevent upward bias in $V_A$.

To illustrate this fit the animal model:

``` r
modelx <- asreml(fixed = LAYDATE ~ AGE, 
                 random =~ vm(ANIMAL, ainv),
                 residual=~idv(units),
                 data=gryphonRM, 
                 na.action = na.method(x="omit", y="omit"))
``` 

``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:33 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -8751.390           1.0   1602 09:47:33    0.0
##  2     -7169.205           1.0   1602 09:47:33    0.0
##  3     -5427.604           1.0   1602 09:47:33    0.0
##  4     -4219.598           1.0   1602 09:47:33    0.0
##  5     -3569.815           1.0   1602 09:47:33    0.0
##  6     -3382.341           1.0   1602 09:47:33    0.0
##  7     -3352.867           1.0   1602 09:47:33    0.0
##  8     -3351.565           1.0   1602 09:47:33    0.0
##  9     -3351.560           1.0   1602 09:47:33    0.0
``` 

Variance components are almost unchanged:

``` r
summary.asreml(modelx)$varcomp
``` 

``` 
##                  component std.error   z.ratio bound %ch
## vm(ANIMAL, ainv)  13.91784  1.443968  9.638607     P   0
## units!units       16.84008  0.707365 23.806768     P   0
## units!R            1.00000        NA        NA     F   0
``` 

This suggests that all of the among-individual variance is – rightly or wrongly – being partitioned as $V_A$ here. To instead obtain an unbiased estimate of $V_A$ we need to allow for both additive genetic and non-genetic sources of individual variation.
We do this by fitting `ANIMAL` twice, once with a pedigree, and once without a pedigree (using `ide()`).

``` r
modely <- asreml(fixed = LAYDATE ~ AGE, 
                 random =~ vm(ANIMAL, ainv) + ide(ANIMAL), 
                 residual=~idv(units),
                 data=gryphonRM, 
                 na.action = na.method(x="omit", y="omit"))
``` 
``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:33 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -7731.394           1.0   1602 09:47:33    0.0
##  2     -6426.548           1.0   1602 09:47:33    0.0
##  3     -4997.252           1.0   1602 09:47:33    0.0
##  4     -4018.486           1.0   1602 09:47:33    0.0
##  5     -3504.988           1.0   1602 09:47:33    0.0
##  6     -3363.160           1.0   1602 09:47:33    0.0
##  7     -3341.611           1.0   1602 09:47:33    0.0
##  8     -3340.682           1.0   1602 09:47:33    0.0
##  9     -3340.679           1.0   1602 09:47:33    0.0
``` 

``` r
summary(modely)$varcomp
``` 

``` 
##                  component std.error   z.ratio bound %ch
## vm(ANIMAL, ainv)  4.876101 1.8087709  2.695809     P   0
## ide(ANIMAL)       7.400983 1.7280113  4.282948     P   0
## units!units      16.380188 0.6866189 23.856300     P   0
## units!R           1.000000        NA        NA     F   0
``` 

The estimate of $V_A$ is now much lower since the additive and permanent environment effects are being properly separated. We can estimate h2h2 and the repeatability from this model:

``` r
vpredict(modely, h2 ~ V1/(V1+V2+V3))
``` 
``` 
##     Estimate         SE
## h2 0.1701523 0.06073974
``` 

``` r
vpredict(modely, repeatability ~ (V1+V2)/(V1+V2+V3))
``` 

``` 
##                Estimate         SE
## repeatability 0.4284108 0.02741602
``` 

### 5) Adding additional effects and testing significance

Models of repeated measures can be extended to include other fixed or random effects. For example try including year of measurement (`YEAR`) and birth year (`BYEAR`) as random effects.

``` r
modelz <- asreml(fixed = LAYDATE ~ AGE,
                 random =~ vm(ANIMAL, ainv) + ide(ANIMAL) +
                 YEAR + BYEAR,
                 residual=~idv(units),
                 data=gryphonRM, 
                 na.action = na.method(x="omit", y="omit"))
``` 

``` 
## Model fitted using the sigma parameterization.
## ASReml 4.1.0 Wed Sep 16 09:47:33 2020
##           LogLik        Sigma2     DF     wall    cpu
##  1     -4650.748           1.0   1602 09:47:33    0.0
##  2     -4088.264           1.0   1602 09:47:33    0.0
##  3     -3494.147           1.0   1602 09:47:33    0.0
##  4     -3127.161           1.0   1602 09:47:33    0.0 (1 restrained)
##  5     -2976.449           1.0   1602 09:47:33    0.0 (1 restrained)
##  6     -2955.785           1.0   1602 09:47:33    0.0 (1 restrained)
##  7     -2955.097           1.0   1602 09:47:33    0.0 (1 restrained)
##  8     -2955.095           1.0   1602 09:47:33    0.0 (1 restrained)
##  9     -2955.095           1.0   1602 09:47:33    0.0
``` 

``` r
summary(modelz)$varcomp
``` 

``` 
##                     component std.error   z.ratio bound %ch
## BYEAR            1.650876e-07        NA        NA     B   0
## YEAR             7.938576e+00 1.9344619  4.103765     P   0
## vm(ANIMAL, ainv) 4.815136e+00 1.6682351  2.886365     P   0
## ide(ANIMAL)      8.433325e+00 1.5495778  5.442337     P   0
## units!units      7.795560e+00 0.3324411 23.449443     P   0
## units!R          1.000000e+00        NA        NA     F   0
``` 

This model will return additional variance components corresponding to variation in lay dates between years of measurement and between birth cohorts of females.
$V_'BYEAR} is very low and if you compare this model to a reduced model with BYEAR excluded the log-likelihood remains unchanged.

`YEAR` effects could alternatively be included as fixed effects (try this!).
This will reduce $V_R$ and increase the estimates of heritability and repeatability, which must now be interpreted as proportions of phenotypic variance after conditioning on both age and year of measurement effects.