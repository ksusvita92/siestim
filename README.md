# siestim
Serial Interval Estimation

Welcome to the introduction of `siestim` project! The project aims to estimate serial interval distribution when some portian of data are sampled. Under partially sampled data, purported infector-infectee pairs may actually be separated by one or more unsampled cases in between. Misunderstanding such pairs as direct transmissions will result in overestimating the length of serial intervals. On the other hand, two cases that are infected by an unseen third case may be classified as a direct transmission pair (known as coprimary transmission), leading to an underestimation of the serial interval. Here, we introduce a method to jointly estimate the distribution of serial intervals factoring in these two sources of error. We simultaneously estimate the distribution of the number of intermediate cases between purported infector/infectee pairs, as well as the fraction of such cases that are coprimary. 


This page is dedicated as a tutorial and complimentary of paper 
**Serial Interval Estimation**, by
Kurnia Susvitasari, Paul Tupper, Caroline Colijn


## Installing *siestim* package
To install the package, you will need `remotes` package:

```r
install.packages("remotes")
remotes::install_github("ksusvita92/siestim")
```


Once installed, the package can be loaded using:

```r
library("siestim")
```

## Tutorial and asking questions
- For the general question and bugs report, send to <ksusvita@gmail.com>
