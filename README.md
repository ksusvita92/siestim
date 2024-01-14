# siestim


Welcome to the `siestim` project! The project aims to estimate the serial interval distribution of an infectious disease even when some portion of the cases are undetected. Under partially sampled data, purported infector-infectee pairs may actually be separated by one or more unsampled cases in between. Misunderstanding such pairs as direct transmissions will result in overestimating the length of serial intervals. On the other hand, two cases that are infected by an unseen third case (known as coprimary transmission) may be classified as a direct transmission pair, leading to an underestimation of the serial interval. Here, we introduce a method to jointly estimate the distribution of serial intervals factoring in these two sources of error. We simultaneously estimate the distribution of the number of intermediate cases between purported infector/infectee pairs, as well as the fraction of such cases that are coprimary. 


This page is dedicated as a tutorial and complimentary of paper 
**A Method to Estimate Serial Interval Distribution**, by
Kurnia Susvitasari, Paul Tupper, Jessica Stockdale, Caroline Colijn


## Installing the *siestim* package
To install the package, you will need the `remotes` package:

```r
install.packages("remotes")
remotes::install_github("ksusvita92/siestim")
```


Once installed, the package can be loaded using:

```r
library("siestim")
```

## Tutorial
The package tutorial can be found [here](https://github.com/ksusvita92/siestim/blob/main/tutorial/tutorial_html.Rmd). Other version (html or pdf) can be found under subfolder **tutorial**.

## Contact
- For general questions and bug reports, please contact <kurniasusvitasari@ui.ac.id>
