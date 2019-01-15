# ImSpectR
R package to quantify spectratypes for immune repertoire diversity

Version: 1.0.0

Author: Martijn Cordes
Maintainer: Martijn Cordes <m.cordesf@lumc.nl>
Description: With this package .fsa files containing spectratype data representing the CDR3 region of B and T-cell receptors can be loaded, visualized and scored for quantification of the immune repertoire. 

## Installation

To install directly from GitHub, run this:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("rstudio/rsconnect")
```

## Usage

To run a preloaded test dataset of 22 Vb families of a B6 mouse run:

#### Preprocessing

```r
data(fsa.loaded)
ladder.info.attach(stored=fsa.loaded, ladder=my.ladder, ladd.init.thresh=1000, draw=F,method="iter2")
all_samples <- get_basepair_positions(fsa.loaded, cols = 1, my.ladder, channel.ladder=NULL,  init.thresh=1750, ladd.init.thresh=1000)
```
#### Score a single sample 

```r
score_sample(all_samples[1], no.peaks = 7)
```
#### Score complete dataset

```r
score_dataset(all_samples, no.peaks = 7)
```