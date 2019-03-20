# ImSpectR
R package to quantify spectratypes for immune repertoire diversity

Version: 1.5.1<br>
Author and Maintainer: [Martijn Cordes](mailto:m.cordes@lumc.nl) <br>
Description: With this package .fsa intensity files or files with sequencing lengths representing the CDR3 region of B and T-cell receptors can be loaded, visualized and scored for quantification of the immune repertoire. 

## Installation

To install directly from GitHub, run this:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("martijn-cordes/ImSpectR")
```

## Analyze Spectratype data

To analyze a preloaded test dataset of 22 Vb families of a single C57Bl/6 mouse please run:

#### Preprocessing

```r
library(ImSpectR)
data(fsa.loaded)
my.ladder <- c(35, 50, 75, 100, 139, 150, 160, 200, 250,300, 340, 350, 400, 450, 490, 500)
ladder.info.attach(stored=fsa.loaded, ladder=my.ladder, ladd.init.thresh=1000, draw=F,method="iter2")
spectratype_dataset <- get_basepair_positions(fsa.loaded, cols = 1, my.ladder, channel.ladder=NULL,  init.thresh=1750, ladd.init.thresh=1000)
```
#### Score a single sample 

```r
score_sample(spectratype_dataset[1], no.peaks = 7)
```
#### Score complete dataset

```r
score_dataset(spectratype_dataset, no.peaks = 7)
```

## Analyze CDR3 sequencing data

To load 

#### Preprocessing

```r
library(ImSpectR)
data(fsa.loaded)
my.ladder <- c(35, 50, 75, 100, 139, 150, 160, 200, 250,300, 340, 350, 400, 450, 490, 500)
ladder.info.attach(stored=fsa.loaded, ladder=my.ladder, ladd.init.thresh=1000, draw=F,method="iter2")
spectratype_dataset <- get_basepair_positions(fsa.loaded, cols = 1, my.ladder, channel.ladder=NULL,  init.thresh=1750, ladd.init.thresh=1000)
```
#### Score a single sample 

```r
score_sample(spectratype_dataset[1], no.peaks = 7)
```
#### Score complete dataset

```r
score_dataset(spectratype_dataset, no.peaks = 7)
```



