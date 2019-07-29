# ImSpectR
R package to quantify CDR3 data for immune repertoire diversity

Version: 1.7.5<br>
Author and Maintainer: [Martijn Cordes](mailto:m.cordes@lumc.nl) <br>
Description: With this package .fsa intensity files or files with sequencing lengths representing the CDR3 region of T-cell receptors can be loaded, visualized and scored for quantification of the immune repertoire. 

## Installation

To install directly from GitHub, run this:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("martijn-cordes/ImSpectR")
```
If the devtools-based approach does not work, you can download one of the built tar-balls from the builds directory and manually install the package from source by executing the following lines (replace DOWNLOAD_FOLDER with the absolute path to the tar-ball and VERSION with the version number of the downloaded package):

```r
install.packages("/DOWNLOAD_FOLDER/ImSpectR_VERSION.tar.gz",
                     repos = NULL,
                     type  = "source")

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

#### Loading .fsa files and analyze dataset

The 22 raw fsa files from the preloaded dataset are also included as external data as an example how to start analyzing Spectratype data from scratch:

```r
#Load raw data
folder <- system.file("extdata/spectratype_data/", package="ImSpectR")
fsa.loaded <- storing.inds(folder)

#Preprocess
my.ladder <- c(35, 50, 75, 100, 139, 150, 160, 200, 250,300, 340, 350, 400, 450, 490, 500)
ladder.info.attach(stored=fsa.loaded, ladder=my.ladder, ladd.init.thresh=1000, draw=F,method="iter2")
spectratype_dataset <- get_basepair_positions(fsa.loaded, cols = 1, my.ladder, channel.ladder=NULL,  init.thresh=1750, ladd.init.thresh=1000)

#Score dataset
score_dataset(spectratype_dataset, no.peaks = 7)
```


## Analyze CDR3 sequencing data

A control data set of healthy mice was included as a test set The data set is comprised of 2 C57Bl/6 mice from Charles River Laboratories. Samples include spleen tissue from each mouse. Tab seperated (tsv) files were dowloaded from the Adaptive Biotech database for immunoSEQ data: https://clients.adaptivebiotech.com/pub/b4ac7a84-1e69-4d60-8254-845720454d7d

#### Preprocessing

```r
library(ImSpectR)
folder <- system.file("extdata/cdr_sequencing_data/", package="ImSpectR")
sequencing_dataset <- load_cdr3_seq(folder, cdr3Length_column = 5, geneFamily_column = 7, sep="")
```
#### Score a single sample 

```r
score_sample(sequencing_dataset[1], no.peaks = 7, peak.window = c(0,200))
```
#### Score complete dataset

```r
score_dataset(sequencing_dataset, no.peaks = 7, peak.window = c(0,200))
```



