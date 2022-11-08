# spicyR

Spatial analysis of in situ cytometry data.

## Overview


**spicyR** provides a series of functions to aid in the analysis of both 
    immunofluorescence and mass cytometry imaging data as well as other assays that 
    can deeply phenotype individual cells and their spatial location such as 
    high-definiton spatial transcriptomics. 

## Installation

For the Bioconductor release version.
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("spicyR")
```

If you would like the most up-to-date features, install the most recent development version.
```r
# Install the development version from Bioconductor:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
# This will update all your Bioconductor packages to devel version
BiocManager::install(version='devel')

BiocManager::install("spicyR")

# Otherwise install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("ellispatrick/spicyR")
library(spicyR)
```

## Submitting an issue or feature request

`spicyR` is still under active development. We would greatly appreciate any and 
all feedback related to the package.

* R package related issues should be raised [here](https://github.com/ellispatrick/spicyR/issues).
* For general questions and feedback, please contact us directly via [ellis.patrick@sydney.edu.au](mailto:ellis.patrick@sydney.edu.au).


## Author

* **Nicolas Canete**
* **Ellis Patrick**  - [@TheEllisPatrick](https://twitter.com/TheEllisPatrick)

