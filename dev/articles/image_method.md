# The original image-level test in spicyR

Abstract

The original spicyR test, available as spicy(method = “image”):
per-image L-function summaries compared between conditions with weighted
linear and mixed-effects models. Since spicyR 2.0 the default is spicyR
Cell; see the main vignette.

## Installation

``` r

if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install("spicyR", version = "devel")   # spicyR 2.0, until the next Bioconductor release
```

``` r

# load required packages
library(SummarizedExperiment)


library(spicyR)
library(ggplot2)
library(SpatialExperiment)
library(SpatialDatasets)
library(imcRtools)
library(dplyr)
library(survival)
```

## Overview

This guide provides step-by-step instructions on how to apply a linear
model to multiple segmented and labelled images to assess how the
localisation of different cell types changes across different disease
conditions.

## Example data

We use the Keren et al. (2018) breast cancer dataset to compare the
spatial distribution of immune cells in individuals with different
levels of tumour infiltration (cold and compartmentalised).

The data is stored as a `SpatialExperiment` object and contains
single-cell spatial data from 41 images.

``` r

kerenSPE <- SpatialDatasets::spe_Keren_2018()
```

The cell types in this dataset includes 11 immune cell types (double
negative CD3 T cells, CD4 T cells, B cells, monocytes, macrophages, CD8
T cells, neutrophils, natural killer cells, dendritic cells, regulatory
T cells), 2 structural cell types (endothelial, mesenchymal), 2 tumour
cell types (keratin+ tumour, tumour) and one unidentified category.

## Linear modelling

To investigate changes in localisation between two different cell types,
we measure the level of localisation between two cell types by modelling
with the L-function. The L-function is a variance-stabilised K-function
given by the equation

``` math
\widehat{L_{ij}} (r) = \sqrt{\frac{\widehat{K_{ij}}(r)}{\pi}}
```

with $`\widehat{K_{ij}}`$ defined as

``` math
\widehat{K_{ij}} (r) = \frac{|W|}{n_i n_j} \sum_{n_i} \sum_{n_j} 1 \{d_{ij} \leq r \} e_{ij} (r)
```

where $`\widehat{K_{ij}}`$ summarises the degree of co-localisation of
cell type $`j`$ with cell type $`i`$, $`n_i`$ and $`n_j`$ are the number
of cells of type $`i`$ and $`j`$, $`|W|`$ is the image area, $`d_{ij}`$
is the distance between two cells and $`e_{ij} (r)`$ is an edge
correcting factor.

Specifically, the mean difference between the experimental function and
the theoretical function is used as a measure for the level of
localisation, defined as

``` math
u = \sum_{r' = r_{\text{min}}}^{r_{\text{max}}} \widehat L_{ij, \text{Experimental}} (r') - \widehat L_{ij, \text{Poisson}} (r')
```

where $`u`$ is the sum is taken over a discrete range of $`r`$ between
$`r_{\text{min}}`$ and $`r_{\text{max}}`$. Differences of the statistic
$`u`$ between two conditions is modelled using a weighted linear model.

### Test for change in localisation for a specific pair of cells

Firstly, we can test whether one cell type tends to be more localised
with another cell type in one condition compared to the other. This can
be done using the
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md)
function, where we specify the `condition` parameter.

In this example, we want to see whether or not neutrophils (`to`) tend
to be found around CD8 T cells (`from`) in compartmentalised tumours
compared to cold tumours. Given that there are 3 conditions, we can
specify the desired conditions by setting the order of our `condition`
factor. `spicy` will choose the first level of the factor as the base
condition and the second level as the comparison condition. `spicy` will
also naturally coerce the `condition` column into a factor if it is not
already a factor. The column containing cell type annotations can be
specified using the `cellType` argument. By default, `spicy` uses the
column named `cellType` in the `SpatialExperiment` object.

``` r

spicyTestPair <- spicy(method = "image", 
  kerenSPE,
  condition = "tumour_type",
  from = "CD8_T_cell",
  to = "Neutrophils"
)

topPairs(spicyTestPair)
#>                         intercept coefficient      p.value   adj.pvalue
#> CD8_T_cell__Neutrophils  -109.081    112.0185 2.166646e-05 2.166646e-05
#>                               from          to
#> CD8_T_cell__Neutrophils CD8_T_cell Neutrophils
```

We obtain a `spicy` object which details the results of the modelling
performed. The
[`topPairs()`](https://sydneybiox.github.io/spicyR/dev/reference/topPairs.md)
function can be used to obtain the associated coefficients and p-value.

As the `coefficient` in `spicyTestPair` is positive, we find that
neutrophils are significantly more likely to be found near CD8 T cells
in the compartmentalised tumours group compared to the cold tumour
group.

### Test for change in localisation for all pairwise cell combinations

We can perform what we did above for all pairwise combinations of cell
types by excluding the `from` and `to` parameters in
[`spicy()`](https://sydneybiox.github.io/spicyR/dev/reference/spicy.md).

``` r

spicyTest <- spicy(method = "image", 
  kerenSPE,
  condition = "tumour_type"
)

topPairs(spicyTest)
#>                             intercept coefficient      p.value   adj.pvalue
#> Macrophages__dn_T_CD3       56.446065   -50.08474 1.080268e-07 3.035554e-05
#> dn_T_CD3__Macrophages       54.987150   -48.38664 2.194026e-07 3.082607e-05
#> Macrophages__DC_or_Mono     73.239408   -59.90362 5.224650e-06 4.893755e-04
#> DC_or_Mono__Macrophages     71.777083   -58.46833 7.431188e-06 5.220409e-04
#> dn_T_CD3__dn_T_CD3         -63.786032   100.61010 2.878802e-05 1.208706e-03
#> Neutrophils__dn_T_CD3      -63.141839    69.64356 2.891869e-05 1.208706e-03
#> dn_T_CD3__Neutrophils      -63.133727    70.15508 3.011011e-05 1.208706e-03
#> DC__Macrophages             96.893239   -92.55112 1.801305e-04 5.758112e-03
#> Macrophages__DC             96.896215   -93.25194 1.844235e-04 5.758112e-03
#> CD4_T_cell__Keratin_Tumour  -4.845036   -22.14995 2.834660e-04 7.409012e-03
#>                                   from             to
#> Macrophages__dn_T_CD3      Macrophages       dn_T_CD3
#> dn_T_CD3__Macrophages         dn_T_CD3    Macrophages
#> Macrophages__DC_or_Mono    Macrophages     DC_or_Mono
#> DC_or_Mono__Macrophages     DC_or_Mono    Macrophages
#> dn_T_CD3__dn_T_CD3            dn_T_CD3       dn_T_CD3
#> Neutrophils__dn_T_CD3      Neutrophils       dn_T_CD3
#> dn_T_CD3__Neutrophils         dn_T_CD3    Neutrophils
#> DC__Macrophages                     DC    Macrophages
#> Macrophages__DC            Macrophages             DC
#> CD4_T_cell__Keratin_Tumour  CD4_T_cell Keratin_Tumour
```

Again, we obtain a `spicy` object which outlines the result of the
linear models performed for each pairwise combination of cell types.

We can also examine the L-function metrics of individual images by using
the convenient
[`bind()`](https://sydneybiox.github.io/spicyR/dev/reference/bind.md)
function on our `spicyTest` results object.

``` r

bind(spicyTest)[1:5, 1:5]
#>   imageID         condition Keratin_Tumour__Keratin_Tumour
#> 1       1             mixed                      -2.300602
#> 2       2             mixed                      -1.989699
#> 3       3 compartmentalised                      11.373530
#> 4       4 compartmentalised                      33.931133
#> 5       5 compartmentalised                      17.922818
#>   dn_T_CD3__Keratin_Tumour B_cell__Keratin_Tumour
#> 1                -5.298543             -20.827279
#> 2               -16.020022               3.025815
#> 3               -21.857447             -24.962913
#> 4               -36.438476             -40.470221
#> 5               -20.816783             -38.138076
```

The results can be represented as a bubble plot using the
[`signifPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/signifPlot.md)
function.

``` r

signifPlot(
  spicyTest,
  breaks = c(-3, 3, 1),
  marksToPlot = c("Macrophages", "DC_or_Mono", "dn_T_CD3", "Neutrophils",
                  "CD8_T_cell", "Keratin_Tumour")
)
```

![](image_method_files/figure-html/chunk-08-1.png)

Here, we can observe that the most significant relationships occur
between macrophages and double negative CD3 T cells, suggesting that the
two cell types are far more dispersed in compartmentalised tumours
compared to cold tumours.

To examine a specific cell type-cell type relationship in more detail,
we can use
[`spicyBoxPlot()`](https://sydneybiox.github.io/spicyR/dev/reference/spicyBoxPlot.md)
and specify either `from = "Macrophages"` and `to = "dn_T_CD3"` or
`rank = 1`.

``` r

spicyBoxPlot(results = spicyTest, rank = 1)
```

![](image_method_files/figure-html/chunk-09-1.png)

## Linear modelling for custom metrics

`spicyR` can also be applied to custom distance or abundance metrics. A
kNN interactions graph can be generated with the function
`buildSpatialGraph` from the `imcRtools` package. This generates a
`colPairs` object inside of the `SpatialExperiment` object.

`spicyR` provides the function `convPairs` for converting a `colPairs`
object into an abundance matrix by calculating the average number of
nearby cells types for every cell type for a given `k`. For example, if
there exists on average 5 neutrophils for every macrophage in image 1,
the column `Neutrophil__Macrophage` would have a value of 5 for image 1.

``` r

kerenSPE <- imcRtools::buildSpatialGraph(kerenSPE, 
                                         img_id = "imageID", 
                                         type = "knn", k = 20,
                                        coords = c("x", "y"))

pairAbundances <- convPairs(kerenSPE,
                  colPair = "knn_interaction_graph")

head(pairAbundances["B_cell__B_cell"])
#>    B_cell__B_cell
#> 1      12.7349608
#> 10      0.2777778
#> 11      0.0000000
#> 12      1.3333333
#> 13      1.2200957
#> 14      0.0000000
```

The custom distance or abundance metrics can then be included in the
analysis with the `alternateResult` parameter. The `Statial` package
contains other custom distance metrics which can be used with `spicy`.

``` r

spicyTestColPairs <- spicy(method = "image", 
  kerenSPE,
  condition = "tumour_type",
  alternateResult = pairAbundances,
  weights = FALSE
)

topPairs(spicyTestColPairs)
#>                            intercept coefficient     p.value adj.pvalue
#> CD8_T_cell__Neutrophils  0.833333333  -0.7592968 0.002645466  0.3291833
#> B_cell__Tumour           0.001937984   0.2602822 0.004872664  0.3291833
#> Other_Immune__NK         0.012698413   0.2612881 0.005673068  0.3291833
#> Unidentified__CD8_T_cell 0.106626794   0.6387339 0.005906526  0.3291833
#> dn_T_CD3__NK             0.004242424   0.2148797 0.006317829  0.3291833
#> CD4_T_cell__Neutrophils  0.036213602   0.2947696 0.007902670  0.3291833
#> Tregs__CD4_T_cell        0.128876212   0.5726201 0.010207087  0.3291833
#> Endothelial__DC          0.008771930   0.3008523 0.011189533  0.3291833
#> Tumour__Neutrophils      0.021638939   0.2529045 0.011388850  0.3291833
#> Mesenchymal__Neutrophils 0.004504505   0.2494301 0.012761315  0.3291833
#>                                  from          to
#> CD8_T_cell__Neutrophils    CD8_T_cell Neutrophils
#> B_cell__Tumour                 B_cell      Tumour
#> Other_Immune__NK         Other_Immune          NK
#> Unidentified__CD8_T_cell Unidentified  CD8_T_cell
#> dn_T_CD3__NK                 dn_T_CD3          NK
#> CD4_T_cell__Neutrophils    CD4_T_cell Neutrophils
#> Tregs__CD4_T_cell               Tregs  CD4_T_cell
#> Endothelial__DC           Endothelial          DC
#> Tumour__Neutrophils            Tumour Neutrophils
#> Mesenchymal__Neutrophils  Mesenchymal Neutrophils
```

``` r

signifPlot(
  spicyTestColPairs,
  breaks = c(-3, 3, 1),
  marksToPlot = c("Macrophages", "dn_T_CD3", "CD4_T_cell", 
                  "B_cell", "DC_or_Mono", "Neutrophils", "CD8_T_cell")
)
```

![](image_method_files/figure-html/chunk-12-1.png)

## Performing survival analysis

`spicy` can also be used to perform survival analysis to asses whether
changes in co-localisation between cell types are associated with
survival probability. `spicy` requires the `SingleCellExperiment` object
being used to contain a column called `survival` as a `Surv` object.

``` r

kerenSPE$event = 1 - kerenSPE$Censored
kerenSPE$survival = Surv(kerenSPE$`Survival_days_capped*`, kerenSPE$event)
```

We can then perform survival analysis using the `spicy` function by
specifying `condition = "survival"`. We can then access the
corresponding coefficients and p-values by accessing the
`survivalResults` slot in the `spicy` results object.

``` r

# Running survival analysis
spicySurvival = spicy(method = "image", kerenSPE,
                      condition = "survival")

# top 10 significant pairs
head(spicySurvival$survivalResults, 10)
#>                      test         coef     se.coef      p.value
#> 1     Other_Immune__Tregs  0.023565206 0.008656342 8.929522e-06
#> 2       CD4_T_cell__Tregs  0.017697194 0.006849369 1.241121e-05
#> 3     Tregs__Other_Immune  0.023714941 0.008733582 1.264807e-05
#> 4       Tregs__CD4_T_cell  0.017081289 0.006758900 2.852449e-05
#> 5  CD8_T_cell__CD8_T_cell  0.006050072 0.002723387 3.316528e-04
#> 6      Tumour__CD8_T_cell -0.030532595 0.011429924 6.168710e-04
#> 7      CD8_T_cell__Tumour -0.030478957 0.011593521 7.214929e-04
#> 8    CD4_T_cell__dn_T_CD3  0.008453662 0.003533034 7.936924e-04
#> 9    dn_T_CD3__CD4_T_cell  0.008398973 0.003530625 9.371907e-04
#> 10       DC__Other_Immune -0.028885515 0.012294396 1.034119e-03
```

## Accounting for tissue inhomogeneity

The `spicy` function can also account for tissue inhomogeneity to avoid
false positives or negatives. This can be done by setting the `sigma =`
parameter within the spicy function. By default, `sigma` is set to
`NULL`, and `spicy` assumes a homogeneous tissue structure.

For example, when we examine the L-function for
`Keratin_Tumour__Neutrophils` when `sigma = NULL` and `Rs = 100`, the
value is positive, indicating attraction between the two cell types.

``` r

# filter SPE object to obtain image 24 data
kerenSubset = kerenSPE[, colData(kerenSPE)$imageID == "24"]

pairwiseAssoc = getPairwise(kerenSubset, 
                            sigma = NULL, 
                            Rs = 100) |>
  as.data.frame()

pairwiseAssoc[["Keratin_Tumour__Neutrophils"]]
#> [1] 10.88892
```

When we specify `sigma = 20` and re-calculate the L-function, it
indicates that there is no relationship between `Keratin_Tumour` and
`Neutrophils`, i.e., there is no major attraction or dispersion, as it
now takes into account tissue inhomogeneity.

``` r

pairwiseAssoc = getPairwise(kerenSubset, 
                            sigma = 20, 
                            Rs = 100) |>
  as.data.frame()

pairwiseAssoc[["Keratin_Tumour__Neutrophils"]]
#> [1] 0.9024836
```

We can use the `plotImage` function to plot any pair of cell types for a
specific image and visually inspect cellular relationships.

``` r

plotImage(kerenSPE, "24", from = "Keratin_Tumour", to = "Neutrophils")
```

![](image_method_files/figure-html/chunk-17-1.png)

Plotting image 24 shows that the supposed co-localisation occurs due to
the dense cluster of cells near the bottom of the image.

## Mixed effects modelling

`spicyR` supports mixed effects modelling when multiple images are
obtained for each subject. In this case, `subject` is treated as a
random effect and `condition` is treated as a fixed effect. To perform
mixed effects modelling, we can specify the `subject` parameter in the
`spicy` function.

``` r

data("diabetesData")
spicyMixedTest <- spicy(method = "image",
  diabetesData,
  condition = "stage",
  subject = "case"
)
topPairs(spicyMixedTest)
#>                          intercept coefficient      p.value adj.pvalue
#> beta__delta           6.081500e+01  -15.433250 0.0006616383 0.09171141
#> delta__beta           6.090722e+01  -15.276363 0.0007164954 0.09171141
#> B__Th                -2.664535e-15   10.480122 0.0127535357 0.42316682
#> delta__delta          7.021913e+01  -16.358091 0.0155005353 0.42316682
#> Th__B                 1.776357e-15   10.003453 0.0173028300 0.42316682
#> B__unknown           -8.881784e-16    4.584273 0.0182158777 0.42316682
#> otherimmune__naiveTc -3.179087e+00   11.944639 0.0199646991 0.42316682
#> unknown__macrophage   4.339586e+00   -5.274674 0.0222619404 0.42316682
#> unknown__B            1.165734e-15    4.680750 0.0244004344 0.42316682
#> macrophage__unknown   4.305429e+00   -4.886438 0.0249693485 0.42316682
#>                             from         to
#> beta__delta                 beta      delta
#> delta__beta                delta       beta
#> B__Th                          B         Th
#> delta__delta               delta      delta
#> Th__B                         Th          B
#> B__unknown                     B    unknown
#> otherimmune__naiveTc otherimmune    naiveTc
#> unknown__macrophage      unknown macrophage
#> unknown__B               unknown          B
#> macrophage__unknown   macrophage    unknown
```

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] survival_3.8-6              dplyr_1.2.1                
#>  [3] imcRtools_1.18.1            SpatialDatasets_1.10.0     
#>  [5] ExperimentHub_3.2.2         AnnotationHub_4.2.2        
#>  [7] BiocFileCache_3.2.0         dbplyr_2.6.0               
#>  [9] SpatialExperiment_1.22.0    SingleCellExperiment_1.34.0
#> [11] ggplot2_4.0.3               spicyR_1.99.4              
#> [13] SummarizedExperiment_1.42.0 Biobase_2.72.0             
#> [15] GenomicRanges_1.64.0        Seqinfo_1.2.0              
#> [17] IRanges_2.46.0              S4Vectors_0.50.3           
#> [19] BiocGenerics_0.58.1         generics_0.1.4             
#> [21] MatrixGenerics_1.24.0       matrixStats_1.5.0          
#> [23] BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.6.1          later_1.4.8            bitops_1.1-0          
#>   [4] filelock_1.0.3         tibble_3.3.1           svgPanZoom_0.3.4      
#>   [7] polyclip_1.10-7        lifecycle_1.0.5        httr2_1.3.0           
#>  [10] sf_1.1-3               vroom_1.7.1            lattice_0.22-9        
#>  [13] MASS_7.3-65            magrittr_2.0.5         sass_0.4.10           
#>  [16] rmarkdown_2.32         jquerylib_0.1.4        yaml_2.3.12           
#>  [19] httpuv_1.6.17          otel_0.2.0             spatstat.sparse_3.2-0 
#>  [22] sp_2.2-3               DBI_1.3.0              RColorBrewer_1.1-3    
#>  [25] abind_1.4-8            purrr_1.2.2            ggraph_2.2.2          
#>  [28] RCurl_1.98-1.20        tweenr_2.0.3           rappdirs_0.3.4        
#>  [31] ggrepel_0.9.8          RTriangle_1.6-0.15     spatstat.utils_3.2-5  
#>  [34] terra_1.9-50           pheatmap_1.0.13        units_1.0-1           
#>  [37] goftest_1.2-3          spatstat.random_3.5-2  pkgdown_2.2.1         
#>  [40] svglite_2.2.2          codetools_0.2-20       DelayedArray_0.38.2   
#>  [43] DT_0.34.0              ggforce_0.5.0          tidyselect_1.2.1      
#>  [46] raster_3.6-32          farver_2.1.2           viridis_0.6.5         
#>  [49] spatstat.explore_3.8-3 jsonlite_2.0.0         BiocNeighbors_2.6.0   
#>  [52] e1071_1.7-17           tidygraph_1.3.1        systemfonts_1.3.2     
#>  [55] tools_4.6.1            ragg_1.5.2             Rcpp_1.1.2            
#>  [58] glue_1.8.1             gridExtra_2.3.1        SparseArray_1.12.3    
#>  [61] xfun_0.61              EBImage_4.54.0         HDF5Array_1.40.0      
#>  [64] shinydashboard_0.7.3   withr_3.0.3            BiocManager_1.30.27   
#>  [67] fastmap_1.2.0          rhdf5filters_1.24.1    digest_0.6.39         
#>  [70] R6_2.6.1               mime_0.13              textshaping_1.0.5     
#>  [73] tensor_1.5.1           spatstat.data_3.1-9    jpeg_0.1-11           
#>  [76] RSQLite_3.53.3         h5mread_1.4.1          tidyr_1.3.2           
#>  [79] data.table_1.18.6.1    class_7.3-23           graphlayouts_1.2.5    
#>  [82] httr_1.4.9             htmlwidgets_1.6.4      S4Arrays_1.12.1       
#>  [85] pkgconfig_2.0.3        gtable_0.3.6           blob_1.3.0            
#>  [88] S7_0.2.2               XVector_0.52.0         htmltools_0.5.9       
#>  [91] bookdown_0.48          fftwtools_0.9-11       scales_1.4.0          
#>  [94] png_0.1-9              spatstat.univar_3.2-0  knitr_1.52            
#>  [97] tzdb_0.5.0             rjson_0.2.23           nlme_3.1-169          
#> [100] curl_8.0.0             proxy_0.4-29           cachem_1.1.0          
#> [103] rhdf5_2.56.1           stringr_1.6.0          KernSmooth_2.23-26    
#> [106] BiocVersion_3.23.1     parallel_4.6.1         vipor_0.4.7           
#> [109] concaveman_1.2.0       AnnotationDbi_1.74.0   desc_1.4.3            
#> [112] pillar_1.11.1          grid_4.6.1             vctrs_0.7.3           
#> [115] promises_1.5.0         distances_0.1.13       beachmat_2.28.0       
#> [118] xtable_1.8-8           beeswarm_0.4.0         evaluate_1.0.5        
#> [121] readr_2.2.0            magick_2.9.1           cli_3.6.6             
#> [124] locfit_1.5-9.12        compiler_4.6.1         rlang_1.3.0           
#> [127] crayon_1.5.3           labeling_0.4.3         classInt_0.4-11       
#> [130] fs_2.1.0               ggbeeswarm_0.7.3       stringi_1.8.9         
#> [133] deldir_2.0-4           viridisLite_0.4.3      BiocParallel_1.46.0   
#> [136] nnls_1.6               cytomapper_1.24.0      Biostrings_2.80.2     
#> [139] tiff_0.1-12            spatstat.geom_3.8-3    scrapper_1.6.3        
#> [142] Matrix_1.7-5           hms_1.1.4              bit64_4.8.6           
#> [145] Rhdf5lib_2.0.0         KEGGREST_1.52.2        shiny_1.14.0          
#> [148] igraph_2.3.4           memoise_2.0.1          bslib_0.12.0          
#> [151] bit_4.6.0
```

## References

Keren, L, M Bosse, D Marquez, et al. 2018. “A Structured Tumor-Immune
Microenvironment in Triple Negative Breast Cancer Revealed by
Multiplexed Ion Beam Imaging.” *Cell* 174 (6): 1373–1387.e19.
