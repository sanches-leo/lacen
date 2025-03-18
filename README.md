# lacen

## Introduction

`lacen` is an R package designed for analyzing gene expression data, with a focus on long non-coding RNAs (lncRNAs) and weighted gene co-expression network analysis (WGCNA).
It facilitates data handling, transformation and network construction, ensuring robust stability assessment through bootstrap analysis.

## Dependencies

    WGCNA,
    ggplot2,
    gprofiler2,
    fastcluster,
    rrvgo,
    limma,
    rtracklayer,
    foreach,
    doParallel,
    Polychrome,
    ggraph,
    graphlayouts,
    igraph,
    scatterpie,
    org.Hs.eg.db,
    methods,
    utils,
    stats,
    parallel


## Installation

`lacen` is not currently available in CRAN. To ensure all Bioconductor dependencies are installed, use BiocManager:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sanches-leo/lacen")
```

## Loading the Data

As input, `lacen`requires five different files: 
    - The count dataframe, composed by the RNA-seq raw counts with the genes in the rows and samples in the columns. The gene name/ID should be in the row names and the sample names or ID in the column names.
    - The traits dataframe, indicating the which class does the sample is.  A two column dataframe, with the sample ID in the "Sample" column and the condition code on the "Trait" column. The conditions must be numeric codified, with each number indicating a different condition.
    - A differential expression dataframe. The user can use the standard Limma differential expression output, or give a dataframe with the columns "ID", "log2FC" and "pvalue", in this order.
    - The gene annotation data, a two columns dataframe with sequence IDs on "gene_id" column and Gene names on "gene_name" column. This dataframe can be obtained using `loadGTF()` to open an GTF file annotation or `downloadGTF()` the annotation using an web link.
    - The non-coding data, a subset of "annotationData" with only the non coding RNAs.
    
Obs: In the case of the user uses directly the gene name/symbol, both gene_id and gene_name columns could repeat the same gene_symbol.
    
    
The package includes a pre-processed TCGA reduced dataset to test, which can be loaded as follows:

```r
library(lacen)

data("annotation_data")
data("expression_DGE")
data("nc_annotation")
data("raw_expression")
data("traits")
```

## Creating a lacen Object

All data is handled within a `lacen`-type S3 object. The `initLacen` function initializes the object by integrating gene expression data, differential expression analysis, sample conditions, and gene annotation.

```r
lacenObject <- initLacen(annotationData = annotation_data,
                         datCounts = raw_expression,
                         datExpression = expression_DGE,
                         datTraits = traits,
                         ncAnnotation = nc_annotation)
```

## Checking Data Format

To ensure data compatibility, the `checkData` function verifies that all input datasets conform to the expected structure.

```r
checkData(lacenObject)
```

## Filtering and Transforming Data

The `filterTransform` function removes genes with low variation to prepare data for WGCNA. It allows filtering based on differentially expressed genes (DEGs) or the most variable genes (MAD) while applying necessary transformations.

```r
lacenObject <- filterTransform(lacenObject = lacenObject,
                               pThreshold = 0.01,
                               fcThreshold = 1,
                               filterMethod = "DEG")
```

Alternatively, the user can also input a personalized filtered and transformed data frame `lacenObject$datExpr <- filteredDataFrame`, with the sample names as column names and gene names as row names.



## Removing outliers

This function relies on WGCNA to plot the samples cluster Tree, helping to find the height value to exclude outlier samples.

```r
selectOutlierSample(lacenObject, height = FALSE)
```

![Figure 1: Cluster Tree](/home/leo/Documents/LACEN/lacen/figures/1a_clusterTree.png)

When height is provided, the function will return the samples to keep.

```r
selectOutlierSample(lacenObject, height = 270)
```

![Figure 2: Cluster Tree with height](/home/leo/Documents/LACEN/lacen/figures/1b_clusterTree.png)


To confirm the choose height, use `cutOutlierSample` to update it in the lacen object.

```r
lacenObject <- cutOutlierSample(lacenObject, height = 270)
```


## Picking the Beta-Value

To ensure the generated network will follow a scale-free topology, the choose soft-threshold should maximize the model fit ($R^2$) while minimizing the number of connections lost. 

```r
plotSoftThreshold(lacenObject,
                  filename = "2_indicePower.png",
                  maxBlockSize = 16000,
                  plot = TRUE
)
```

![Figure 3: Soft Threshold](/home/leo/Documents/LACEN/lacen/figures/2_indicePower.png)

In this example, the value that maximizes the model fit to a scale-free topology is close to 15.

```r
lacenObject <- selectSoftThreshold(lacenObject = lacenObject,
                                   indicePower = 15)
```

After choosing the soft threshold value, confirm it using `selectSoftThreshold`.

## Bootstrapping

The bootstrap function makes the network n times, removing 1/n genes each repetition. Consequently, the less stable genes can be removed from the analysis.

```lacenObject <- lacenBootstrap(lacenObject = lacenObject)```


![Figure 4: Bootstrap - Module Groups](/home/leo/Documents/LACEN/lacen/figures/3_modgroups.png)

![Figure 5: Bootstrap - Stability](/home/leo/Documents/LACEN/lacen/figures/4_stability_bootstrap.png)


After analyzing the figures, the user can set the value of `setBootstrap` to remove genes less stable than this threshold.

```
lacenObject <- setBootstrap(lacenObject = lacenObject, cutBootstrap = 10)
```

## Final network

The lacen core function automatically generates the final network, enrich the modules based on the given database, and returns the reduced enrichment as a png.

```
lacenObject <- summarizeAndEnrichModules(lacenObject = lacenObject,
                                           #WGCNA parameters
                                           maxBlockSize = 11000,
                                           TOMType = "unsigned",
                                           minModuleSize = 30,
                                           reassignThreshold = 0,
                                           mergeCutHeight = 0.3,
                                           pamRespectsDendro = FALSE,
                                           corType = "bicor",
                                           
                                           #Enrichment analysis parameters
                                           userThreshold = 0.05,
                                           ontology = "BP",
                                           organism = "hsapiens",
                                           orgdb="org.Hs.eg.db",
                                           reducedTermsThreshold=0.7,
                                           filename = "5_enrichedgraph.png",
                                           mod_path = ".", # hidden, no doc
                                           
                                           #Log parameters
                                           log = TRUE, # hidden, no doc
                                           log_path = "log.txt" # hidden, no doc
)
```

![Figure 6: Enriched Modules](/home/leo/Documents/LACEN/lacen/figures/5_enrichedgraph.png)


## Modules Summarization

The modules are summarized in a barplot as its protein-coding/non-coding genes and the trait-module correlation.

```
stackedBarplot(lacenObject = lacenObject,
               filename = "6_stackedplot_desk.png",
               plot = TRUE)
```

## Enriched modules heatmap

This function generates a heatmap of a given module or module/submodule. The figure will exhibits the interconnectedness based on the topological overlap matrix between the lncRNAs and the protein coding genes.

![Figure 7: Module Composition](/home/leo/Documents/LACEN/lacen/figures/6_stackedplot_desk.png)

```
heatmapTopConnectivity(lacenObject = lacenObject,
                       module = 3,
                       submodule = FALSE,
                       hmDimensions = FALSE,
                       filename = "7_heatmap.png", #if false, path.file
                       removeNonDEG = FALSE,
                       outTSV = FALSE,
                       plothm = TRUE)
```

## lncRNA centered network

Based on a lncRNA identifier and a LACEN object containing coexpression network data, and performs functional enrichment analysis based on the genes most highly correlated with the given lncRNA. It identifies the most correlated genes with the specified lncRNA using the topological overlap matrix (TOM), enriches this group for functional terms such as Gene Ontology (GO) biological processes, molecular functions, cellular components, KEGG pathways, and Reactome pathways, and visualizes the resulting network and enrichment graph.

![Figure 8: Module Composition](/home/leo/Documents/LACEN/lacen/figures/6_stackedplot_desk.png)


```
lncRNAEnrich(lncName = "LINC00665",
             lacenObject = lacenObject,
             nGenes = 300,
             nHighlight = 32,
             nGenesNet = 50,
             lncHighlight = FALSE)
```

# Issues to fix

# Fix the name of the variables
# Use only ID or transcript ID, remove gene_ID
# Fix the names of transcript IDs
