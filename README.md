# lacen: An R package for integrative analysis of lncRNA Expression Networks

## Introduction

`lacen` is a dedicated R software package developed to facilitate integrative analysis of gene expression data, focusing specifically on the elucidation of functional roles for long non-coding RNAs (lncRNAs) through weighted gene co-expression network analysis (WGCNA). 

The package offers an intuitive, structured, and robust workflow that seamlessly integrates data handling, transformation, and network construction steps and includes advanced stability checks via bootstrap analyses. `lacen` is designed to be user-friendly and accessible, even for researchers without extensive computational or programming experience, guiding users through a clear and consistent analysis pathway.

Although initially optimized using human cancer RNA-seq datasets, the package is not limited to oncological studies and can handle data from numerous experimental contexts and organisms. The only strict requirements are that the input RNA-seq expression profiles must represent at least two sample conditions (for instance, two different biological or clinical groups), and input files must conform to the required format. Beyond human datasets, `lacen` supports most widely-researched organisms included in [Bioconductor's OrgDb annotation packages](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb).

---

## Installation

`lacen` is currently not available in the Comprehensive R Archive Network (CRAN). Instead, the development version can be installed directly via BiocManager, ensuring the appropriate installation of all Bioconductor dependencies:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sanches-leo/lacen")
```

---

## Input Data Requirements

The `lacen` package requires five input datasets to commence analysis:

1. **Expression Counts Dataframe** (`datCounts`): A matrix/dataframe composed of raw RNA-seq counts, structured with genes as rows and samples as columns. Gene identifiers (gene symbols or IDs) should be used as row names, and sample identifiers should appear as column names.

2. **Experimental Conditions Dataframe (Traits)** (`datTraits`): A two-column dataframe specifying sample identifiers and their respective experimental condition labels. Columns must be explicitly named `"Sample"` (sample identifier) and `"Trait"` (numeric condition/group identifier), respectively. Experimental conditions must be numerically encoded (e.g., "0" and "1" for two different groups).

3. **Differential Expression Results** (`datExpression`): A dataframe summarizing differential expression analysis results. By default, standard outputs from tools such as Limma can be directly used. Alternatively, the dataframe can include three required columns: `"ID"` (gene identifier), `"log2FC"` (log2 fold-change between conditions), and `"pvalue"` (statistical significance), arranged in this exact order.

4. **Gene Annotation Data** (`annotationData`): A two-column dataframe linking gene identifiers to gene names, with columns specifically labeled `"gene_id"` and `"gene_name"`. Annotation data can be directly loaded from a GTF file using included utility functions (`loadGTF()` or `downloadGTF()` functions). 
   - **Note**: if analysis is based solely on gene symbols, both columns might replicate the same symbol.

5. **Non-coding RNA Annotation Dataframe** (`ncAnnotation`): A specific subset of the annotation data (`annotationData`) containing only long non-coding RNA entries.


---

## Example Data (Test Dataset)

The `lacen` package provides a pre-processed and reduced subset of TCGA breast carcinoma (BRCA) RNA-seq data to facilitate familiarity and testing. Load example datasets conveniently using:

```r
library(lacen)
data("annotation_data")
data("expression_DGE")
data("nc_annotation")
data("raw_expression")
data("traits")
```

**Note**: The provided example datasets represent randomly selected paired tumor-adjacent tissue samples derived from the TCGA BRCA project, processed using the STAR alignment pipeline available through the TCGAbiolinks package. Procedures to obtain and preprocess the TCGA data can be provided by request.

---

## Creating the `lacen` Object

All analysis steps are encapsulated within a convenient S3 object. Begin by consolidating your data into one structured object using `initLacen`:

```r
lacenObject <- initLacen(annotationData = annotation_data,
                         datCounts = raw_expression,
                         datExpression = expression_DGE,
                         datTraits = traits,
                         ncAnnotation = nc_annotation)
```

---

## Data Integrity Check

Ensure your data is compatible with `lacen` standards by running:

```r
checkData(lacenObject)
```

---

## Data Filtering and Transformation

The `filterTransform` function prepares the data for network construction by filtering low-variation genes. Users can select filtering criteria based on either differentially expressed genes (DEGs) or the median absolute deviation (MAD):

```r
lacenObject <- filterTransform(lacenObject = lacenObject,
                               pThreshold = 0.01,
                               fcThreshold = 1,
                               filterMethod = "DEG")
```

Alternatively, manually filtered data frames can be directly integrated:

```r
lacenObject$datExpr <- filteredDataFrame
```

Ensure your data frame maintains genes as rows and samples as columns.

---

## Outlier Detection and Removal

To visualize and identify outlier samples, generate a hierarchical cluster dendrogram:

```r
selectOutlierSample(lacenObject, height = FALSE)
```

![Figure 1: Cluster Tree](figures/1a_clusterTree.png)

Adjust the `height` parameter to define a threshold for removing outlier samples. For example, if `height = 270` identifies suitable outliers:

```r
selectOutlierSample(lacenObject, height = 270)
```

![Figure 2: Cluster Tree with height](figures/1b_clusterTree.png)

To apply this threshold and update your object:

```r
lacenObject <- cutOutlierSample(lacenObject, height = 270)
```

---

## Selecting the Optimal Soft-Threshold (Beta)

Selecting an appropriate soft-threshold power ensures a scale-free topology for your network:

```r
plotSoftThreshold(lacenObject,
                  filename = "soft_threshold_plot.png")
```

![Figure 3: Soft Threshold](figures/2_indicePower.png)

Inspect the resulting plot and choose the power value maximizing the model fit ($R^2$). Suppose the optimal value is 15:

```r
lacenObject <- selectSoftThreshold(lacenObject, indicePower = 15)
```

---

## Bootstrap Stability Analysis

Assess the robustness of the gene network using the bootstrap method. This approach repeatedly builds networks, omitting subsets of genes to identify consistently stable modules:

```r
lacenObject <- lacenBootstrap(lacenObject = lacenObject)
```

![Figure 4: Bootstrap - Module Groups](figures/3_modgroups.png)

![Figure 5: Bootstrap - Stability](figures/4_stability_bootstrap.png)

Review the provided stability plots, then specify a stability cutoff to retain robust genes:

```r
lacenObject <- setBootstrap(lacenObject = lacenObject, cutBootstrap = 0.1)
```

---

## Final Network Construction and Enrichment Analysis

The core function of LACEN automatically constructs the final network, performs functional enrichment of the identified modules using the specified database, and generates a visual representation of the reduced enrichment results as a PNG file.

```r
lacenObject <- summarizeAndEnrichModules(lacenObject = lacenObject)
```

The generated tree plot provides a structured summary of the functional enrichment of the network modules, containing the following key information:

- Module Number – The identifier for the coexpression module.
- Module-Trait Association – The correlation value between the module and the phenotype/trait of interest, along with its associated p-value.
- Submodule Number – Identifies distinct submodules within a module.
- Reduced Enrichment Description – The summarized functional categories significantly enriched in each submodule after redundancy reduction.
- lncRNA Composition – The number of differentially expressed lncRNAs (DEGs) among the highly connected lncRNAs within the submodule.

![Figure 6: Enriched Modules](figures/5_enrichedgraph.png)

> **Note:** Since lncRNAs often lack functional annotations in standard enrichment databases, LACEN applies a guilt-by-association approach. This assumes that an lncRNA may share similar biological roles with its coexpressed protein-coding genes within the same module or submodule. Given that highly connected lncRNAs can be associated with multiple submodules, there may be overlaps between lncRNAs across different submodules. This reflects the complexity of lncRNA functional associations and provides additional insights into their potential biological roles.


---

## Module Composition Visualization

Visualize module composition by protein-coding and lncRNA content alongside trait-module relationships:

```r
stackedBarplot(lacenObject, filename = "modules_summary.png", plot = TRUE)
```

![Figure 7: Module Composition](figures/6_stackedplot_desk.png)

---

## Module-specific Connectivity Heatmaps

This function generates a heatmap for a selected module or module/submodule, visualizing the interconnectedness of genes based on the topological overlap matrix (TOM). The heatmap highlights the coexpression relationships between lncRNAs and protein-coding genes, helping to identify highly connected lncRNAs within the network.

By analyzing these highly connected lncRNAs, users can investigate their co-expressed protein-coding genes and infer potential functional roles using a guilt-by-association approach.

For example, to examine Module 6, a small module with a strong correlation to the trait of interest, run:

```r
heatmapTopConnectivity(lacenObject = lacenObject,
                       module = 6,
                       submodule = 1)
                       
```

![Figure 8: Connectivity Heatmap](figures/7_heatmap.png)

---

## lncRNA-Centric Network Analysis

For cases where an lncRNA does not cluster into any module or is filtered out by the standard module selection criteria, LACEN provides an alternative approach to analyze its connectivity. This function identifies the most connected genes associated with a specific lncRNA, generates a subnetwork visualization, and performs functional enrichment analysis.

Using a given lncRNA identifier and a LACEN object containing the coexpression network data, the function:

- Identifies the most highly correlated genes with the selected lncRNA using the topological overlap matrix (TOM).
- Performs functional enrichment analysis of this gene set, identifying enriched terms from:
    - Gene Ontology (GO): Biological Processes, Molecular Functions, Cellular Components.
    - Pathways: KEGG and Reactome.
- Generates a network plot showing the strongest coexpression relationships.
- Visualizes the enrichment results, helping to infer potential functions of the lncRNA.

To analyze the lncRNA "LINC01055", run:

```r
lncRNAEnrich(lncName = "LINC01055",
             lacenObject = lacenObject)
```

![Figure 9: LINC01055 network](figures/LINC01055_net.png)

![Figure 10: LINC01055 enrichment](figures/LINC01055_enr.png)


## R Package Dependencies

`lacen` depends on several widely accessible bioinformatics and visualization packages:

- **Core Bioinformatics and analysis**: `WGCNA`, `limma`, `rtracklayer`, `gprofiler2`, `rrvgo`, `org.Hs.eg.db`
- **Parallelization and Speedup**: `foreach`, `doParallel`, `fastcluster`
- **Visualization**: `ggplot2`, `ggraph`, `igraph`, `scatterpie`

---

## Session and Environment Information

```
R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Sao_Paulo
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] lacen_0.0.3

loaded via a namespace (and not attached):
  [1] splines_4.4.3               later_1.4.1                 BiocIO_1.14.0               bitops_1.0-9               
  [5] tibble_3.2.1                R.oo_1.27.0                 polyclip_1.10-7             preprocessCore_1.66.0      
  [9] XML_3.99-0.18               rpart_4.1.24                lifecycle_1.0.4             httr2_1.1.0                
 [13] fastcluster_1.2.6           doParallel_1.0.17           NLP_0.3-2                   lattice_0.22-5             
 [17] MASS_7.3-64                 backports_1.5.0             magrittr_2.0.3              limma_3.60.6               
 [21] Hmisc_5.2-2                 plotly_4.10.4               rmarkdown_2.29              yaml_2.3.10                
 [25] httpuv_1.6.15               askpass_1.2.1               reticulate_1.40.0           DBI_1.2.3                  
 [29] RColorBrewer_1.1-3          abind_1.4-8                 zlibbioc_1.50.0             GenomicRanges_1.56.2       
 [33] purrr_1.0.4                 R.utils_2.12.3              ggraph_2.2.1                BiocGenerics_0.50.0        
 [37] RCurl_1.98-1.16             yulab.utils_0.2.0           nnet_7.3-20                 tweenr_2.0.3               
 [41] rappdirs_0.3.3              GenomeInfoDbData_1.2.12     rrvgo_1.17.0                IRanges_2.38.1             
 [45] S4Vectors_0.42.1            tm_0.7-15                   ggrepel_0.9.6               pheatmap_1.0.12            
 [49] umap_0.2.10.0               RSpectra_0.16-2             codetools_0.2-20            DelayedArray_0.30.1        
 [53] xml2_1.3.6                  ggforce_0.4.2               tidyselect_1.2.1            UCSC.utils_1.1.0           
 [57] farver_2.1.2                viridis_0.6.5               matrixStats_1.5.0           stats4_4.4.3               
 [61] dynamicTreeCut_1.63-1       base64enc_0.1-3             GenomicAlignments_1.40.0    jsonlite_1.8.9             
 [65] tidygraph_1.3.1             Formula_1.2-5               survival_3.8-3              iterators_1.0.14           
 [69] systemfonts_1.2.1           foreach_1.5.2               tools_4.4.3                 ragg_1.3.3                 
 [73] Rcpp_1.0.14                 glue_1.8.0                  gridExtra_2.3               SparseArray_1.4.8          
 [77] xfun_0.50                   MatrixGenerics_1.16.0       GenomeInfoDb_1.40.1         dplyr_1.1.4                
 [81] withr_3.0.2                 fastmap_1.2.0               openssl_2.3.2               digest_0.6.37              
 [85] R6_2.6.0                    mime_0.12                   textshaping_1.0.0           colorspace_2.1-1           
 [89] GO.db_3.19.1                RSQLite_2.3.9               R.methodsS3_1.8.2           tidyr_1.3.1                
 [93] generics_0.1.3              data.table_1.16.4           rtracklayer_1.64.0          graphlayouts_1.2.2         
 [97] httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.4.1              scatterpie_0.2.4           
[101] scatterplot3d_0.3-44        pkgconfig_2.0.3             gtable_0.3.6                blob_1.2.4                 
[105] impute_1.78.0               XVector_0.44.0              htmltools_0.5.8.1           scales_1.3.0               
[109] Biobase_2.64.0              png_0.1-8                   wordcloud_2.6               ggfun_0.1.8                
[113] knitr_1.49                  rstudioapi_0.17.1           rjson_0.2.23                checkmate_2.3.2            
[117] curl_6.2.0                  org.Hs.eg.db_3.19.1         cachem_1.1.0                Polychrome_1.5.1           
[121] stringr_1.5.1               parallel_4.4.3              foreign_0.8-88              AnnotationDbi_1.66.0       
[125] treemap_2.4-4               restfulr_0.0.15             pillar_1.10.1               grid_4.4.3                 
[129] vctrs_0.6.5                 slam_0.1-55                 promises_1.3.2              xtable_1.8-4               
[133] cluster_2.1.8               htmlTable_2.4.3             evaluate_1.0.3              cli_3.6.4                  
[137] compiler_4.4.3              Rsamtools_2.20.0            rlang_1.1.5                 crayon_1.5.3               
[141] labeling_0.4.3              gprofiler2_0.2.3            fs_1.6.5                    stringi_1.8.4              
[145] viridisLite_0.4.2           WGCNA_1.73                  gridBase_0.4-7              BiocParallel_1.38.0        
[149] munsell_0.5.1               Biostrings_2.72.1           lazyeval_0.2.2              GOSemSim_2.30.2            
[153] Matrix_1.7-2                bit64_4.6.0-1               ggplot2_3.5.1               KEGGREST_1.44.1            
[157] statmod_1.5.0               shiny_1.10.0                SummarizedExperiment_1.34.0 igraph_2.1.4               
[161] memoise_2.0.1               bit_4.5.0.1   
```

