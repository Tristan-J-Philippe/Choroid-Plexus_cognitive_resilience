#######Installation
#snRNA-sequencing, snATAC-sequencing, and proteomics analyses were performed in R (v4.2) and Seurat v4.2.0. Spatial transcriptomcis which was performed in R 4.4 and Seurat v5.1.0. Differential expression used limma/edgeR; gene-set enrichment used fgsea with gene sets from MSigDB; single-cell and spatial analyses used Seurat, Signac, singleCellTK, CellChat, and CellularNeighborhoods.
#######
install.packages(c(
  "reticulate",
  "Matrix",
  "data.table",
  "plyr",
  "dplyr",
  "future",
  "ggplot2",
  "cowplot",
  "SoupX",
  "DoubletFinder",
  "qlcMatrix",
  "FastCAR",
  "gridExtra",
  "dsb",
  "cluster",
  "fitdistrplus",
  "harmony",
  "tidyverse",
  "ggVennDiagram",
  "janitor",
  "clustree",
  "igraph",
  "igraphdata",
  "export",
  "lattice",
  "viridis",
  "ggrepel",
  "splines",
  "stringr",
  "patchwork",
  "tibble",
  "tidyr",
  "sf",
  "purrr",
  "readr",
  "lme4",
  "emmeans",
  "devtools",
  "pheatmap",
  "rlang",
  "forcats",
  "circlize",
  "RColorBrewer",
  "future.apply",
  "Seurat",
  "Signac",
  "singleCellTK",
  "CellChat"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GenomicRanges",
  "GenomeInfoDb",
  "org.Hs.eg.db",
  "limma",
  "edgeR",
  "ComplexHeatmap",
  "fgsea",
  "msigdbr",
  "clusterProfiler"
))

devtools::install_github(
  "Nanostring-Biostats/CosMx-Analysis-Scratch-Space",
  subdir = "_code/cellular-neighborhoods"
)

####Software versions for snRNAseq, snATACseq, and proteomics data analysis
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 11 (bullseye)

Matrix products: default
BLAS:   /usr/local/lib/R/lib/libRblas.so
LAPACK: /usr/local/lib/R/lib/libRlapack.so

locale:
  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C             
[9] LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
  [1] splines   grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.66.3                   rtracklayer_1.58.0                Biostrings_2.66.0                 XVector_0.38.0                   
[6] EnsDb.Hsapiens.v86_2.99.0         ensembldb_2.22.0                  AnnotationFilter_1.22.0           GenomicFeatures_1.50.4            TFBSTools_1.36.0                 
[11] JASPAR2020_0.99.10                motifmatchr_1.1.1                 ggsignif_0.6.4                    patchwork_1.1.3                   CellChat_1.6.1                   
[16] GSVA_1.46.0                       clusterProfiler_4.6.2             viridis_0.6.4                     viridisLite_0.4.2                 DESeq2_1.38.3                    
[21] circlize_0.4.16                   msigdbr_7.5.1                     fgsea_1.25.1                      devtools_2.4.5                    usethis_2.2.2                    
[26] ggrepel_0.9.4                     ComplexHeatmap_2.14.0             edgeR_3.40.2                      limma_3.54.2                      lattice_0.21-8                   
[31] ggVennDiagram_1.2.3               lubridate_1.9.2                   forcats_1.0.0                     stringr_1.5.1                     purrr_1.0.2                      
[36] readr_2.1.4                       tidyr_1.3.0                       tidyverse_2.0.0                   janitor_2.2.0                     clustree_0.5.0                   
[41] ggraph_2.1.0                      fitdistrplus_1.1-11               survival_3.5-7                    MASS_7.3-60                       cluster_2.1.4                    
[46] export_0.3.0                      igraphdata_1.0.1                  igraph_1.5.1                      dsb_1.0.3                         gridExtra_2.3                    
[51] FastCAR_0.1.0                     qlcMatrix_0.9.7                   sparsesvd_0.2-2                   slam_0.1-50                       SoupX_1.6.2                      
[56] org.Hs.eg.db_3.16.0               AnnotationDbi_1.60.2              reticulate_1.44.1                 googleCloudStorageR_0.7.0         cowplot_1.1.1                    
[61] ggplot2_3.4.4                     future_1.33.0                     harmony_1.2.4                     Rcpp_1.0.11                       Signac_1.10.0                    
[66] SeuratObject_4.1.3                Seurat_4.3.0.1                    singleCellTK_2.8.0                DelayedArray_0.24.0               SingleCellExperiment_1.20.1      
[71] SummarizedExperiment_1.28.0       Biobase_2.58.0                    MatrixGenerics_1.10.0             matrixStats_1.0.0                 GenomicRanges_1.50.2             
[76] GenomeInfoDb_1.34.9               IRanges_2.32.0                    S4Vectors_0.36.2                  BiocGenerics_0.44.0               tibble_3.2.1                     
[81] dplyr_1.1.4                       plyr_1.8.8                        data.table_1.14.8                 Matrix_1.5-4.1                   

loaded via a namespace (and not attached):
  [1] graphlayouts_1.0.0          pbapply_1.7-2               fontquiver_0.2.1            vctrs_0.6.5                 GSVAdata_1.34.0             blob_1.2.4                  spatstat.data_3.0-1        
[8] later_1.3.1                 DBI_1.1.3                   R.utils_2.12.2              gfonts_0.2.0                rappdirs_0.3.3              uwot_0.1.16                 dqrng_0.3.0                
[15] zlibbioc_1.44.0             htmlwidgets_1.6.2           GlobalOptions_0.1.2         leiden_0.4.3                parallel_4.2.2              irlba_2.3.5.1               tidygraph_1.2.3            
[22] KernSmooth_2.23-22          promises_1.2.1              pkgload_1.3.2.1             graph_1.76.0                RSpectra_0.16-1             fs_1.6.3                    textshaping_0.3.6          
[29] googleAuthR_2.0.1           fastmatch_1.1-3             digest_0.6.33               png_0.1-8                   StratigrapheR_1.2.4         sctransform_0.3.5           scatterpie_0.2.4           
[36] DOSE_3.24.2                 crul_1.4.0                  pkgconfig_2.0.3             GO.db_3.16.0                docopt_0.7.1                gridBase_0.4-7              spatstat.random_3.1-5      
[43] DelayedMatrixStats_1.20.0   ggbeeswarm_0.7.2            iterators_1.0.14            statnet.common_4.9.0        DropletUtils_1.18.1         network_1.18.1              beeswarm_0.4.0             
[50] GetoptLong_1.0.5            xfun_0.40                   zoo_1.8-15                  tidyselect_1.2.0            reshape2_1.4.4              ica_1.0-3                   gson_0.1.0                 
[57] pkgbuild_1.4.2              rlang_1.1.4                 RVenn_1.1.0                 glue_1.6.2                  gdtools_0.3.3               RColorBrewer_1.1-3          registry_0.5-1             
[64] CNEr_1.34.0                 labeling_0.4.3              fontLiberation_0.1.0        httpuv_1.6.11               BiocNeighbors_1.16.0        seqLogo_1.64.0              webshot_0.5.5              
[71] annotate_1.76.0             jsonlite_1.8.7              bit_4.0.5                   fontBitstreamVera_0.1.1     mime_0.12                   systemfonts_1.0.4           Rsamtools_2.14.0           
[78] stringi_1.8.3               processx_3.8.2              RcppRoll_0.3.0              spatstat.sparse_3.0-2       scattermore_1.2             spatstat.explore_3.2-1      yulab.utils_0.2.0          
[85] bitops_1.0-7                cli_3.6.2                   rhdf5filters_1.10.1         RSQLite_2.3.1               rvg_0.3.3                   timechange_0.2.0            officer_0.6.2              
[92] rstudioapi_0.15.0           GenomicAlignments_1.34.1    nlme_3.1-163                qvalue_2.30.0               locfit_1.5-9.8              listenv_0.9.0               miniUI_0.1.1.1             
[99] gridGraphics_0.5-1          R.oo_1.25.0                 urlchecker_1.0.1            httpcode_0.3.0              ggnetwork_0.5.12            dbplyr_2.3.3                sessioninfo_1.2.2          
[106] lifecycle_1.0.4             munsell_0.5.0               R.methodsS3_1.8.2           caTools_1.18.2              ggalluvial_0.12.5           coda_0.19-4                 codetools_0.2-19           
[113] vipor_0.4.5                 lmtest_0.9-40               xtable_1.8-4                ROCR_1.0-11                 BiocManager_1.30.22         abind_1.4-5                 FNN_1.1.3.2                
[120] farver_2.1.1                parallelly_1.36.0           RANN_2.6.1                  aplot_0.2.0                 askpass_1.1                 poweRlaw_1.0.0              BiocIO_1.8.0               
[127] ggtree_3.6.2                RcppAnnoy_0.0.21            goftest_1.2-3               rgl_1.2.1                   profvis_0.3.8               future.apply_1.11.0         svMisc_1.2.3               
[134] tidytree_0.4.5              ellipsis_0.3.2              prettyunits_1.1.1           ggridges_0.5.4              mclust_6.0.0                flextable_0.9.2             diagram_1.6.5              
[141] remotes_2.4.2.1             gargle_1.5.2                spatstat.utils_3.1-2        BiocFileCache_2.6.1         htmltools_0.5.6             NMF_0.28                    yaml_2.3.7                 
[148] utf8_1.2.4                  plotly_4.10.2               XML_3.99-0.14               ggpubr_0.6.0                withr_2.5.2                 scuttle_1.8.4               BiocParallel_1.32.6        
[155] bit64_4.0.5                 rngtools_1.5.2              foreach_1.5.2               ProtGenerics_1.30.0         progressr_0.14.0            GOSemSim_2.24.0             ragg_1.2.5                 
[162] rsvd_1.0.5                  ScaledMatrix_1.6.0          memoise_2.0.1               evaluate_0.21               geneplotter_1.76.0          tzdb_0.4.0                  callr_3.7.3                
[169] ps_1.7.5                    curl_5.0.2                  fansi_1.0.6                 GSEABase_1.60.0             tensor_1.5                  cachem_1.0.8                deldir_1.0-9               
[176] HDO.db_0.99.1               babelgene_22.9              rjson_0.2.21                openxlsx_4.2.5.2            rstatix_0.7.2               stargazer_5.2.3             clue_0.3-64                
[183] tools_4.2.2                 magrittr_2.0.3              RCurl_1.98-1.12             TFMPvalue_0.0.9             car_3.1-2                   ape_5.7-1                   ggplotify_0.1.2            
[190] xml2_1.3.5                  httr_1.4.7                  assertthat_0.2.1            rmarkdown_2.24              globals_0.16.2              R6_2.5.1                    Rhdf5lib_1.20.0            
[197] progress_1.2.2              DirichletMultinomial_1.40.0 KEGGREST_1.38.0             treeio_1.22.0               gtools_3.9.4                shape_1.4.6                 beachmat_2.14.2            
[204] sna_2.7-1                   HDF5Array_1.26.0            BiocSingular_1.14.0         ggrastr_1.0.2               rhdf5_2.42.1                carData_3.0-5               snakecase_0.11.0           
[211] ggfun_0.1.2                 colorspace_2.1-0            generics_0.1.3              base64enc_0.1-3             pillar_1.9.0                tweenr_2.0.2                sp_2.0-0                   
[218] uuid_1.1-0                  GenomeInfoDbData_1.2.9      gtable_0.3.6                zip_2.3.0                   devEMF_4.4-2                restfulr_0.0.15             knitr_1.43                 
[225] biomaRt_2.54.1              shadowtext_0.1.4            fastmap_1.1.1               doParallel_1.0.17           ComplexUpset_1.3.3          fishpond_2.4.1              broom_1.0.5                
[232] filelock_1.0.2              openssl_2.1.0               scales_1.3.0                backports_1.4.1             enrichplot_1.18.4           hms_1.1.3                   ggforce_0.4.1              
[239] Rtsne_0.17                  shiny_1.7.5                 polyclip_1.10-4             lazyeval_0.2.2              crayon_1.5.2                downloader_0.4              sparseMatrixStats_1.10.0   
[246] svglite_2.1.1               reshape_0.8.9               compiler_4.2.2              spatstat.geom_3.2-4   





#####Software versions for Spatial Transcriptomic data analysis
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Debian GNU/Linux 12 (bookworm)

Matrix products: default
attached base packages:
[1] splines   grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggalluvial_0.12.5                NMF_0.28                         rngtools_1.5.2                   registry_0.5-1                   presto_1.0.0                    
[6] CellularNeighborhoods_0.0.0.9000 ggsignif_0.6.4                   circlize_0.4.16                  CellChat_2.1.2                   ggrastr_1.0.2                   
[11] GSVA_1.46.0                      msigdbr_10.0.1                   fgsea_1.30.0                     ggrepel_0.9.5                    edgeR_4.2.1                     
[16] limma_3.60.4                     lattice_0.22-6                   ggVennDiagram_1.5.2              lubridate_1.9.3                  forcats_1.0.0                   
[21] stringr_1.5.1                    purrr_1.0.2                      readr_2.1.5                      tidyr_1.3.1                      tibble_3.2.1                    
[26] tidyverse_2.0.0                  janitor_2.2.0                    igraphdata_1.0.1                 igraph_2.0.3                     clustree_0.5.1                  
[31] ggraph_2.2.1                     harmony_1.2.0                    Rcpp_1.0.13                      fitdistrplus_1.2-1               survival_3.6-4                  
[36] MASS_7.3-60.2                    cluster_2.1.6                    ComplexHeatmap_2.20.0            patchwork_1.2.0                  tonsilref.SeuratData_2.0.0      
[41] stxBrain.SeuratData_0.1.2        pbmcref.SeuratData_1.0.0         pancreasref.SeuratData_1.0.0     mousecortexref.SeuratData_1.0.0  lungref.SeuratData_2.0.0        
[46] kidneyref.SeuratData_1.0.2       humancortexref.SeuratData_1.0.0  heartref.SeuratData_1.0.0        fetusref.SeuratData_1.0.0        bonemarrowref.SeuratData_1.0.0  
[51] adiposeref.SeuratData_1.0.0      SeuratData_0.2.2.9001            Azimuth_0.5.0                    shinyBS_0.61.1                   export_0.3.0                    
[56] dsb_1.0.4                        gridExtra_2.3                    FastCAR_0.1.0                    qlcMatrix_0.9.8                  sparsesvd_0.2-2                 
[61] slam_0.1-51                      DoubletFinder_2.0.4              SoupX_1.6.2                      org.Hs.eg.db_3.19.1              AnnotationDbi_1.66.0            
[66] reticulate_1.45.0                googleCloudStorageR_0.7.0        cowplot_1.1.3                    ggplot2_3.5.1                    future_1.33.2                   
[71] Signac_1.13.0                    SeuratWrappers_0.3.5             Seurat_5.1.0                     SeuratObject_5.0.2               sp_2.1-4                        
[76] singleCellTK_2.14.0              DelayedArray_0.30.1              SparseArray_1.4.8                S4Arrays_1.4.1                   abind_1.4-5                     
[81] SingleCellExperiment_1.26.0      SummarizedExperiment_1.34.0      Biobase_2.64.0                   MatrixGenerics_1.16.0            matrixStats_1.3.0               
[86] GenomicRanges_1.56.1             GenomeInfoDb_1.40.1              IRanges_2.38.1                   S4Vectors_0.42.1                 BiocGenerics_0.50.0             
[91] dplyr_1.1.4                      plyr_1.8.9                       data.table_1.15.4                Matrix_1.7-0                    

loaded via a namespace (and not attached):
  [1] graph_1.82.0                      ica_1.0-3                         plotly_4.10.4                     zlibbioc_1.50.0                   tidyselect_1.2.1                 
[6] bit_4.0.5                         doParallel_1.0.17                 clue_0.3-65                       rjson_0.2.21                      blob_1.2.4                       
[11] parallel_4.4.1                    seqLogo_1.70.0                    png_0.1-8                         cli_3.6.3                         ProtGenerics_1.36.0              
[16] askpass_1.2.0                     openssl_2.2.0                     goftest_1.2-3                     gargle_1.5.2                      textshaping_0.4.0                
[21] BiocIO_1.14.0                     officer_0.6.6                     ggnetwork_0.5.13                  BiocNeighbors_1.22.0              uwot_0.2.2                       
[26] curl_5.2.1                        mime_0.12                         evaluate_0.24.0                   devEMF_4.4-2                      leiden_0.4.3.1                   
[31] stringi_1.8.4                     backports_1.5.0                   XML_3.99-0.17                     httpuv_1.6.15                     magrittr_2.0.3                   
[36] rappdirs_0.3.3                    RcppRoll_0.3.1                    mclust_6.1.1                      GSVAdata_1.40.0                   webshot_0.5.5                    
[41] DT_0.33                           sctransform_0.4.1                 ggbeeswarm_0.7.2                  DBI_1.2.3                         HDF5Array_1.32.0                 
[46] jquerylib_0.1.4                   withr_3.0.0                       systemfonts_1.1.0                 lmtest_0.9-40                     GSEABase_1.66.0                  
[51] tidygraph_1.3.1                   rtracklayer_1.64.0                BiocManager_1.30.27               htmlwidgets_1.6.4                 fs_1.6.4                         
[56] statnet.common_4.9.0              labeling_0.4.3                    cellranger_1.1.0                  annotate_1.82.0                   zoo_1.8-14                       
[61] JASPAR2020_0.99.10                XVector_0.44.0                    knitr_1.48                        network_1.18.2                    TFBSTools_1.42.0                 
[66] UCSC.utils_1.0.0                  docopt_0.7.1                      timechange_0.3.0                  TFMPvalue_0.0.9                   foreach_1.5.2                    
[71] fansi_1.0.6                       caTools_1.18.2                    rhdf5_2.48.0                      pwalign_1.0.0                     rvg_0.3.3                        
[76] R.oo_1.26.0                       poweRlaw_0.80.0                   RSpectra_0.16-2                   irlba_2.3.5.1                     fastDummies_1.7.3                
[81] lazyeval_0.2.2                    yaml_2.3.9                        scattermore_1.2                   crayon_1.5.3                      RcppAnnoy_0.0.22                 
[86] tweenr_2.0.3                      RColorBrewer_1.1-3                progressr_0.14.0                  later_1.3.2                       gfonts_0.2.0                     
[91] ggridges_0.5.6                    codetools_0.2-20                  base64enc_0.1-3                   GlobalOptions_0.1.2               KEGGREST_1.44.1                  
[96] Rtsne_0.17                        shape_1.4.6.1                     gdtools_0.3.7                     Rsamtools_2.20.0                  pkgconfig_2.0.3                  
[101] ggpubr_0.6.0                      xml2_1.3.6                        spatstat.univar_3.0-0             GenomicAlignments_1.40.0          gridBase_0.4-7                   
[106] spatstat.sparse_3.1-0             BSgenome_1.72.0                   viridisLite_0.4.2                 xtable_1.8-4                      car_3.1-2                        
[111] httr_1.4.7                        tools_4.4.1                       eds_1.6.0                         globals_0.16.3                    beeswarm_0.4.0                   
[116] broom_1.0.6                       nlme_3.1-164                      googleAuthR_2.0.2                 hdf5r_1.3.11                      shinyjs_2.1.0                    
[121] ComplexUpset_1.3.3                rgl_1.3.1                         assertthat_0.2.1                  digest_0.6.36                     farver_2.1.2                     
[126] tzdb_0.4.0                        AnnotationFilter_1.28.0           reshape2_1.4.4                    viridis_0.6.5                     DirichletMultinomial_1.46.0      
[131] crul_1.5.0                        glue_1.7.0                        cachem_1.1.0                      polyclip_1.10-7                   generics_0.1.3                   
[136] Biostrings_2.72.1                 googledrive_2.1.1                 parallelly_1.37.1                 statmod_1.5.0                     RcppHNSW_0.6.0                   
[141] ragg_1.3.2                        ScaledMatrix_1.12.0               carData_3.0-5                     fontBitstreamVera_0.1.1           pbapply_1.7-2                    
[146] spam_2.10-0                       dqrng_0.4.1                       utf8_1.2.4                        stargazer_5.2.3                   graphlayouts_1.1.1               
[151] gtools_3.9.5                      httpcode_0.3.0                    shiny_1.8.1.1                     GenomeInfoDbData_1.2.12           R.utils_2.12.3                   
[156] rhdf5filters_1.16.0               RCurl_1.98-1.16                   memoise_2.0.1                     rmarkdown_2.27                    scales_1.3.0                     
[161] R.methodsS3_1.8.2                 googlesheets4_1.1.1               svglite_2.1.3                     RANN_2.6.1                        fontLiberation_0.1.0             
[166] spatstat.data_3.1-2               rstudioapi_0.16.0                 spatstat.utils_3.0-5              hms_1.1.3                         munsell_0.5.1                    
[171] colorspace_2.1-0                  FNN_1.1.4                         rlang_1.1.4                       DelayedMatrixStats_1.26.0         sparseMatrixStats_1.16.0         
[176] dotCall64_1.1-1                   shinydashboard_0.7.2              ggforce_0.4.2                     scuttle_1.14.0                    xfun_0.46                        
[181] coda_0.19-4.1                     sna_2.7-2                         CNEr_1.40.0                       remotes_2.5.0                     iterators_1.0.14                 
[186] EnsDb.Hsapiens.v86_2.99.0         Rhdf5lib_1.26.0                   bitops_1.0-7                      promises_1.3.0                    RSQLite_2.3.7                    
[191] openxlsx_4.2.6.1                  GO.db_3.19.1                      compiler_4.4.1                    beachmat_2.20.0                   listenv_0.9.1                    
[196] BSgenome.Hsapiens.UCSC.hg38_1.4.5 fontquiver_0.2.1                  BiocSingular_1.20.0               tensor_1.5                        uuid_1.2-0                       
[201] BiocParallel_1.38.0               babelgene_22.9                    spatstat.random_3.3-1             R6_2.5.1                          rstatix_0.7.2                    
[206] fastmap_1.2.0                     fastmatch_1.1-4                   vipor_0.4.7                       ensembldb_2.28.0                  ROCR_1.0-11                      
[211] SeuratDisk_0.0.0.9021             rsvd_1.0.5                        gtable_0.3.5                      KernSmooth_2.23-24                miniUI_0.1.1.1                   
[216] deldir_2.0-4                      htmltools_0.5.8.1                 bit64_4.0.5                       spatstat.explore_3.3-1            lifecycle_1.0.4                  
[221] msigdbdf_24.1.0                   zip_2.3.1                         restfulr_0.0.15                   sass_0.4.9                        vctrs_0.6.5                      
[226] snakecase_0.11.1                  flextable_0.9.6                   spatstat.geom_3.3-2               future.apply_1.11.2               pracma_2.4.4                     
[231] bslib_0.7.0                       pillar_1.9.0                      GenomicFeatures_1.56.0            DropletUtils_1.24.0               locfit_1.5-9.10                  
[236] jsonlite_1.8.8                    GetoptLong_1.0.5 