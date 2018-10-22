# http://127.0.0.1:20734/help/library/clusterProfiler/doc/clusterProfiler.R

library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)
library("pathview")

x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(ids)

data(gcSample)
hg <- gcSample[[1]]
head(hg)

eg2np <- bitr_kegg(hg, fromType='kegg', toType='ncbi-proteinid', organism='hsa')
head(eg2np)

bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)

head(ego2)

ego2b <- setReadable(ego2, OrgDb = org.Hs.eg.db)

head(ego2b)

ego3 <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                nPerm        = 1000,
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)

head(ego3)

ego3b <- setReadable(ego3, OrgDb = org.Hs.eg.db)

head(ego3b)

# kegg

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

mkk <- enrichMKEGG(gene = gene,
                  organism = 'hsa')

mkk

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa')
                 
head(mkk2)

david <- enrichDAVID(gene = gene,
                       idType = "ENTREZ_GENE_ID",
                       listType = "Gene",
                       annotation = "KEGG_PATHWAY",
                       david.user = "clusterProfiler@hku.hk")

gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

barplot(ggo, drop=TRUE, showCategory=12)
barplot(ego, showCategory=8)

dotplot(ego)

emapplot(ego)

cnetplot(ego, categorySize="pvalue", foldChange=geneList)

goplot(ego)

gseaplot(kk2, geneSetID = "hsa04145")

browseKEGG(kk, 'hsa04110')
library("pathview")

hsa04110 <- pathview(gene.data  = geneList,
                       pathway.id = "hsa04110",
                       species    = "hsa",
                       limit      = list(gene=max(abs(geneList)), cpd=1))

data(gcSample)
lapply(gcSample, head)

ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))

mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))

dotplot(ck)

dotplot(formula_res)

dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)

sessionInfo()

# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pathview_1.20.0       GSEABase_1.42.0       graph_1.58.2          annotate_1.58.0       XML_3.98-1.16         BiocInstaller_1.30.0  GO.db_3.6.0          
# [8] DOSE_3.6.1            clusterProfiler_3.8.1 org.Hs.eg.db_3.6.0    AnnotationDbi_1.42.1  IRanges_2.14.12       S4Vectors_0.18.3      Biobase_2.40.0       
# [15] BiocGenerics_0.26.0   RColorBrewer_1.1-2    pheatmap_1.0.10       DT_0.4                openxlsx_4.1.0        shinyBS_0.61          shiny_1.1.0          
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-6        enrichplot_1.0.2    bit64_0.9-7         httr_1.3.1          UpSetR_1.3.3        Rgraphviz_2.24.0    tools_3.5.1         R6_2.3.0           
# [9] DBI_1.0.0           lazyeval_0.2.1      colorspace_1.3-2    tidyselect_0.2.5    gridExtra_2.3       bit_1.1-14          compiler_3.5.1      labeling_0.3       
# [17] KEGGgraph_1.40.0    scales_1.0.0        ggridges_0.5.1      stringr_1.3.1       digest_0.6.18       XVector_0.20.0      pkgconfig_2.0.2     htmltools_0.3.6    
# [25] htmlwidgets_1.3     rlang_0.2.2         rstudioapi_0.8      RSQLite_2.1.1       bindr_0.1.1         farver_1.0          jsonlite_1.5        crosstalk_1.0.0    
# [33] BiocParallel_1.14.2 GOSemSim_2.6.2      dplyr_0.7.7         zip_1.0.0           RCurl_1.95-4.11     magrittr_1.5        Matrix_1.2-14       Rcpp_0.12.19       
# [41] munsell_0.5.0       viridis_0.5.1       stringi_1.2.4       yaml_2.2.0          ggraph_1.0.2        zlibbioc_1.26.0     MASS_7.3-51         plyr_1.8.4         
# [49] qvalue_2.12.0       grid_3.5.1          blob_1.1.1          promises_1.0.1      ggrepel_0.8.0       DO.db_2.9           crayon_1.3.4        lattice_0.20-35    
# [57] Biostrings_2.48.0   cowplot_0.9.3       splines_3.5.1       KEGGREST_1.20.2     pillar_1.3.0        fgsea_1.6.0         igraph_1.2.2        reshape2_1.4.3     
# [65] fastmatch_1.1-0     glue_1.3.0          data.table_1.11.8   png_0.1-7           tweenr_1.0.0        httpuv_1.4.5        gtable_0.2.0        purrr_0.2.5        
# [73] tidyr_0.8.1         assertthat_0.2.0    ggplot2_3.0.0       ggforce_0.1.3       mime_0.6            xtable_1.8-3        later_0.7.5         rsconnect_0.8.8    
# [81] viridisLite_0.3.0   tibble_1.4.2        rvcheck_0.1.1       memoise_1.1.0       units_0.6-1         bindrcpp_0.2.2    
