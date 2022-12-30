###########################
## PCA ##
###########################

#install.package("devtools",dependencies=T)
#install_github("vqv/ggbiplot")
library(devtools)
library(ggbiplot)
library(ggcorrplot)
library(gplots)
library(readxl)
library(readr)
library(marray)
library(RColorBrewer)

setwd("~/Documents/GitHub/Nucleomics-VIB/Shiny-apps/counts2pca/Data")
ExpNr <- "3949"
RNAseqCounts <- read_excel("exp3949-RNAseqCounts.xlsx")
sample.groups <- read_delim("~/Documents/GitHub/Nucleomics-VIB/Shiny-apps/fpkm2pca/Data/3949-sample-groups.txt", 
                            show_col_types = FALSE,
                            delim = "\t", 
                            escape_double = FALSE, 
                            col_names = FALSE, 
                            trim_ws = TRUE)
colnames(sample.groups) <- c("name", "group")

sample.groups$labels <- sapply(sample.groups$name, function(strings){
  ind = unlist(gregexpr(pattern = "@", text = strings))
  if (length(ind) < 3){NA}
  else{substr(strings, 1, ind[length(ind) - 2] - 1)}
})

sample.groups$labels <- gsub("@", "_", sample.groups$labels)

chromosome.col <- which(colnames(RNAseqCounts)==as.vector("Chromosome"))
counts <- RNAseqCounts[,3:(chromosome.col-1)]
counts <- counts[apply(counts, 1, var) != 0,]

# sample correlation
colnames(counts) <- sapply(colnames(counts), function(strings){
  ind = unlist(gregexpr(pattern = "@", text = strings))
  if (length(ind) < 3){NA}
  else{substr(strings, 1, ind[length(ind) - 2] - 1)}
})
colnames(counts) <-  gsub("@", "_", colnames(counts))

# NC plot (gplots)
Cor <- cor(counts, use="pairwise.complete.obs", method="spearman")
my.pal <- maPalette(low="green", high="red", mid="yellow", k=69)

heatmap.2(Cor,
          col=my.pal,
          cexRow=0.7,cexCol=0.7,
          trace="none",scale="none")

#corr <- round(cor(counts), 1)
# Compute a matrix of correlation p-values
#p.mat <- cor_pmat(counts)
#head(p.mat[, 1:4])

# print correlation heatmap
# Add correlation significance level

#g <- ggcorrplot(corr, hc.order = TRUE, 
#                type = "lower",
#                p.mat = p.mat,
#                outline.col = "white",
#                colors = c("#6D9EC1", "white", "#E46726"),
#                tl.cex = 8)
#print(g)

# PCA output table
counts.pca <- prcomp(t(counts), scale. = TRUE)
pca.table <- as.data.frame(counts.pca$x[,1:4])
rownames(pca.table) <- sample.groups$labels
pca.table<-cbind(pca.table,"Group"=sample.groups$group)

# ggbiplot can only draw ellipses if there are at least 5 samples
if (length(rownames(counts.pca$x)) > 4){
  drawEllipses <- TRUE
}else{
  drawEllipses <- FALSE
}

# print PCA with points
g <- ggbiplot(counts.pca, 
              choices = c(1,2), 
              obs.scale = 1, 
              var.scale = 1, 
              var.axes = F,
              groups = sample.groups$group,
              ellipse = drawEllipses, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
if (length(levels(sample.groups))>7){
  g <- g + guides(col = guide_legend(ncol = 5, byrow = TRUE))
}
print(g)

# Print PCA with labels
g <- ggbiplot(counts.pca, 
              choices = c(1,2), 
              obs.scale = 1, 
              var.scale = 1, 
              var.axes = F,
              groups = sample.groups$group,
              ellipse = drawEllipses, 
              circle = TRUE,
              labels = sample.groups$labels
              )
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
if (length(levels(sample.groups))>7){
  g <- g + guides(col = guide_legend(ncol = 5, byrow = TRUE))
}
print(g)
