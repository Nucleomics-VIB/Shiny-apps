###########################
## PCA ##
###########################

#install.package("devtools",dependencies=T)
#install_github("ggbiplot", "vqv")    
library(devtools)
library(ggbiplot)

#counts.pca <- prcomp(t(counts(PpData.forExcel)), scale. = TRUE)
#sample.info.file<-file.path(DIR,"Data",paste(ExpNr,"-sample-groups.txt",sep=""))

plot.pca <- function(counts.pca="", 
                     sample.groups="",
                     SessionID="",
                     ExpNr="",
                     LogFile="logfile.log"
){
  
  # simplify sample names, preserving unique names
  sample.labels<-c() 
  i<-1
  while (length(unique(sample.labels))!=length(unique(rownames(counts.pca$x)))){
    i<-i+1
    sample.labels<-as.character(sapply(rownames(counts.pca$x),function(x){as.character(paste(strsplit(x,split="@")[[1]][1:i],collapse="@"))}))
  }       
  sample.labels<-gsub(paste("@",ExpNr,sep=""),"",sample.labels)
  
  # Output files
  pca.output.file<-file.path("Images",paste(SessionID,"-summarization-PCA.png",sep=""))
  pca.output.lbl.file<-file.path("Images",paste(SessionID,"-summarization-PCA-labels.png",sep=""))
  pca.output.table.file<-file.path("Robjects",paste(SessionID,"-summarization-PCA.tab",sep=""))
  
  # output table
  pca.table <- as.data.frame(counts.pca$x[,1:4])
  rownames(pca.table)<-sample.labels
  
  # Plot expected groups of samples, when available
  if (is.factor(sample.groups)){
    
    pca.table<-cbind(pca.table,"Group"=as.vector(sample.groups))
    
    # ggbiplot can only draw ellipses if there are at least 5 samples
    if (length(rownames(counts.pca$x)) > 4){
      drawEllipses <- TRUE
    }else{
      drawEllipses <- FALSE
    }
    
    # print PCA with points
    png(file = pca.output.file,
        units="cm",
        width = 15, height = 15,
        res = 300, bg = "white")
    g <- ggbiplot(counts.pca, obs.scale = 1, var.scale = 1, var.axes = F,
                  groups = sample.groups,
                  ellipse = drawEllipses, circle = TRUE)
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal',
                   legend.position = 'top')
    if (length(levels(sample.groups))>7){
      g <- g + guides(col = guide_legend(ncol = 5, byrow = TRUE))
    }
    print(g)
    dev.off()
    # Print PCA with labels
    png(file = pca.output.lbl.file,
        units="cm",
        width = 20, height = 20,
        res = 300, bg = "white")
    g <- ggbiplot(counts.pca, obs.scale = 1, var.scale = 1, var.axes = F,
                  groups = sample.groups,
                  ellipse = drawEllipses, circle = TRUE,
                  labels = sample.labels
    )
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal',
                   legend.position = 'top')
    if (length(levels(sample.groups))>7){
      g <- g + guides(col = guide_legend(ncol = 5, byrow = TRUE))
    }
    print(g)
    dev.off()
  }else{
    # print PCA with points
    png(file = pca.output.file,
        units="cm",
        width = 15, height = 15,
        res = 300, bg = "white")
    g <- ggbiplot(counts.pca, obs.scale = 1, var.scale = 1, var.axes = F,
                  ellipse = TRUE, circle = TRUE)
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal', 
                   legend.position = 'top')
    print(g)
    dev.off()
    # Print PCA with labels
    png(file = pca.output.lbl.file,
        units="cm",
        width = 15, height = 15,
        res = 300, bg = "white")
    g <- ggbiplot(counts.pca, obs.scale = 1, var.scale = 1, var.axes = F,
                  ellipse = TRUE, circle = TRUE,
                  labels = sample.labels
    )
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal', 
                   legend.position = 'top')
    
    print(g)
    dev.off()
  }
  
  # Convert format of output figures from png to eps
  if(bmeps){
    system(paste("/usr/local/bin/bmeps -leps ",pca.output.file,sep=""))
    system(paste("/usr/local/bin/bmeps -leps ",pca.output.lbl.file,sep=""))  
  }else{
    system(paste("sam2p",pca.output.file,gsub("\\.png$",".eps",pca.output.file),sep=" "))
    system(paste("sam2p",pca.output.lbl.file,gsub("\\.png$",".eps",pca.output.lbl.file),sep=" "))
  }
  
  # export PC table
  write.table(pca.table,
              file = pca.output.table.file,
              quote=FALSE,
              sep="\t",
              dec=",",
              row.names=TRUE,
              col.names=NA)
  
  return(TRUE)
  
}
