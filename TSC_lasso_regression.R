install.packages("caret", dependencies=T)
install.packages(c("randomForest", "glmnet"), dependencies=T)


require(caret)
require(randomForest)
require(glmnet)

fullFeatureSet <- read.table("http://seandavi.github.io/ITR/expression-prediction/features.txt");
target <- scan(url("http://seandavi.github.io/ITR/expression-prediction/target.txt"), skip=1)


#https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("sva")
#https://www.biostars.org/p/452858/

library(EnrichmentBrowser)
library("data.table")
library(DESeq2)
library("tidyverse")
library(statmod)
library(ggrepel)


setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata")

#grab a list of files that end with .htseqcounts
#in all subdirectories - so run this script at the top of results
#add in *PB* when running on the cluster
file_list <- list.files(pattern="*.htSeqCounts", recursive = T)

#exclude ones that come from subdirectory "TODELETE" 
#change "unfinished_stuff" to "to_DELETE"
file_list <- file_list[ !grepl("subdirectory", file_list) ]
#omit fibor data
file_list <- file_list[ !grepl("fibroblast", file_list) ]
file_list <- file_list[ !grepl("N029", file_list) ]
file_list <- file_list[ !grepl("N036", file_list) ]


#convert the second column to the filename itself 
imported_files <- lapply(file_list, function(x) {
  DT <- fread(x)
  new_colname <- x %>%
    basename %>%
    str_replace_all(c("^output_"="", "\\.htSeqCounts$"="")) %>%
    str_c("", .)
  setnames(DT, old="V2", new=new_colname)
  return(DT)
})

#merge and reduce the individual dataframes 
merged_data <- reduce(imported_files, merge, by="V1", all=TRUE)
#convert to dataframe to fix rownames:NULL error after DESeqDataSetFromMatrix
mergeded=as.data.frame(merged_data)


#create metadata file 
#ah leave it for later quick and dirty excel test
#lapply(file_list, read.table, sep="", header = TRUE) %>%
#  set_names(file_list) %>%
#  bind_rows(.id = 'grp')
#metadata=read.csv("metadata.csv",sep=",")
metadata=read.csv("metadata_extended.csv",sep=",")

#try putting it into deseq
#will not run without metadata file 
#will need to account for ~family and ~treatment affect 
dds <- DESeqDataSetFromMatrix(countData=mergeded, 
                              colData=metadata, 
                              design=~condition + family + TSC1, tidy = TRUE)


#exploratory analysis and visualization
#drop rows where everything is zero 
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)



#testing different transformations for the data
#still recccomended to do deseqs interanl normalization tho
#vst transformation 
vsd <- vst(dds, blind = FALSE)