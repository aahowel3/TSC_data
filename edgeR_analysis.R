library(edgeR)
library(EnrichmentBrowser)
library("data.table")
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
file_list <- file_list[ !grepl("fibroblast", file_list) ]

#only do per-family analysis13
file_list <- file_list[ grepl("N036", file_list) ]

#N036 quickly drop sibling to look at just child-father
file_list <- file_list[ !grepl("P04", file_list) ]

#WONT WORK W PREDOD BC ALL FILES NAMED DIFFERENTLY
#family 
#file_list <- file_list[ grepl("LNTS_0001|LNTS_0002|LNTS_0003", file_list) ]
#family5
file_list <- file_list[ grepl("LNTS_0004|LNTS_0005", file_list) ]
#file_list <- file_list[ grepl("LNTS_0010|LNTS_0011", file_list) ]
#family 3
#file_list <- file_list[ grepl("LNTS_0041|LNTS_0043|LNTS_0045", file_list) ]
#family 6
#file_list <- file_list[ grepl("LNTS_0049|LNTS_0052|LNTS_0055", file_list) ]
#family 8 
#file_list <- file_list[ grepl("LNTS_0059|LNTS_0061|LNTS_0062", file_list) ]
#file_list <- file_list[ !grepl("fibroblast", file_list) ]
#file_list <- file_list[ grepl("36", file_list) ]


#cant just remove the unaffected from the analysis - changes the DE
#sure you can
#file_list <- file_list[ !grepl("N015P02_", file_list) ]

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
metadata=read.csv("metadata_LNTS_preDOD_assumptions_old.csv",sep=",")
#only do per-family analysis
#CHANGE PER FAMILY 
metadata <- metadata[ (metadata$family=="N036"), ]

#drop N036 sibling
row_names_df_to_remove = c("21")
metadata=metadata[!(row.names(metadata) %in% row_names_df_to_remove),]


#drop unaffected
#row_names_df_to_remove = c("8")
#metadata=metadata[!(row.names(metadata) %in% row_names_df_to_remove),]

#fibro meta with edgeR
#metadata=read.csv("fibroblast_metadata.csv",sep=",")
#metadata <- metadata[ (metadata$family=="N036"), ]


#row.names.remove <- c("8")
#metadata=metadata[!(row.names(metadata) %in% row.names.remove), ]

#most N0 fams 
groups<- c( "severe","mild")

#n036
groups<- c( "severe","mild","mild")


#family2
#groups<- c( "mild","severe","mild")
#fam5 
groups<- c("mild","severe")
#family 3 
#fam5 
#groups<- c("mild","severe","severe")
#fam6
#groups<- c("mild","mild","severe")

#fam8 
#groups<- c( "mild","severe","mild")



gene_names = c(subset(mergeded, select = c(V1) ))
mergeded = subset(mergeded, select = -c(V1) )


dgList =DGEList(counts=mergeded,group=factor(groups),genes=gene_names)

#filter out low count genes
countsPerMillion <- cpm(dgList)
summary(countsPerMillion)
#'summary' is a useful function for exploring numeric data; eg. summary(1:100)
countCheck <- countsPerMillion > 1
head(countCheck)

#this does not mean #transcripts > 2 it means log transformed cpm >2 - reasonable filter 
keep <- which(rowSums(countCheck) >= 2)


dgList <- dgList[keep,]
summary(cpm(dgList)) 

#perform normalization
dgList <- calcNormFactors(dgList, method="TMM")

#set up model desing 
#taken from more complicated tutorial
designMat <- model.matrix(~groups)
designMat


#estimate disperions 
#will throw a warning when you onyl have 1 sample per condition 
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
#dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
#dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)

#the function above estmates disperesion parameters - which you cant do w one sample
#so you can instead just make-up a BCV 
bcv <- 0.35

#lets keep all comparisons as GLM not exact test
#et <- exactTest(dgList, dispersion=bcv^2)
#et_results=topTags(et)
#summary(decideTests(et))

#will consistently use glm test 
fit <- glmFit(dgList, dispersion=bcv^2,designMat)
lrt <- glmLRT(fit,coef=ncol(fit$design))
summary(decideTests(lrt))

#if more than 2 conditions use these comparisons
#################################################
#3 v 1 
#lrt <- glmLRT(fit, coef=3)
#summary(decideTests(qlf.3vs1))
#2 v. 3  contrast=c(0,-1,1)
#lrt <- glmLRT(fit, contrast=c(0,-1,1))
#summary(decideTests(qlf.2v3))
#####################################################

#this should be the same number as summary(decidetestDEG) 
out <- topTags(lrt, n=Inf, adjust.method="BH")
res <- topTags(lrt, n = nrow(dgList), sort.by = "none")

keep <- out$table$FDR <= 0.05 & abs(out$table$logFC) >= 1
keep2=out[keep,]
keep2=as.data.frame(keep2)

#RUN THIS FOR EACH RESSIG COMBINATION - will work on getting rid of the NAs later 
ens.str <- keep2$V1

library("AnnotationDbi")
library("org.Hs.eg.db")
keep2$symbol <- mapIds(org.Hs.eg.db,
                        keys=ens.str,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
keep2$entrez <- mapIds(org.Hs.eg.db,
                        keys=ens.str,
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

resOrderedDF <- as.data.frame(keep2)
write.csv(resOrderedDF, file = "results2.csv")



#dont really need to run this part if you use upsetR_plots.R
############################
#manually putting it all togeher 
results_grouped=read.csv("within-family-analysis.csv")
df=as.data.frame(results_grouped)
genes=as.list(unique(df$all))

cat(unique(df$all),sep="\n")

for (i in genes){
  cols <- colSums(mapply('==', i, df))
  new.df <- df[,which(cols > 0)]
  print(colnames(new.df))
} 


venn=read.csv("overlapping.csv", header = FALSE)
#######################################
#check none of the 1:1 family genes contained Cdk
ens.str <- df$all

df$symbol <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
df$entrez <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

df[df$symbol %like% "CDK", ]
###############useless

###still do this for each - less overlap easier to interpert with hack method 
#for each run lets do an ORA and a GGEA
go <- goana(lrt)
out=topGO(go)
frame=as.data.frame(out$Term)

keg <- kegga(lrt)
out=topKEGG(keg)
frame$kegg=out$Pathway
write.csv(frame, file = "results3.csv")








###once its all manually complied look here 
data=read.csv("go_overlap_withinfamily_correct.csv",header=FALSE)

df=as.data.frame(data)
genes=as.list(unique(df$V1))

cat(unique(df$V1),sep="\n")

for (i in genes){
  cols <- colSums(mapply('==', i, df))
  new.df <- df[,which(cols > 0)]
  print(colnames(new.df))
} 

###once its all manually complied look here 
data=read.csv("kegg_overlap_withinfamily_correct.csv",header=FALSE)

df=as.data.frame(data)
genes=as.list(unique(df$V1))

cat(unique(df$V1),sep="\n")

for (i in genes){
  cols <- colSums(mapply('==', i, df))
  new.df <- df[,which(cols > 0)]
  print(colnames(new.df))
} 

