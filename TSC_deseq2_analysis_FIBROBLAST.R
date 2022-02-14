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


setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata/fibroblast/")

#grab a list of files that end with .htseqcounts
#in all subdirectories - so run this script at the top of results
#add in *PB* when running on the cluster
file_list <- list.files(pattern="*.htSeqCounts", recursive = T)

#exclude ones that come from subdirectory "TODELETE" 
#change "unfinished_stuff" to "to_DELETE"
#file_list <- file_list_raw[ !grepl("subdirectory", file_list_raw) ]

file_list <- file_list[ !grepl("N030", file_list) ]


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

write.csv(mergeded, "fibro_exp_out.csv")
#create metadata file 
#ah leave it for later quick and dirty excel test
#lapply(file_list, read.table, sep="", header = TRUE) %>%
#  set_names(file_list) %>%
#  bind_rows(.id = 'grp')
#metadata=read.csv("metadata.csv",sep=",")
metadata=read.csv("fibroblast_metadata_updated_2_10_21.csv",sep=",")

#try putting it into deseq
#will not run without metadata file 
#will need to account for ~family and ~treatment affect 
dds <- DESeqDataSetFromMatrix(countData=mergeded, 
                              colData=metadata, 
                              design=~condition, tidy = TRUE)


#exploratory analysis and visualization
#drop rows where everything is < 10 - keep in line with blood plasma filtering 
nrow(dds)
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
nrow(dds)


#testing different transformations for the data
#still recccomended to do deseqs interanl normalization tho
#vst transformation 
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

#plot those different transformations
dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


#checking sample disntances 
sampleDists <- dist(t(assay(vsd)))
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$ï..id_short, vsd$condition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
#checking sample disntances 
sampleDists <- dist(t(assay(rld)))
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$family, rld$condition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


###USE THE DESEQSETFROMMATRIX OBJECT TO PLUG INTO VST NOT THE DDX FROM DESEQ OBJECT
names = colData(dds)$ï..id_short
#plotPCA(vsdata, intgroup=c("family","condition")) + geom_text_repel(aes(label=names))

#try rlog and vsd data
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(vsd, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", 37, "% variance")) +
  ylab(paste0("PC2: ", 21, "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data") + 
  geom_text_repel(aes(label=names))


#no idea why the above code is printing NA variance but if you use this itll print it auto for you
#and you can just manually plug in those numbers for the axis 
plotPCA(vsd, intgroup="condition")

##trying different PC's was worth a shot
vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(metadata, pca$x)
ggplot(df, aes(PC2, PC7, label=ï..id_short)) + geom_point() +
  geom_text(aes(colour = factor(condition)))

#actually running the differnetial expression 
ddx = DESeq(dds)
#default results summary is last 2 variaibles (in this case - family)
res_basic=results(ddx)
summary(res_basic)
#CAN FURHTER reduce junk by filtering on lfcThreshold = 1 - this means instead of a log fold change just greater 
#than ANYTHIGN (0) it has to be greater than 1 
#res <- results(ddx, contrast=c("condition","severe","mild"),pAdjustMethod = "BH", lfcThreshold = 1)
#summary (res)


#print out annotated names of the DE genes between each treatment + control 
#specify what constrasts you want to look at 
res <- results(ddx, contrast=c("condition","control","severe"))
summary(res)
resSig <- subset(res, padj < 0.1)
#summary(resSig)

####other combinations/corrections 
#print out annotated names of the DE genes between each treatment + control 
#specify what constrasts you want to look at 
res <- results(ddx, contrast=c("condition","mild","unaffected"),alpha = 0.3)
resSig=res
resSig <- subset(res, padj < 0.2)


resOrdered <- res[order(res$padj),]

summary(res)
#multiple corrections with increased lfc fold change 
#than ANYTHIGN (0) it has to be greater than 1 
res <- results(ddx, contrast=c("condition","mild","unaffected"),pAdjustMethod = "BH", lfcThreshold = 1)
summary (res)
#no multiple corrections 
#wont let you change the reg pvalue cutoff? whatever its a bad idea not to do multiple corrections anyway
res <- results(ddx, contrast=c("condition","mild","unaffected"), pAdjustMethod = "none")
summary(res)



#specify what constrasts you want to look at 
#res <- results(ddx, contrast=c("condition","severe","unaffected"))
#summary(res)
#resSig <- subset(res, padj < 0.1)


#specify what constrasts you want to look at 
#res <- results(ddx, contrast=c("condition","mild","unaffected"))
#summary(res)
#resSig <- subset(res, padj < 0.1)


##############################################################################
#put top genes in a human-readable list for presentation
#ANNOTATING AND EXPORTING RESULTS
library("AnnotationDbi")
library("org.Hs.eg.db")
library("EnsDb.Hsapiens.v75")
#RUN THIS FOR EACH RESSIG COMBINATION - will work on getting rid of the NAs later 
ens.str <- substr(rownames(resSig), 1, 15)
resSig$symbol <- mapIds(EnsDb.Hsapiens.v75,
                        keys=ens.str,
                        column="SYMBOL",
                        keytype="GENEID",
                        multiVals="first")
resSig$entrez <- mapIds(EnsDb.Hsapiens.v75,
                        keys=ens.str,
                        column="ENTREZID",
                        keytype="GENEID",
                        multiVals="first")

resOrdered <- resSig[order(resSig$padj),]

resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "results.csv")
##############################################################################

#start functional enrichment on the severe v mild genes - use 6_30_2021.R script
se <- import(ddx, res)
#remaining steps in the enrichmentbrowser are: 
#looks like idMap works just as well as AnnotationDBI 
try <- idMap(se, org = "hsa", from = "ENSEMBL", to = "ENTREZID")
#choose a gene set - could be kegg, GO, etc 
kegg.gs <- getGenesets(org = "hsa", db = "kegg")
go.gs <- getGenesets(org = "hsa", db = "go")
#perform which set-based enrichment analysis (sbea) you choose - ora = overrep analysis
air.res <- sbea(method = "ora", se = try, gs = kegg.gs, perm = 0, alpha = 0.1)
#try with GO too see results diff
######try increasing memory 
air.res <- sbea(method = "ora", se = try, gs = go.gs, perm = 0, alpha = 0.1)

ggout=gsRanking(air.res)
ggout=as.data.frame(ggout)
write.csv(ggout, file = "ggout.csv")


memory.limit(1000)

gsea.all <- sbea(method="gsea", se= try, gs=kegg.gs, perm=100)  

#check if any genes from blood sets overlap with fibro sets 
#using same hack as with checking edgeR overlap 
library(arsenal)
results_grouped=read.csv("../DEresults.csv")
df=as.data.frame(results_grouped)
genes=as.list(unique(df$ï..all))

cat(unique(df$ï..all),sep="\n")

for (i in genes){
  cols <- colSums(mapply('==', i, df))
  new.df <- df[,which(cols > 0)]
  print(colnames(new.df))
} 


venn=read.csv("overlapping_fibrovblood.csv", header = FALSE)



intersect(results_grouped$ï..sev_mild_BH_0.1_nolifc.ENSG, results_grouped$sev_unaff_BH_0.1_nolifc.ENSG)
intersect(results_grouped$ï..sev_mild_BH_0.1_nolifc.ENSG, results_grouped$mild_unaff_BH_0.5_nolifc.ENSG)
intersect(results_grouped$sev_unaff_BH_0.1_nolifc.ENSG, results_grouped$mild_unaff_BH_0.5_nolifc.ENSG)

#check interaction between all 3 sets 
intersect(intersect(results_grouped$ï..sev_mild_BH_0.1_nolifc.ENSG, results_grouped$sev_unaff_BH_0.1_nolifc.ENSG),results_grouped$mild_unaff_BH_0.5_nolifc.ENSG)


results_grouped=read.csv("osa_gsearesults.csv")
intersect(results_grouped$sev_mild_BH_0.1_nolifc.ova.go, results_grouped$sev_mild_BH_0.1_nolifc.gsea.go)


#network analysis ggea
#frequires add gs set - basucally kegg set tho?
data.dir <- system.file("extdata", package="EnrichmentBrowser")
gmt.file <- file.path(data.dir, "hsa_kegg_gs.gmt")
hsa.gs <- getGenesets(gmt.file)
#create grn file
hsa.grn <- compileGRN(org="hsa", db="kegg",kegg.native=TRUE)
#run ggea
ggea.all <- nbea(method="ggea", se=try, gs=hsa.gs, grn=hsa.grn, alpha=0.2)
ggeaout=gsRanking(ggea.all)
as.data.frame(ggeaout)



#graph signifigant results - use diff tool
par(mfrow=c(1,2))
ggeaGraph( gs=hsa.gs[["hsa04514_Cell_adhesion_molecules_(CAMs)"]], grn=hsa.grn, se=try)
ggeaGraphLegend()   