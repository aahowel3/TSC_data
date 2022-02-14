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
#DO NOT CHANGE TO RECURSIVE=F IT WILL CHANGE INPUT ORDER OF FILES 

#exclude ones that come from subdirectory "TODELETE" 
#change "unfinished_stuff" to "to_DELETE"
#file_list <- file_list_raw[ !grepl("subdirectory", file_list_raw) ]
#omit fibor data
file_list <- file_list[ !grepl("fibroblast", file_list) ]


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


#try removing unaffected indivs - see if this affects DEs between mild and severe 
mergeded = subset(mergeded, select=-c(ND_N015P02_3_PB_Whole_FN015G2S1MP_TSRGD))


#others to drop 
#remove 39
mergeded = subset(mergeded, select=-c(ND_N039P01_1_PB_Whole_C2,
                                      ND_N039P03_1_PB_Whole_C2))


###remove 6 
mergeded = subset(mergeded, select=-c(LNTS_0049_1_PB_Whole_C1_OHFSK_L13685,
                                      LNTS_0052_1_PB_Whole_C1_OHFSK_L13686,
                                      LNTS_0055_1_PB_Whole_C1_OHFSK_L13687))
                                      
#remove 8                                     
mergeded = subset(mergeded, select=-c(LNTS_0059_1_PB_Whole_C1_OHFSK_L13688,
                                      LNTS_0061_1_PB_Whole_C1_OHFSK_L13689,
                                      LNTS_0062_1_PB_Whole_C1_OHFSK_L13690))

write.csv(mergeded, "expressiondata_minus8.csv")

#create metadata file 
#ah leave it for later quick and dirty excel test
#lapply(file_list, read.table, sep="", header = TRUE) %>%
#  set_names(file_list) %>%
#  bind_rows(.id = 'grp')
#metadata=read.csv("metadata.csv",sep=",")
#old just means I removed 6 and 39
metadata=read.csv("metadata_LNTS_preDOD_assumptions_old.csv",sep=",")

#you dont have to do this if theyre are filtered now
#try removing unaffected indivs - see if this affects DEs between mild and severe 
row_names_df_to_remove = c("8")
#remove 39 
row_names_df_to_remove = c("8","23","24")
#remove 6 and 8 and 39 and 15
row_names_df_to_remove = c("8","23","24","35","36","37","38","39","40")

#remove new 8
#try removing unaffected indivs - see if this affects DEs between mild and severe 
row_names_df_to_remove = c("32","33","34")
metadata=metadata[!(row.names(metadata) %in% row_names_df_to_remove),]

#metadata['condition:batch'] <- paste(metadata$condition, ":", metadata$batch)


#try putting it into deseq
#will not run without metadata file 
#will need to account for ~family and ~treatment affect 
dds <- DESeqDataSetFromMatrix(countData=mergeded, 
                              colData=metadata, 
                              design=~condition + family, tidy = TRUE)


#exploratory analysis and visualization
#drop rows where everything is zero 
nrow(dds)
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
nrow(dds)


#even more stringent filterirng 
# at least 3 samples with a count of 10 or higher (not just the row total)
keep <- rowSums(counts(dds) >= 10) >= 3
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
rownames(sampleDistMatrix) <- paste( vsd$family, vsd$condition, sep = " - " )
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
names = colData(dds)$family

#try rlog and vsd data
pcaData <- plotPCA(vsd, intgroup = c("family","condition"), returnData = TRUE)
percentVar <- round(100 * attr(vsd, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", "43%", "% variance")) +
  ylab(paste0("PC2: ", "30%", "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data") + 
  geom_text_repel(aes(label=names))


#no idea why the above code is printing NA variance but if you use this itll print it auto for you
#and you can just manually plug in those numbers for the axis 
plotPCA(vsd, intgroup="condition")


rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(metadata, pca$x)
ggplot(df, aes(PC12, PC13, label=id_short)) + geom_point() +
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
res <- results(ddx, contrast=c("condition","severe","mild"),alpha = 0.05) 

#independent filtering
#alpha is what you can tweak to change sig level
res <- results(ddx, contrast=c("condition","severe","mild"),alpha = 0.05,lfcThreshold = 0) 
               
#print out annotated names of the DE genes between each treatment + control 
#specify what constrasts you want to look at 
res <- results(ddx, contrast=c("condition","unaffected","severe"))
summary(res)
resSig=res


#resSig <- subset(res, padj < 0.1)
#summary(resSig)

####other combinations/corrections 
#print out annotated names of the DE genes between each treatment + control 
#specify what constrasts you want to look at 
res <- results(ddx, contrast=c("condition","mild","unaffected"),alpha = 0.3)
resSig=res
resSig <- subset(res, padj < 0.1)


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
#RUN THIS FOR EACH RESSIG COMBINATION - will work on getting rid of the NAs later 
resSig=res

ens.str <- substr(rownames(resSig), 1, 15)
resSig$symbol <- mapIds(org.Hs.eg.db,
                        keys=ens.str,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
resSig$entrez <- mapIds(org.Hs.eg.db,
                        keys=ens.str,
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

#checking CKD outputs - have to do this for each cobbination set 
resOrdered <- resSig[order(resSig$padj),]
resOrderedDF <- as.data.frame(resOrdered)
resOrderedDF[resOrderedDF$symbol %like% "CDK", ]

df=resOrderedDF[resOrderedDF$symbol %like% "CDK", ]


resOrdered <- resSig[order(resSig$pvalue),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results.csv")
write.csv(resOrderedDF, file = "results2.csv")

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


ggeaout=gsRanking(air.res)
write.csv(ggeaout, file = "out.csv")
eaBrowse(air.res)


gsea.all <- sbea(method="gsea", se= try, gs=kegg.gs, perm=100)  

library(arsenal)
results_grouped=read.csv("DEresults.csv")
intersect(results_grouped$ï..sev_mild_BH_0.1_nolifc.ENSG, results_grouped$sev_unaff_BH_0.1_nolifc.ENSG)
intersect(results_grouped$ï..sev_mild_BH_0.1_nolifc.ENSG, results_grouped$mild_unaff_BH_0.5_nolifc.ENSG)
intersect(results_grouped$sev_unaff_BH_0.1_nolifc.ENSG, results_grouped$mild_unaff_BH_0.5_nolifc.ENSG)

#check interaction between all 3 sets 
intersect(intersect(results_grouped$ï..sev_mild_BH_0.1_nolifc.ENSG, results_grouped$sev_unaff_BH_0.1_nolifc.ENSG),results_grouped$mild_unaff_BH_0.5_nolifc.ENSG)

results_grouped=read.csv("osa_gsearesults.csv")
intersect(results_grouped$sev_mild_BH_0.1_nolifc.ova.go, results_grouped$sev_mild_BH_0.1_nolifc.gsea.go)

#look at intersection between 
ens.str=rownames(vsd)
resSig=as.data.frame(ens.str)
#^have to run this through annote first 
results_grouped=read.csv("DEresults.csv")


mtorgenes=read.csv("mTor_gene_list.csv")
length(intersect(mtorgenes$ï..Gene..153.genes., resSig$symbol))


intersect(results_grouped$sev_mild_BH_0.1_nolifc.Ensembl, mtorgenes$ï..Gene..153.genes.)
intersect(results_grouped$sev_unaff_BH_0.1_nolifc.Ensembl, mtorgenes$ï..Gene..153.genes.)
intersect(results_grouped$mild_unaff_BH_0.1_nolifc.Ensembl, mtorgenes$ï..Gene..153.genes.)

#just quick check all DE genes between fam and child and check if any of those appear in mtor set
#if any are wil figure out which parent-child it came from later 
within=read.csv("within-family-analysis.csv")
ens.str=within$all
intersect(within$symbol, mtorgenes$ï..Gene..153.genes.)


#hm check if any DE genes between parent-child are DE between severe-mild-unaffected
intersect(results_grouped$sev_mild_BH_0.1_nolifc.Ensembl, within$symbol)
intersect(results_grouped$sev_unaff_BH_0.1_nolifc.Ensembl, within$symbol)
intersect(results_grouped$mild_unaff_BH_0.1_nolifc.Ensembl, within$symbol)

overlap=read.csv("overlapping.csv",header = FALSE)
ens.str=overlap$V1

overlap$symbol <- mapIds(org.Hs.eg.db,
                        keys=ens.str,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
overlap$entrez <- mapIds(org.Hs.eg.db,
                        keys=ens.str,
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

write.csv(overlap, file = "overlap_withensembl.csv")



counts(ddx, normalized=T)

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
write.csv(ggeaout, file = "ggeaout.csv")

#graph signifigant results - use diff tool
par(mfrow=c(1,2))
ggeaGraph( gs=hsa.gs[["hsa04390_Hippo_signaling_pathway"]], grn=hsa.grn, se=try)
ggeaGraphLegend()   




results1=read.csv("results.csv")
results2=read.csv("results2.csv")

intersect(results1$X, results2$X)


###########################
#########################################rerun with batch effects s
#check for batch effects 
#https://biodatascience.github.io/compbio/dist/sva.html 
library("sva")
dat  <- counts(ddx, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition, colData(ddx))
mod0 <- model.matrix(~   1, colData(ddx))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
norm.cts <- dat[rowSums(dat) > 0,]
fit <- svaseq(norm.cts, mod=mod, mod0=mod0, n.sv=2)


par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ ddx$family, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

#i like this graphial representation better
library(rafalib)
bigpar()
plot(fit$sv[,1:2], col=dds$family, cex=2,
     xlab="SV1", ylab="SV2")
legend("top", levels(dds$family), pch=16,
       col=1:3, cex=.8, ncol=3, title="family")

options(ggrepel.max.overlaps = Inf)
ggplot(as.data.frame(fit$sv[,1:2]), aes(x = V1, y = V2, color = dds$batch, shape = dds$condition)) +
  geom_point(size =3) +
  xlab(paste0("SV1")) +
  ylab(paste0("SV2")) +
  coord_fixed() +
  geom_text_repel(aes(label=dds$id_short))


#rerunning deseq with additional surrgate variables to account for batch effects
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]

design(ddssva) <- ~ SV1 + condition + family 
#design(ddssva) <- ~ condition + SV1 



ddssva2 = DESeq(ddssva)

vsd <- vst(ddssva, blind = FALSE)


#PCA the surrgate variable ones 
#try rlog and vsd data
pcaData <- plotPCA(vsd, intgroup = c("family","condition"), returnData = TRUE)
percentVar <- round(100 * attr(vsd, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = family, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data") + 
  geom_text_repel(aes(label=names))


#specify what constrasts you want to look at 
res <- results(ddssva2, contrast=c("condition","severe","mild"),alpha=0.1)
res
summary(res)
#run with BH adjusted p value
resSig <- subset(res, padj < 0.05)
summary(resSig)


#BiocManager::install("biomaRt")
library(biomaRt)
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
filterType <- "ensembl_gene_id"
filterValues <- rownames(resSig)[1:1000]

attributeNames <- c('ensembl_gene_id','external_gene_name','entrezgene_accession')

annot <- getBM(attributes=attributeNames, 
               filters = filterType, 
               values = filterValues, 
               mart = ensembl)



#ANNOTATING AND EXPORTING RESULTS
library("AnnotationDbi")
library("org.Hs.eg.db")

ens.str <- substr(rownames(resSig), 1, 15)
resSig$symbol <- mapIds(org.Hs.eg.db,
                        keys=ens.str,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
resSig$entrez <- mapIds(org.Hs.eg.db,
                        keys=ens.str,
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

resOrdered <- resSig[order(resSig$pvalue),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results.csv")




#does not work 
#enrichment browser only allows for 2 fixed vaariables at a time - condition + group 
#and by adding the SVA variableds in - +condition + group + SV1 + SV2 
#will not take it 
#but i dont know why enrichment browser needs the design of the study since it only looks at genes
#in DE and not in DE set compared to those in gene set and those DE in gene set 
#se <- import(ddssva2, res)

#have to use this one with new res object but old deseq object 
se <- import(ddx, res)


#remaining steps in the enrichmentbrowser are: 
#looks like idMap works just as well as AnnotationDBI 
try <- idMap(se, org = "hsa", from = "ENSEMBL", to = "ENTREZID")
#choose a gene set - could be kegg, GO, etc 
kegg.gs <- getGenesets(org = "hsa", db = "kegg")
go.gs <- getGenesets(org = "hsa", db = "go")
#perform which set-based enrichment analysis (sbea) you choose - ora = overrep analysis
air.res <- sbea(method = "ora", se = try, gs = kegg.gs, perm = 0, alpha = 0.2)

ggeaout=gsRanking(air.res)
write.csv(ggeaout, file = "out3.csv")
#try with GO too see results diff
######try increasing memory 
air.res <- sbea(method = "ora", se = try, gs = go.gs, perm = 0, alpha = 0.2)
ggeaout=gsRanking(air.res)
write.csv(ggeaout, file = "out4.csv")




mtorgenes=read.csv("mTor_gene_list.csv")
intersect(resSig$symbol, mtorgenes$ï..Gene..153.genes.)


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
write.csv(ggeaout, file = "ggeaout.csv")

