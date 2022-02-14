library(EnrichmentBrowser)
library("data.table")
library(DESeq2)
library("tidyverse")
library(statmod)
library(ggrepel)

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata")

mergeded=read.csv("expressiondata_minus8.csv",sep = ",")
metadata=read.csv("metadata_LNTS_preDOD_assumptions_minus8.csv",sep=",")

#try putting it into deseq
#will not run without metadata file 
#will need to account for ~family and ~treatment affect 
dds <- DESeqDataSetFromMatrix(countData=mergeded, 
                              colData=metadata, 
                              design=~condition + family, tidy = TRUE)


#even more stringent filterirng 
# at least 3 samples with a count of 10 or higher (not just the row total)
#keep <- rowSums(counts(dds) >= 10) >= 3
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]
nrow(dds)

#actually running the differnetial expression 
ddx = DESeq(dds)


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


#rerunning deseq with additional surrgate variables to account for batch effects
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]

design(ddssva) <- ~ SV1 + SV2 + condition + family 
#design(ddssva) <- ~ condition + SV1 


ddssva2 = DESeq(ddssva)


#specify what constrasts you want to look at 
res <- results(ddssva2, contrast=c("condition","severe","mild"),alpha=0.05)
summary(res)

#############################################################plotting DE genes as a heatmap 
deseq2Results <- res
deseq2ResDF <- as.data.frame(deseq2Results)
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

#deseq2VST <- vst(ddssva)
#trying diff transforms - not helpful
#deseq2VST <- normTransform(ddssva)
#deseq2VST <- rlog(ddssva)
#deseq2VST=deseq2ResDF
deseq2VST <- (ddssva)


# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)

#reorder column names based on order in list 
level_order=metadata[order(metadata$condition),]
deseq2VST=deseq2VST[,level_order$id]
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)



# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 0,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]



# Convert the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
#here is where the column "value " is automatically introduced with the shifting 
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

library(scales) # needed for oob parameter
library(viridis)
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

###WOOF heatmap representation is aboslutely terrible - really highlights how little of a difference there was 
#try volcano plot instead 
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlim=c(-3,3)))

# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
library(calibrate)
#with(subset(res, padj<.05), textxy(log2FoldChange, -log10(pvalue), labs=symbol, pos=4, cex=.8))


library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=log2FoldChange, label=symboll)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) 

















#tried again with relative expression values - still not any better 

library(ComplexHeatmap)
library(circlize)

#####after these steps you ahve the same format as example expr 
#deseq2VST <- (ddssva)
# Convert the DESeq transformed object to a data frame
#deseq2VST <- assay(deseq2VST)
#deseq2VST <- as.data.frame(deseq2VST)
####need to get it down to this step
#deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]
#to reduce number of genes it has to plot or itll overload R

#https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html 

mat = as.matrix(deseq2VST[, grep("N", colnames(deseq2VST))])
mat_scaled = t(apply(mat, 1, scale))

colnames(mat_scaled) = level_order$id

library(ComplexHeatmap)
library(circlize)


ht=Heatmap(mat_scaled, name = "expression", row_km = 0, 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
           column_order = level_order$id,
           
                 row_title = NULL, show_row_dend = FALSE,show_column_dend = FALSE)
  
draw(ht, padding = unit(c(100, 1, 1, 1), "mm"))

