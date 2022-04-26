#it DOES have the chromosme and pos info its just hiding it - not under info(vcf) but rowranges(vcf)
#https://rdrr.io/bioc/VariantAnnotation/man/readVcf-methods.html
#setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata/")
library(VariantAnnotation)
library(biomaRt)
library(dplyr)
library(tidyr)

library(foreach)
library(doParallel)
library(iterators)
#ensembl = useEnsembl(biomart='ensembl', 
#                     dataset="hsapiens_gene_ensembl") 

#gr37 seems the better one to use
ensembl = useMart(biomart='ensembl', host="grch37.ensembl.org",
                     dataset="hsapiens_gene_ensembl") 

#combine these two
vcf <- readVcf("N013_annotated.vcf", "hg19")
df = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))

total <- merge(df,df2, by="row.names")

#https://support.bioconductor.org/p/127035/
product = function(x, output){
  
  # accessing elements from first column
  chromosome = x$seqnames
  
  # accessing elements from second column
  start = x$start
  
  # accessing elements from third column
  end = x$end
  
  positions <- data.frame(chromosome,
                          start,
                          end)
  
  results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                   filters = c("chromosome_name", "start", "end"),
                   values = list(positions[,1], positions[,2], positions[,3]),
                   mart = ensembl)
  
  # return product
  return(results$hgnc_symbol)
}

numCores <- detectCores()
numCores
#registerDoParallel(numCores)
registerDoParallel(10)

total$genename = foreach(row = iter(total, by = 'row'), .packages='biomaRt') %dopar% {
  product(row, 0)
}

mtorgenes=read.csv("mTor_gene_list.csv")
total$mtor=ifelse(total$genename %in% mtorgenes$gene, 1, 0)

total=dplyr::filter(total, total$mtor=="1")

##When the qualitative prediction tools (except for MutationAssessor and Aloft) had D the variant gained one score
total$score_qual = rowSums(sapply(total[c("dbNSFP_SIFT_pred",
                                        "dbNSFP_Polyphen2_HDIV_pred",
                                        "dbNSFP_Polyphen2_HVAR_pred",
                                        "dbNSFP_LRT_pred",
                                        "dbNSFP_MutationTaster_pred",
                                        "dbNSFP_FATHMM_pred",
                                        "dbNSFP_PROVEAN_pred",
                                        "dbNSFP_MetaSVM_pred")], grepl, pattern = "D")) 
total=total %>% mutate_at(vars("score_qual"), ~replace_na(.,0))

#An indicator for a deleterious variant of MutationAssessor was "H" and of Aloft was "R" or "D" with high confidence.
#MutationAssessor
#Aloft
#can only use rowsums if greater than one field
#Some ALoFT scores/information are missing in dbNSFP
total$score_mut = ifelse(total$dbNSFP_MutationAssessor_pred == "H", 1, 0)
total=total %>% mutate_at(vars("score_mut"), ~replace_na(.,0))


#baseline 
baseline=quantile(as.numeric(total$dbNSFP_CADD_phred), probs = c(.90),na.rm=TRUE)
#variant indicators were higher than 90% of predicted variants
total$score_quant = ifelse(total$dbNSFP_CADD_phred > baseline, 1, 0)
total=total %>% mutate_at(vars("score_quant"), ~replace_na(.,0))



total$score_final = total$score_qual + total$score_mut + total$score_quant
final_table = total[c("genename", "score_final")]

ordered_df <- final_table[order(-final_table$score_final),]
ordered_df <- apply(ordered_df,2,as.character)

write.csv(ordered_df, file="dbnsfp_allN013_annotations.csv")





