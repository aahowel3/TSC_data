#it DOES have the chromosme and pos info its just hiding it - not under info(vcf) but rowranges(vcf)
#https://rdrr.io/bioc/VariantAnnotation/man/readVcf-methods.html
#setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata/")
library(VariantAnnotation)
library(biomaRt)
library(dplyr)
library(tidyr)
library(tidyverse)

#combine these two
vcf <- readVcf("N013_annotated.snpeff.snpsift.vcf", "hg19")
df = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df3=as.data.frame((geno(vcf)$GT))

total1 <- merge(df,df2, by="row.names")
total1 = total1 %>% remove_rownames %>% column_to_rownames(var="Row.names")

total=merge(total1,df3,by="row.names")

total=within(total, GENEINFO<-data.frame(do.call('rbind', strsplit(as.character(GENEINFO), ":", fixed=TRUE))))


mtorgenes=read.csv("mTor_gene_list.csv")
total$mtor=ifelse(total$GENEINFO$X1 %in% mtorgenes$gene, 1, 0)

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
baseline=quantile(as.numeric(unlist(total$dbNSFP_CADD_phred)), probs = c(.90),na.rm=TRUE)
#variant indicators were higher than 90% of predicted variants
baseline_num=unname(baseline)
total$dbNSFP_CADD_phred2=ifelse(grepl("NA",total$dbNSFP_CADD_phred),0,total$dbNSFP_CADD_phred)
total$score_quant = ifelse(as.numeric(unlist(total$dbNSFP_CADD_phred2)) > baseline_num, 1, 0)

#total=total %>% mutate_at(vars("score_quant"), ~replace_na(.,0))


total$genename = total$GENEINFO$X1 
total$score_final = total$score_qual + total$score_mut + total$score_quant
final_table = total[c("genename", "seqnames", "start", "Row.names", "ND_N013P01", "ND_N013P02","score_final")]

ordered_df <- final_table[order(-final_table$score_final),]
ordered_df <- apply(ordered_df,2,as.character)

write.csv(ordered_df, file="dbnsfp_allN013_annotations_snp.csv")





