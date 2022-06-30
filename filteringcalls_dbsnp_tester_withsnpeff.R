#it DOES have the chromosme and pos info its just hiding it - not under info(vcf) but rowranges(vcf)
#https://rdrr.io/bioc/VariantAnnotation/man/readVcf-methods.html
#setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata/")
library(VariantAnnotation)
library(biomaRt)
library(dplyr)
library(tidyr)
library(tidyverse)

#combine these two
vcf <- readVcf("N013_annotated.snpeff.snpsift.dbnsfp4.1_withnarg_1000.vcf", "hg19")
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

total$dbNSFP
total$dbNSFP_gnomAD_genomes_EAS_AF

##When the qualitative prediction tools (except for MutationAssessor and Aloft) had D the variant gained one score
total$score_qual = rowSums(sapply(total[c("dbNSFP_SIFT_pred",
                                          "dbNSFP_SIFT4G_pred",
                                          "dbNSFP_Polyphen2_HDIV_pred",
                                          "dbNSFP_Polyphen2_HVAR_pred",
                                          "dbNSFP_LRT_pred",
                                          "dbNSFP_MutationTaster_pred",
                                          "dbNSFP_FATHMM_pred",
                                          "dbNSFP_PROVEAN_pred",
                                          "dbNSFP_MetaSVM_pred",
                                          "dbNSFP_MetaLR_pred",
                                          "dbNSFP_M_CAP_score",
                                          "dbNSFP_PrimateAI_pred",
                                          "dbNSFP_DEOGEN2_pred",
                                          "dbNSFP_BayesDel_addAF_pred",
                                          "dbNSFP_BayesDel_noAF_pred",
                                          "dbNSFP_ClinPred_pred",
                                          "dbNSFP_LIST_S2_pred",
                                          "dbNSFP_fathmm_MKL_coding_pred",
                                          "dbNSFP_fathmm_XF_coding_pred"
                                          
                                          )], grepl, pattern = "D")) 
total=total %>% mutate_at(vars("score_qual"), ~replace_na(.,0))

#An indicator for a deleterious variant of MutationAssessor was "H" and of Aloft was "R" or "D" with high confidence.
#MutationAssessor
#Aloft
#can only use rowsums if greater than one field
#Some ALoFT scores/information are missing in dbNSFP
total$score_mut = ifelse(total$dbNSFP_MutationAssessor_pred == "H", 1, 0)
total=total %>% mutate_at(vars("score_mut"), ~replace_na(.,0))

total$score_mut2 = ifelse(total$dbNSFP_Aloft_pred == "R" | total$dbNSFP_Aloft_pred == "D", 1, 0)
total=total %>% mutate_at(vars("score_mut2"), ~replace_na(.,0))


insert_zero <- function(smpl_vec){
  smpl_vec[is.na(smpl_vec)] <- 0
  baseline=quantile(as.numeric(unlist(smpl_vec)), probs = c(.90),na.rm=TRUE)
  baseline_num=unname(baseline)
  totals=ifelse(grepl("NA",smpl_vec),0,smpl_vec)
  score_quant = ifelse(as.numeric(unlist(totals)) > baseline_num, 1, 0)
  return(score_quant)
}

total$score_quant1=apply(total[c('dbNSFP_VEST4_score')], 2, insert_zero)
total$score_quant2=apply(total[c('dbNSFP_REVEL_score')], 2, insert_zero)
total$score_quant3=apply(total[c('dbNSFP_MutPred_score')], 2, insert_zero)
total$score_quant4=apply(total[c('dbNSFP_MVP_score')], 2, insert_zero)
total$score_quant5=apply(total[c('dbNSFP_MPC_score')], 2, insert_zero)
total$score_quant6=apply(total[c('dbNSFP_DANN_score')], 2, insert_zero)
total$score_quant7=apply(total[c('dbNSFP_CADD_phred')], 2, insert_zero)
total$score_quant8=apply(total[c('dbNSFP_Eigen_phred_coding')], 2, insert_zero)
total$score_quant9=apply(total[c('dbNSFP_Eigen_PC_phred_coding')], 2, insert_zero)


#total=total %>% mutate_at(vars("score_quant"), ~replace_na(.,0))
total$genename = total$GENEINFO$X1 
total$score_final = total$score_qual + total$score_mut + total$score_mut2 + total$score_quant1 + total$score_quant2 +
  total$score_quant3 + total$score_quant4 +total$score_quant5 + total$score_quant6 + total$score_quant7 + total$score_quant8 + total$score_quant9 
final_table = total[c("genename", "seqnames", "start", "Row.names", "ND_N013P01", "ND_N013P02","score_final")]

ordered_df <- final_table[order(-final_table$score_final),]
ordered_df <- apply(ordered_df,2,as.character)

write.csv(ordered_df, file="dbnsfp_allN013_annotations_snp.csv")




