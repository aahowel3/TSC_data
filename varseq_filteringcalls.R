#filter varseq calls 
#if column 
#this is to access the Varseq calls - they have diff col than TIGR and shouldnt be used 
#even if it supposed to be more up to date 

#just drop RNA calls pull them back for the mtor pathway variants bc Im fucking sick of this 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata/tygir_calls/")
#use LNTS and ND_Haplotypecaller datasets
#cols that dont match are 30,36,29
N030_tygir=read.csv("ND_N036_HaplotypeCaller.csv")
mtorgenes=read.csv("../mTor_gene_list.csv")


#toggle == and != to look for within mtor pathway vs outside mtor pathway
#not all have this use ifelse apply filter
#N030_tygir_mtor=filter(N030_tygir, mTor.pathway.gene. == "#N/A")

N030_tygir$mtor=ifelse(N030_tygir$Gene..HUGO. %in% mtorgenes$ï..Gene..153.genes., 1, 0)
N030_tygir=dplyr::filter(N030_tygir, N030_tygir$mtor=="1")




#for outside mTOR pathway only apply the 2 filters that were originally applied 
# == rare and not shared by proband and parent 
N030_tygir_mtor_rare1=dplyr::filter(N030_tygir, frequencyCategory == "rare" | frequencyCategory == "private")


#if there is the string "rna" in column "filename" - perform an action on column ND_N007P01_1_PB_Whole_FN007G1S2MN007P02P_KHTSC
#N030_tygir_mtor_rare1=subtract_match("filename", 14, "rna", N030_tygir_mtor_rare1)
write.csv(N030_tygir, "N013_mtor_nofil.csv")


write.csv(N030_tygir_mtor_rare1, "N036_mtor.csv")


#just need these like 40 lines 















setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents")
N030=read.csv("N030_varseqcalls.csv")
mtorgenes=read.csv("TSCdata/mTor_gene_list.csv")
N030$mtor=ifelse(N030$Gene.Names %in% mtorgenes$ï..Gene..153.genes., 1, 0)
#this WILL filer out the candidate TSC mutation since they should both share it - be forewarned

N030_unique=subset(N030, (as.character(X0.1.Genotypes..GT..Proband) != as.character(X0.1.Genotypes..GT..Father)))
       
N030_unique=filter(N030_unique, pLI > 0.9)
N030_unique=filter(N030_unique, mtor=="0")
write.csv(N030_unique, "N030_varseqcalls_outsidemtor_filtered.csv")

##################################################################################################START HERE
#compare with TyGIR calls - which is what I think this file is that I found in /LNTS/DOD/Analysis/Family_Excels 
#SOME of these are from /Family_Excels - some I had to rerun in tygir bc they were already filtered by mtor 
##tygir datasets have more info than varseq FATTHAM and CADD scores 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata/tygir_calls/")
##all the variables are N030 bc thats what I practiced with but need to cycle throug each one 
####CHANGE VAR NAME HERE 

#atch for 3 family sets with siblings - only F2 and F3
#30,36,29 are all diff format use next section for them
#N030_tygir=read.csv("ND_N035_HaplotypeCaller.csv")
N030_tygir=read.csv("F10.csv")

mtorgenes=read.csv("../mTor_gene_list.csv")


#toggle == and != to look for within mtor pathway vs outside mtor pathway
#not all have this use ifelse apply filter
#N030_tygir_mtor=filter(N030_tygir, mTor.pathway.gene. == "#N/A")

N030_tygir$mtor=ifelse(N030_tygir$Gene..HUGO. %in% mtorgenes$ï..Gene..153.genes., 1, 0)
N030_tygir=filter(N030_tygir, N030_tygir$mtor=="0")


#for outside mTOR pathway only apply the 2 filters that were originally applied 
# == rare and not shared by proband and parent 
N030_tygir_mtor_rare1=filter(N030_tygir, frequencyCategory == "rare")






#TWO possible freq actegories - frequencyCategory and dbFreqCat
N030_tygir_mtor_rare1=filter(N030_tygir, frequencyCategory == "rare")
N030_tygir_mtor_rare1=filter(N030_tygir_mtor_rare1, pLI > 0.9)
#ok so maybe you couldve filtered it more by looking at just the DNA calls
#not the repeat RNA lines 
#ALSO - this spreadsheet design? HORRIBLE
#have one line for No30 family with 2 cols for DNA genotype per indiv
#and MULTIPLE indiivudal lines for specific family memebers for RNA genotype 
#N030_tygir_mtor_rare1=filter(N030_tygir_mtor_rare1, Type == "DNA")
#for ND's
#N030_tygir_mtor_rare1=subset(N030_tygir_mtor_rare1, (as.character(N030_tygir_mtor_rare1[, 14]) != as.character(N030_tygir_mtor_rare1[, 15])))
#for Fs
N030_tygir_mtor_rare1=subset(N030_tygir_mtor_rare1, (as.character(N030_tygir_mtor_rare1[, 13]) != as.character(N030_tygir_mtor_rare1[, 16])))

#cant use this as col names are different each one 
#N030_tygir_mtor_rare1=subset(N030_tygir_mtor_rare1, (as.character(ND_N030P01) != as.character(ND_N030P03)))
#https://docs.varsome.com/meaning-of-the-mutationtaster-score 
#A means Autoamtically disease causing 
N030_tygir_mtor_rare1=dplyr::filter(N030_tygir_mtor_rare1, grepl('D|A',MutationTaster..A.Automatic..D.Deleterious.P.Polymorphism.N.Neutral.))
N030_tygir_mtor_rare1=dplyr::filter(N030_tygir_mtor_rare1, grepl('D',FATHMM..D..Damaging..T.Tolerated.))

####CHANGE VAR NAME HERE 
write.csv(N030_tygir_mtor_rare1, "F10_tygir_mtor_rare1.csv")



#this one does not help
N030_tygir_mtor_rare2=filter(N030_tygir_mtor, dbFreqCat == "rare") 






#For 30,36,29 - these are all in slightly different format - but 29 and 36 are not available on tygir 
#36 also requries some tweaking - cant match sibling either

##tygir datasets have more info than varseq FATTHAM and CADD scores 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata/tygir_calls/")
##all the variables are N030 bc thats what I practiced with but need to cycle throug each one 
####CHANGE VAR NAME HERE 
N030_tygir=read.csv("ND_N036_HaplotypeCaller.csv")
mtorgenes=read.csv("../mTor_gene_list.csv")


#toggle == and != to look for within mtor pathway vs outside mtor pathway
#not all have this use ifelse apply filter
#N030_tygir_mtor=filter(N030_tygir, mTor.pathway.gene. == "#N/A")

N030_tygir$mtor=ifelse(N030_tygir$Gene..HUGO. %in% mtorgenes$ï..Gene..153.genes., 1, 0)
N030_tygir=filter(N030_tygir, mtor=="0")

#TWO possible freq actegories - frequencyCategory and dbFreqCat
N030_tygir_mtor_rare1=filter(N030_tygir, frequencyCategory == "rare")
N030_tygir_mtor_rare1=filter(N030_tygir_mtor_rare1, pLI > 0.9)
#ok so maybe you couldve filtered it more by looking at just the DNA calls
#not the repeat RNA lines 
#ALSO - this spreadsheet design? HORRIBLE
#have one line for No30 family with 2 cols for DNA genotype per indiv
#and MULTIPLE indiivudal lines for specific family memebers for RNA genotype 
N030_tygir_mtor_rare1=filter(N030_tygir_mtor_rare1, Type == "DNA")

N030_tygir_mtor_rare1=subset(N030_tygir_mtor_rare1, (as.character(N030_tygir_mtor_rare1[, 16]) != as.character(N030_tygir_mtor_rare1[, 17])))
#cant use this as col names are different each one 
#N030_tygir_mtor_rare1=subset(N030_tygir_mtor_rare1, (as.character(ND_N030P01) != as.character(ND_N030P03)))
#https://docs.varsome.com/meaning-of-the-mutationtaster-score 
#A means Autoamtically disease causing 
N030_tygir_mtor_rare1=dplyr::filter(N030_tygir_mtor_rare1, grepl('D|A',MutationTaster..A.Automatic..D.Deleterious.P.Polymorphism.N.Neutral.))
N030_tygir_mtor_rare1=dplyr::filter(N030_tygir_mtor_rare1, grepl('D',FATHMM..D..Damaging..T.Tolerated.))

####CHANGE VAR NAME HERE 
write.csv(N030_tygir_mtor_rare1, "N036_tygir_mtor_rare1.csv")


#manually gluing all together and looking for repeats
all=read.csv("groupedall.csv")
all_sub=subset(all,duplicated(Gene..HUGO.))


all[duplicated(all[,all$Gene..HUGO.]),]




