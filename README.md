# Understanding the gene expression differences between severe and mild forms of Tuberous Sclerosis Complex (TSC) and the role of genetic modifiers to identify novel drug targets 

### Overview

Tuberous sclerosis complex (TSC) is an autosomal dominant neurogenetic disorder, in which affected individuals are heterozygous for a mutation in either the TSC1 or TSC2 gene, thereby causing constitutive activation of the mTOR signaling pathway. Clinically, TSC patients present with a variety of symptoms including skin lesions, renal angiomyolipomas (AML), seizures, and cognitive delay. However, the severity of symptoms in TSC is variable, although mental retardation and intractable epilepsy are common. Phenotypic variability (PV) in TSC is well recognized, but its molecular basis is not understood. Variability in phenotype can be seen even within a single family (intrafamilial phenotypic variability [IPV]), where all affected individuals have the same TSC gene mutation. IPV clearly indicates that TSC1 or TSC2 mutation alone does not account for the observed phenotype and suggests a potential role for modifying factors. Mechanisms currently proposed for PV in monogenic disorders include modifying effects of unlinked genes (genetic modifiers), epigenetic factors, mosaicism, allelic skewing of gene expression, and environmental or other stochastic factors.  

Previously, we described the first report of comprehensive exome and transcriptome sequencing analysis of four TSC families demonstrating phenotypic variability. Here, we evaluated additional familial TSC cases comprised of severely affected children and mildly affected parents and/or siblings, performing whole exome and RNA sequencing of blood leukocytes. Utilizing all combined familial datasets, we aimed to (1) identify aberrant gene expression, aberrant splicing, and mono-allelic expression in severe vs mildly afflicted patients and (2) identify primary TSC mutations and potential modifying aberrations in the mTOR pathway that may drive disease variability. 

Final presentation slides with figures can be found [here](https://github.com/aahowel3/TSC_data/blob/main/Tgen_heliosscholars_finalpresentation_forgithub.pdf).

# Data Processing 
`LC_scratch.sh` includes the important regtools junctions extract step which creates a .junc file for each bam and mashes it altogether into test_juncfiles.txt
`LC_script_PB_LNTS.sh` includes the leafcutter_cluster_regtools.py step which creates several testSEVvsMILD files
really can just combine the two of these into 1

To visualize results with leafviz: 
perl /home/ahowell/leafcutter/leafviz/gtf2leafcutter.pl -o /home/ahowell/leafcutter/leafviz/annotation_codes/gencode_hg19/gencode_hg19 Homo_sapiens.GRCh37.87.gtf.gz

Rscript /home/ahowell/leafcutter/leafviz/prepare_results.R -m groups_file_nopaths_batches.txt ../testSEVvsMILD_perind_numers.counts.gz leafcutter_ds_cluster_significance.txt leafcutter_ds_effect_sizes.txt /home/ahowell/leafcutter/leafviz/annotation_codes/gencode_hg19/gencode_hg19

Final visualization (must be run in leafviz folder) 
ahowell@aahowel3-OptiPlex-7440-AIO:~/leafcutter/leafviz$ ./run_leafviz.R /home/ahowell/Documents/leafcutter_analysis/with_batches_and_exon/leafviz.RData

# March 2022 data rerun
Using the library EnsDb.Hsapiens.v75 vs the org.Hs.eg.db for AnnotationDBI does NOT change the output of the ORA, GGEA, SPIA - only cleans up NAs in the DGE list

# ASGAL for N039 
Creating a fasta of only the target chromosome - https://bioboot.github.io/web-2016/class-material/day3-fasta-practice.html 
The fasta name comes BEFORE the chromosme name argument! Also head the .fai file to see how to call the chromosome names - not necssarily the long thing you see when you grep "^>" 

Used HOMO_SAPIENS from desktop Documents/asgal_test .gtf NOT the long ERCV.III ONE 

Using all the DNA and RNA fastqs combined - from Daniel "The FASTQ files with "_TSRGD_" in their name are RNA and those with "_KHTSC_" are DNA (exome)."

# July 2022 data rerun
Family classifications have changed - N007 both parent and child are severe, N036 both children are now severe - details outlined in TSC_DOD_Meeeting_2022_F powerpoint, N029 seems to still have 1 severe child, 1 mild parent

Families N025, N027, N030, N036, N039 needed to be rerun through the alignment pipeline as they were mislabeled in their strandedness - They are all located here /labs/C4RCD/NarayananLAB/SampathRangasamy/TSC_Fibroblast/results/. They will have "KHRRE_RERUN" within their name. (not all have both an SK or PB rerun - just check which ones have rerun)

Primary scripts are:
`TSC_deseq2_analysis_withpreDod` - original differential expression - with batch effect control, functional annotation 
`filteringcalls_dbsnp_tester_withsnpeff` - functional version on hines
`edgeR_analysis.R` - not actually in TSCdata folder - in main Documents folder 

Update with Sampath 7/13/21:
We are including N029 – so 13 families plus N029 – there is a deletion in one copy of the TSC1 gene – from the 5’UTR to exon 21 
