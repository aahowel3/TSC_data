# TSC_data

From tutorial: https://davidaknowles.github.io/leafcutter/articles/Usage.html 
ON ISILON
LC_scratch.sh includes the important regtools junctions extract step which creates a .junc file for each bam and mashes it altogether into test_juncfiles.txt
LC_script_PB_LNTS.sh includes the leafcutter_cluster_regtools.py step which creates several testSEVvsMILD files
really can just combine the two of these into 1

ON LOCAL DESKTOP (not laptop) in ahowell@aahowel3-OptiPlex-7440-AIO:~/Documents/leafcutter_analysis/with_batches_and_exon
(had to run leafcutter R steps on seperate computer that would allow R installation - moved _perind_numers.counts.gz and groups.txt file) 
Final step for leafcutter: 
to get the .exons.txt.gz file run gtf_to_exons.R documented in the leafcutter Differential Splicing tab 
Rscript /home/ahowell/Documents/leafcutter/scripts/leafcutter_ds.R --num_threads 4 -e Homo_sapiens.GRCh37.87.exons.txt.gz ../testSEVvsMILD_perind_numers.counts.gz groups_file_nopaths_batches.txt

To visualize results with leafviz: 
Preperation steps: 
SO - issue here - you just gave it a made annotation code when they presupply you with actual annotation codes - had to download them using ./download_human_annotation_codes.sh
Commands have been updated accordingly
perl /home/ahowell/leafcutter/leafviz/gtf2leafcutter.pl -o /home/ahowell/leafcutter/leafviz/annotation_codes/gencode_hg19/gencode_hg19 Homo_sapiens.GRCh37.87.gtf.gz

Rscript /home/ahowell/leafcutter/leafviz/prepare_results.R -m groups_file_nopaths_batches.txt ../testSEVvsMILD_perind_numers.counts.gz leafcutter_ds_cluster_significance.txt leafcutter_ds_effect_sizes.txt /home/ahowell/leafcutter/leafviz/annotation_codes/gencode_hg19/gencode_hg19

Final visualization (must be run in leafviz folder) 
ahowell@aahowel3-OptiPlex-7440-AIO:~/leafcutter/leafviz$ ./run_leafviz.R /home/ahowell/Documents/leafcutter_analysis/with_batches_and_exon/leafviz.RData

# March 2022 data rerun
Using the library EnsDb.Hsapiens.v75 vs the org.Hs.eg.db for AnnotationDBI does NOT change the output of the ORA, GGEA, SPIA - only cleans up NAs in the DGE list

# ASGAL for N039 
Creating a fasta of only the target chromosome - https://bioboot.github.io/web-2016/class-material/day3-fasta-practice.html 
The fasta name comes BEFORE the chromosme name argument! Also head the .fai file to see how to call the chromosome names - not necssarily the long thing you see when you grep "^>" 

USED HOMO_SAPIENS from desktop Documents/asgal_test .gtf NOT the long ERCV.III ONE 

Might have messed this up by using all the DNA and RNA fastqs combined - from Daniel "The FASTQ files with "_TSRGD_" in their name are RNA and those with "_KHTSC_" are DNA (exome)." but they're going to do a TSC1 and TSC2 panel on N039 anyway so it doesnt matter

# July 2022 data rerun
Family classifications have changed - N007 both parent and child are severe, N036 both children are now severe - details outlined in TSC_DOD_Meeeting_2022_F powerpoint, N029 seems to still have 1 severe child, 1 mild parent, will see if anything changes 

Families N025, N027, N030, N036, N039 needed to be rerun through the alignment pipeline as they were mislabeled in their strandedness - They are all located here /labs/C4RCD/NarayananLAB/SampathRangasamy/TSC_Fibroblast/results/. They will have "KHRRE_RERUN" within their name. (not all have both an SK or PB rerun - just check which ones have rerun)

Primary scripts are:
TSC_deseq2_analysis_withpreDod - original differential expression - with batch effect control, functional annotation 
filteringcalls_dbsnp_tester_withsnpeff - functional version on hines
edgeR_analysis.R - not actually in TSCdata folder - in main Documents folder 


Update with Sampath 7/13/21
-we are including N029 – so 13 families plus N029 – there is a deletion in one copy of the TSC1 gene – from the 5’UTR to exon 21 

# Important 
location of exome snpEFF.vcf list <br />
/labs/C4RCD/NarayananLAB/SampathRangasamy/LNTS_DoD/Analysis/dbnsfp_annotations/list_of_snpeffvcs_locs.txt <br />
location of script that runs those files through dbnsfp 4.3a <br />
/labs/C4RCD/NarayananLAB/SampathRangasamy/LNTS_DoD/Analysis/dbnsfp_annotations/dbnsfp_annotations.sh -> output files end in .dbnsfp4.3a_withnarg.vcf <br />
location of script that iterates dbsnfp files through the rscript to calculate their deleterious score <br />
(.dbnsfp4.3a_withnarg.vcf moved to hines) /work/aahowel3/TSC/july_2022/dnsfp_annotations/filteringcalls_dbsnp_tester_withsnpeff.sh <br />
Command line code is bash filteringcalls_dbsnp_tester_withsnpeff.sh - iterates through filteringcalls_dbsnp_tester_withsnpeff.R 

location of new htseq files that go into the differential expression
isilon - /home/ahowell/list_of_htseqs.txt

Remember – there are 3 different datasets that variant calling can be done on – exome, PB rna, and SK rna – I don’t think Daniel reran any of the exome variant calls – because the issue was with the rna strandedness data 
I recently noticed that the strand configuration for the RNA library prep kit Kapa RNA Hyper with RiboErase "KHRRE'' within KBase was incorrect. This issue would affect RNA quantification since we hardset the strandeness for the quantification tools like Salmon and HTSeq. This would affect 9 subjects (5 families -- N027 and N025 only have one subject so technically 3 families).
 
So you only need to go rooting around in the RERUN folders when you’re getting the htseq files – or if/when you run the rna variants through dbnsfp – 
Usually the exome vcfs are in a folder labelle “hc” whereaes the rna vcfs are in an “rnaHC” folder 

