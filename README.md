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
