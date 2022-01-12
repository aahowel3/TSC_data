# TSC_data

From tutorial: https://davidaknowles.github.io/leafcutter/articles/Usage.html 
ON ISILON
LC_scratch.sh includes the important regtools junctions extract step which creates a .junc file for each bam and mashes it altogether into test_juncfiles.txt
LC_script_PB_LNTS.sh includes the leafcutter_cluster_regtools.py step which creates several testSEVvsMILD files
really can just combine the two of these into 1

ON LOCAL DESKTOP (not laptop) in ahowell@aahowel3-OptiPlex-7440-AIO:~/Documents/leafcutter_analysis/with_batches_and_exon
(had to run leafcutter R steps on seperate computer that would allow R installation - moved _perind_numers.counts.gz and groups.txt file) 
Final step for leafcutter: 
Rscript /home/ahowell/Documents/leafcutter/scripts/leafcutter_ds.R --num_threads 4 -e Homo_sapiens.GRCh37.87.exons.txt.gz ../testSEVvsMILD_perind_numers.counts.gz groups_file_nopaths_batches.txt

To visualize results with leafviz: 
Preperation steps: 
perl /home/ahowell/leafcutter/leafviz/gtf2leafcutter.pl -o leafcutter-annot Homo_sapiens.GRCh37.87.gtf.gz
Rscript /home/ahowell/leafcutter/leafviz/prepare_results.R -m groups_file_nopaths_batches.txt ../testSEVvsMILD_perind_numers.counts.gz leafcutter_ds_cluster_significance.txt leafcutter_ds_effect_sizes.txt leafcutter-annot

Final visualization (must be run in leafviz folder) 
ahowell@aahowel3-OptiPlex-7440-AIO:~/leafcutter/leafviz$ ./run_leafviz.R /home/ahowell/Documents/leafcutter_analysis/with_batches_and_exon/leafviz.RData
