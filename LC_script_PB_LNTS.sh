#!/usr/bin/env bash
#SBATCH --job-name="leafcutter"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=aahowel3@asu.edu
#SBATCH --mail-type=FAIL

module load R/4.0.3

#find /labs/C4RCD/NarayananLAB/SampathRangasamy/LNTS_DoD/Analysis/results/ \
#-path /labs/C4RCD/NarayananLAB/SampathRangasamy/LNTS_DoD/Analysis/results/to_DELETE -prune -false -o -name *starAligned*.bam > log

#grep -v "N039" log > log1.txt 
#grep -v "SK" log1.txt > log2.txt
#double check this produces what you expect it to - left N015 in on accident manually edited it out of log2
#grep -v "N015P02" log2.txt > log2.1.txt 

#for line in $(<log2.txt); do 
#	samtools index $line
#	/home/ahowell/regtools/build/regtools junctions extract -a 8 -m 50 -M 500000 $line -o ${line}.junc
#	echo ${line}.junc >> test_juncfiles.txt
#done


#python2 ../leafcutter/clustering/leafcutter_cluster_regtools.py -j test_juncfiles.txt -m 50 -o testSEVvsMILD -l 500000 --checkchrom

Rscript ../leafcutter/scripts/leafcutter_ds.R --num_threads 4 testSEVvsMILD_perind_numers.counts.trunc.gz groups_file_nopaths.txt




#remove below this later
#for line in $(<log); do scp ahowell@dback-data1.tgen.org:$line ./dir/$line; done



#for bamfile in ahowell@dback-data1.tgen.org:/labs/C4RCD/NarayananLAB/SampathRangasamy/LNTS_DoD/Analysis/results/ND_N007_ps201810031700/TSRGD/*; do
    #echo $bamfile
    #echo Converting $bamfile to $bamfile.junc
    #samtools index $bamfile
    #regtools junctions extract -a 8 -m 50 -M 500000 $bamfile -o $bamfile.junc
    #echo $bamfile.junc >> test_juncfiles.txt
#done


