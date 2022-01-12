#!/usr/bin/env bash
#SBATCH --job-name="leafcutter"
#SBATCH --time=0-48:00:00
#SBATCH --mail-user=aahowel3@asu.edu
#SBATCH --mail-type=FAIL
#SBATCH -p isilon

module load samtools/1.9 

find /labs/C4RCD/NarayananLAB/SampathRangasamy/LNTS_DoD/Analysis/results/ -path /labs/C4RCD/NarayananLAB/SampathRangasamy/LNTS_DoD/Analysis/results/to_DELETE -prune -false -o -name *starAligned*.bam > log
for line in $(<log); do
        samtools index $line
	if [[ $line == *"KHRRE"* ]]; then 
	/home/ahowell/regtools/build/regtools junctions extract -s 2 -a 8 -m 50 -M 500000 $line -o ${line}.junc
	fi
	if [[ $line == *"TSRGD"* ]]; then
        /home/ahowell/regtools/build/regtools junctions extract -s 1 -a 8 -m 50 -M 500000 $line -o ${line}.junc
	fi 
        echo ${line}.junc >> test_juncfiles.txt
done


