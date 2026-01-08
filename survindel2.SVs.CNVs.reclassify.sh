#!/bin/bash 


refSeq="/home/andreas/Bos_taurus.30bulls.illuminaPE/data/refseq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fasta"
bamPath="/home/andreas/Bos_taurus.30bulls.illuminaPE/results/MQ.tag.bam"


cd /home/andreas/Bos_taurus.30bulls.illuminaPE/results/SVs.CNVs.survindel2.vcf
for bamFile in `realpath /home/andreas/Bos_taurus.30bulls.illuminaPE/results/MQ.tag.bam/*.bam`; do 
	baseName=`basename ${bamFile} | awk -F"." '{print $1}'`

	mkdir -p $baseName 
	cd $baseName

	# Run SV detection 
	# python survindel2.py --threads N_THREADS BAM_FILE WORKDIR REFERENCE_FASTA
	singularity exec /home/andreas/container/survindel2_8753dd8036f6d9df.sif python /opt/conda/bin/survindel2.py --threads 32 $bamFile . $refSeq  
	gunzip -c out.vcf.gz > out.vcf

	# reclassify output vcf with ML model
	# python run_classifier.py WORKDIR/out.vcf WORKDIR/out.pass-ml.vcf.gz WORKDIR/stats.txt ALL ml-model/
	singularity exec /home/andreas/container/survindel2_8753dd8036f6d9df.sif \
	python /opt/conda/bin/survindel2_run_classifier.py \
	out.vcf ${baseName}.OUT.PASS.ML.model.vcf.gz stats.txt ALL /home/andreas/git/SurVIndel2/ml-model 

	cd /home/andreas/Bos_taurus.30bulls.illuminaPE/results/SVs.CNVs.survindel2.vcf
done


