#!/bin/bash -l

################################################################################################################################
###		Bos taurus 30 samples PE illumina 2024 
###		Andreas E Lundberg; andreas.e.lundberg@slu.se; anlu4433@student.uu.se 
### 	
### 	DeepVariant
###		SV/INDEL caller 
###		ARS-UCD1.2_bosTau9, dna.toplevel
###		
###		
###		
################################################################################################################################

# exit 1 
# REMOVE BEFORE RUN SCRIPT


####### Expected paths under project root 
#			.
#			├── bin
#			├── data
#			│   ├── fastq
#			│   └── refseq
#			├── results
#			│   ├── bam
#			│   ├── stats
#			│   └── vcf
#			└── tmp
#
# symlink tmp directory to fast storage SSD/NVMe disk


####### Variables & paths
# absolute path to base directory of project
projectPath="/home/andreas/Bos_taurus.30bulls.illuminaPE"
binPath="/home/andreas/bin"
containerPath="/home/andreas/container"

# make nicer and use 1.5x threads of node 
# threadsNode="60" 	# 1.5x 32 threads
# fastpThreads="6"
# sortThreads="14"
# bwaThreads="40" 		# "((${threadsNode}-${fastpThreads}-${sortThreads}))"

# fasta refeq file 
refSeq="$projectPath/data/refseq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fasta"
			# refSeq="$projectPath/data/bwa-meme.refseq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fasta"


# indexPrefix="Bos_taurus.ARS-UCD1.2"

# fastqPath="$projectPath/data/fastq" 
bamPath="$projectPath/results/bam" 
vcfPath="$projectPath/results/deepvariant.vcf"
logsPath="$projectPath/results/stats/deepvariant.vcf.logs"
tmpPath="$projectPath/tmp" # use faster NVMe/SSD disk if available, needs lots of storage

mkdir -p $vcfPath $logsPath

# DeepVariant variable
# export BIN_VERSION="1.6.1"
export TMPDIR="${projectPath}/tmp"

################################################################################################################################
## MarkDuplicates - Spark parallelization 


##	For loop commands per Sample
for sampleName in `ls ${bamPath}/*.QC.SORTED.BWA-MEME.MarkDuplicatesSpark.DEDUP.bam | awk -F"/" '{print $NF}' | awk -F"." '{print $1}'`; do
	# sampleName variable is shortened basename, i.e.: SH-2270-SE-402_S8_L003

	# echo -e "${sampleName}
	# debug log output of variables

	# singularity run -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"${BIN_VERSION}" /opt/deepvariant/bin/run_deepvariant \

	singularity run -B /usr/lib/locale/:/usr/lib/locale/ ${containerPath}/deepvariant_1.6.1.sif run_deepvariant \
		--model_type=WGS \
		--ref=${refSeq} \
		--reads=${bamPath}/${sampleName}.QC.SORTED.BWA-MEME.MarkDuplicatesSpark.DEDUP.bam \
		--output_vcf=${vcfPath}/${sampleName}.DeepVariant.vcf.gz \
		--output_gvcf=${vcfPath}/${sampleName}.DeepVariant.g.vcf.gz \
		--logging_dir "${logsPath}" \
		--runtime_report true \
		--num_shards=32

		# --intermediate_results_dir "/home/andreas/Bos_taurus.30bulls.illuminaPE/deepvariant.vcf/intermediate_results_dir" \

done


