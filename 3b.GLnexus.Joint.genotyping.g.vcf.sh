#!/bin/bash -l

################################################################################################################################
###		Bos taurus 30 samples PE illumina 2024 
###		Andreas E Lundberg; andreas.e.lundberg@slu.se; anlu4433@student.uu.se 
### 	
### 	GLnexus 
###		Joint genotyping of .g.vcf files
###		ARS-UCD1.2_bosTau9
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

# fasta refeq file 
# refSeq="$projectPath/data/refseq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fasta"
	# refSeq="$projectPath/data/bwa-meme.refseq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fasta"


vcfPath="$projectPath/results/deepvariant.vcf"
statsPath="$projectPath/results/stats"
tmpPath="$projectPath/tmp" # use faster NVMe disk 
populationVcfPath="$projectPath/results/glnexus.joint.genotyping.g.vcf.stats"

# mkdir -p ${populationVcfPath}

################################################################################################################################
# samplesGVCF=`realpath $vcfPath 
realpath ${vcfPath}/*.g.vcf.gz > ${populationVcfPath}/samples.files.DeepVariant.g.vcf.txt
${binPath}/glnexus_cli \
		--config DeepVariantWGS \
		--list ${populationVcfPath}/samples.files.DeepVariant.g.vcf.txt \
		--mem-gbytes 350 \
	| ${binPath}/bcftools view - \
	| bgzip -@ 4 -c > ${populationVcfPath}\30samples.Bos_taurus.GLnexus.joint.genotyping.vcf.gz 
		# -d $tmpPath \


${binPath}/bcftools index \
	--threads 32 \
	--csi \
	--output ${populationVcfPath}\30samples.Bos_taurus.GLnexus.joint.genotyping.vcf.gz


