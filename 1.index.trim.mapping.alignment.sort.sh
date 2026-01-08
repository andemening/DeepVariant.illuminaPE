#!/bin/bash -l

################################################################################################################################
###		Bos taurus 30 samples PE illumina 2024 
###		Andreas E Lundberg; andreas.e.lundberg@slu.se; anlu4433@student.uu.se 
### 	
### 	Index. Mapping, alignment. Samtools sort. 
###		fastp, bwa-mem2, samtools
###		
###		ARS-UCD1.2_bosTau9
###		
###		
###		
################################################################################################################################


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
# symlink tmp directory to fast storage SSD/NVMe disk if availible 


####### Variables & paths
# absolute path to base directory of project
projectPath="/home/andreas/Bos_taurus.30bulls.illuminaPE"
binPath="/home/andreas/bin"
containerPath="/home/andreas/container"

		# use 1.5x threads of node, 32 threads => 48 threads
		threadsNode="48" 	# 1.5x 32 threads
		fastpThreads="6"
		sortThreads="10"
		bwaThreads="32" 		# "((${threadsNode}-${fastpThreads}-${sortThreads}))"


# fasta refeq file 
refSeq="${projectPath}/data/bwa-meme.refseq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fasta"


# indexPrefix="Bos_taurus.ARS-UCD1.2"



fastqPath="${projectPath}/data/fastq"  				# how to handle all fastq files when they're not yet accessible inside singularity container ? 
													# symlink in preparation before source files are "present" in $fastqPath maybe, from /mnt/fastq 
													# --bind src:dest = 
													# $SINGULARITY_BIND variable 
													# export SINGULARITY_BIND="/opt,/data:/mnt" 
													# $binPath/fastp -i $fastqPath/${sampleName}_R1_001.fastq.gz -I $fastqPath/${sampleName}_R2_001.fastq.gz
													# SH-2272-SE-423*.fastq.gz 
													# for sampleName in `
													# `cat metadata.samples.txt | cut -d',' -f 4`
													# for sampleName in `ls ${fastqPath}/*.fastq.gz | awk -F"/" '{print $NF}' | awk -F"_" '{print $1 "_" $2 "_" $3}' | uniq`; do 




bamPath="${projectPath}/results/bam" 
vcfPath="${projectPath}/results/vcf"
statsPath="${projectPath}/results/stats"
tmpPath="${projectPath}/tmp" # use faster NVMe disk 

mkdir -p $bamPath $vcfPath $statsPath $tmpPath

## IMPORTANT
# Bind Singularity directories for access inside container, only needed for fastp then piped output ...
# ... otherwise the fastq files are accessible for Singularity program  
export SINGULARITY_BIND="/home/martin/fastq/scilifelab_30_bulls_2019/ALL,/home/martin/fastq/scilifelab_30_bulls_2019/ALL_SH2271,/home/martin/fastq/scilifelab_30_bulls_2019/ALL_SH2272"

################################################################################################################################
# Declaration of (empty) Read Group variables 
readgroupHeader=""		
idRG=""					# from fastq header: Instrument + run number + Flowcell id + lane number = A00181.113.HLGKLDSXX.3 				for fastq SH-2272-SE-423 = FJLSWEF423
smRG=""					# Sample name, id1 + id2 from metadata file = FJLSWEF423.SKBSWEF020225000186 																
lbRG=""					# library identifier: idRG + sample name = A00181.113.HLGKLDSXX.3.BOHSWEF458
puRG=""					# platform unit: flowcell + lane + (dual index, i.e. + sign) barcode = HLGKLDSXX.CCGTGAAG+CAGTGGAT.3 			{FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}

################################################################################################################################

## Indexing 
# samtools faidx $refSeq

# bwa-meme index -a meme -t 32 $refSeq
# Run other command (build_rmis_dna.sh) inside bwa-meme container, instead of the standard linked container named command (bwa-meme) 
# singularity exec /home/andreas/container/bwa-meme_1.0.6--06fa6e8cadb53114.sif build_rmis_dna.sh $refSeq 

# gatk CreateSequenceDictionary R=$refSeq


################################################################################################################################
## Trim/QC, Map/align, Sort pipe 

##	For loop commands per Sample
for sampleName in `ls ${fastqPath}/*.fastq.gz | awk -F"/" '{print $NF}' | awk -F"_" '{print $1 "_" $2 "_" $3}' | uniq`; do
	# sampleName variable example is: SH-2270-SE-402_S8_L003


	# debug log output of variables
	# echo -e "${sampleName}\t$bamPath/${sampleName}" 2>&1 | tee -a $statsPath/${sampleName}.debug.log
	# echo -e "\t$fastqPath/${sampleName}_R1_001.fastq\t$bamPath/${sampleName}.sam" 2>&1 | tee -a $statsPath/${sampleName}.debug.log

	
	# use current sample name to match and find unique id1 in metadata
	shortName=`echo ${sampleName} | cut -f 1 -d"_" `
	idRG=`zcat ${fastqPath}/${sampleName}_R1_001.fastq.gz | head -n1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/./g'` 							# id from fastq header, e.g. A00181_113_HLGKLDSXX_3
	smRG=`awk -F"," '$4==sampleAwk {print $1 "." $2}' sampleAwk="$shortName" ${projectPath}/data/metadata/metadata.samples.txt` 				# sample names from metadata, e.g. FJLSWEF423.SKBSWEF020225000186
	lbRG=`echo "${idRG}.${smRG}"` 																												# DNA library prep, e.g. A00181_113_HLGKLDSXX_3_BOHSWEF458
	puRG=`zcat ${fastqPath}/${sampleName}_R1_001.fastq.gz | head -n1 | awk -v OFS="." -F":" '{print $3, $4, $10}'  | sed 's/@//'`				# flowcell_barcode + lane + sample_barcode = HN7GGDSXX.4.AATCCGGA+CTACAGTT
	readgroupHeader="@RG\tID:${idRG}\tLB:${lbRG}\tPL:ILLUMINA\tSM:${smRG}\tPU:${puRG}"
	
	# -R "@RG\tID:$idRG\tLB:$lbRG\tPL:ILLUMINA\tSM:$smRG\tPU:$puRG" \
	#	-R "@RG'\t'ID:${idRG}'\t'LB:${lbRG}'\t'PL:ILLUMINA'\t'SM:${smRG}'\t'PU:${puRG}"

	
	singularity exec $containerPath/fastp_0.23.4--4ea6310369653ec7.sif fastp -i $fastqPath/${sampleName}_R1_001.fastq.gz -I $fastqPath/${sampleName}_R2_001.fastq.gz \
	    --stdout --thread $fastpThreads \
	    -j "${statsPath}/fastp.${sampleName}.json" \
	    -h "${statsPath}/fastp.${sampleName}.html" \
	    2> "${statsPath}/fastp.${sampleName}.log" \
	| $binPath/bwa-meme mem \
		-7 \
		-v 2 \
		-M \
		-p \
		-Y \
		-t $bwaThreads \
	    -R $readgroupHeader \
	    $refSeq - 2> $statsPath/bwa-meme.${sampleName}.log \
	| ${binPath}/samtools sort -T $tmpPath \
	    -m 10G \
	    -@ $sortThreads \
	    -O BAM \
	    -o ${bamPath}/${sampleName}.bam \
	    2> ${statsPath}/samtools.${sampleName}.log

done
		# samtools -m per thread memory allocation, set to 75% of absolute max RAM limit (~300GB ?) 

# Runtime memory is peaking (for first sample SH-2270-SE-402) at 235 237  GB  

 

################################################################################################################################
# REMOVE BEFORE RUN SCRIPT 

## Mapping pipe 
# 	bwa-mem2 mem -7 -t $threadsNode $refSeq $fastqPath/${sampleName}_R1_001.fastq.gz $fastqPath/${sampleName}_R2_001.fastq.gz > $bamPath/${sampleName}.sam

# -7 = use option to deploy bwa-meme (i.e. use Learned index for seeding) 
# -v = verbose level 2 
# -p = Smart pairing 
# -M = Mark shorter split hits as secondary
# -Y = soft clipping for CIGAR string 


