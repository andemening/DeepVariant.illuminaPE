#!/bin/bash 

# this script seems to take ~15min per sample from start to finish
# ... but all 30 samples took 15h, so something is not is not computing, maybe samtools index is a real time consumer (which is strange that I didn't know this)  

for bamFile in `realpath /home/andreas/Bos_taurus.30bulls.illuminaPE/results/bam/*.bam`; do 
	baseName=`basename ${bamFile} | awk -F".bam" '{print $1}'`
	cd /home/andreas/Bos_taurus.30bulls.illuminaPE/results/MQ.tag.bam
	sambamba sort $bamFile --sort-by-name --match-mates --memory-limit=340G --nthreads=48 --tmpdir=tmp --out=temp.${baseName}.query.sort.sambamba.bam
	
	sambamba view -h temp.${baseName}.query.sort.sambamba.bam --nthreads=8 | \
	samblaster --addMateTags --acceptDupMarks | \
	samtools view -b -S -u - | \
	sambamba sort /dev/stdin --memory-limit=3400G --nthreads=36 --tmpdir=tmp --out=${baseName}.samblaster.MQ.sambamba.sorted.bam

	# no need to index, sambamba sort does this on the fly 
	samtools index ${baseName}.samblaster.MQ.sambamba.sorted.bam

	rm temp.${baseName}.query.sort.sambamba.bam
done

