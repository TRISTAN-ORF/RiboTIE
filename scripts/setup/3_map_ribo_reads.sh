#!/bin/bash

## prefetch --option-file sralist.txt
while read -r line
do
	# init
	dataset=$(echo $line | cut -f1)
	adapter=$(echo $line | cut -f2)	
	echo $dataset
	mkdir $dataset
	cd $dataset
	mkdir out 
	mkdir out/temp
	trimmed="out/temp/${dataset}_trimmed.fq"
	cleaned="out/${dataset}_trimmed_noXRNAs.fastq"

	# trim files and perform fastqc
	cutadapt -j 20 -m 20 -a $adapter ${dataset}.fastq -o "out/temp/${dataset}_trimmed.fq" > "out/temp/${dataset}_trimmed_report.txt"

	# remove rRNA/tRNA/smRNA/smoRNA
	STAR --genomeLoad NoSharedMemory --seedSearchStartLmaxOverLread .5   --genomeDir /data/jimc/genomes/GRCh38/v107/STAR/excl_RNA --readFilesIn $trimmed --outFilterMultimapNmax 1000 --outFilterMismatchNmax 2 --outFileNamePrefix out/temp/ --runThreadN 20 --outReadsUnmapped Fastx
	mv out/temp/Unmapped.out.mate1 $cleaned 
	fastqc -t 10 $cleaned
	
	# align to genome, output mapping to transcriptome as well
	STAR   --runThreadN 20   --genomeDir /data/jimc/genomes/GRCh38/v107/STAR   --genomeLoad NoSharedMemory   --readFilesIn $cleaned --outFileNamePrefix out/ --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outSAMattributes MD NH --outFilterMultimapNmax 10  --outMultimapperOrder Random --outFilterMismatchNmax 2   --seedSearchStartLmaxOverLread 0.5   --alignEndsType EndToEnd --outWigType bedGraph

	# cleanup
	mv out/Aligned.sortedByCoord.out.bam out/${dataset}_aligned.bam
	mv out/Aligned.toTranscriptome.out.bam out/${dataset}_aligned_tran.bam
	fastqc -t 10 out/${dataset}_aligned.bam
	samtools view -h -o out/${dataset}_aligned.sam out/${dataset}_aligned.bam
	samtools view -h -o out/${dataset}_aligned_tran.sam out/${dataset}_aligned_tran.bam
	cd ..
done < "metadata.txt"

