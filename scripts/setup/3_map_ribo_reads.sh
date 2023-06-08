#!/bin/bash
cd ../../data/ribo

while read -r line
do
	# init
	dataset=$(echo $line | cut -f 1 -d ' ')
	adapter=$(echo $line | cut -f 2 -d ' ')	
	echo $dataset
	mkdir $dataset
	cd $dataset
	mkdir out 
	mkdir out/temp
	mkdir out/genome
	trimmed="out/temp/${dataset}_trimmed.fq"
	cleaned="out/${dataset}_trimmed_noXRNAs.fastq"

	# download files (from ncbi repo)
	fasterq-dump $dataset
	fastqc -t 20 $dataset.fastq

	# trim files and perform fastqc
	cutadapt -j 20 -m 20 -a $adapter ${dataset}.fastq -o "out/temp/${dataset}_trimmed.fq" > "out/temp/${dataset}_trimmed_report.txt"
	fastqc -t 20 out/temp/${dataset}_trimmed.fq

	# remove rRNA/tRNA/smRNA/smoRNA
	STAR --genomeLoad NoSharedMemory --seedSearchStartLmaxOverLread .5   --genomeDir '../../genome/STAR/excl_RNA' --readFilesIn $trimmed --outFilterMultimapNmax 1000 --outFilterMismatchNmax 2 --outFileNamePrefix out/temp/ --runThreadN 20 --outReadsUnmapped Fastx
	mv out/temp/Unmapped.out.mate1 $cleaned 
	fastqc -t 20 $cleaned
	
	# align to genome, output mapping to transcriptome as well
	STAR   --runThreadN 20   --genomeDir '../../genome/STAR'   --genomeLoad NoSharedMemory   --readFilesIn $cleaned --outFileNamePrefix out/ --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outSAMattributes MD NH --outFilterMultimapNmax 10  --outMultimapperOrder Random --outFilterMismatchNmax 2   --seedSearchStartLmaxOverLread 0.5   --alignEndsType EndToEnd --outWigType bedGraph
		
	# cleanup
	mv out/Aligned.sortedByCoord.out.bam out/genome/${dataset}_aligned.bam
	mv out/Aligned.toTranscriptome.out.bam out/${dataset}_aligned_tran.bam
	fastqc -t 20 out/genome/${dataset}_aligned.bam
	fastqc -t 20 out/${dataset}_aligned_tran.bam
	samtools view -h -o out/genome/${dataset}_aligned.sam out/genome/${dataset}_aligned.bam
	samtools view -h -o out/${dataset}_aligned_tran.sam out/${dataset}_aligned_tran.bam
	cd ..
	
done < "metadata.txt"
