#!/bin/bash -e
#first create a folder in experiments with experiment name, in that folder create a folder called "raw"
#paste the raw data to be analyzed in the raw folder

#experimentName="CV001N"
#experimentName="CV001L3"
experimentName="CV001L5S3"

for e in $(ls "../data/experiments/${experimentName}/raw/" | sed "s/\.fastq//")
do
	# e here needs to be of the form 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample1_B_uniformis_atcc_8492
	# Specifically, the part 'B_uniformis_atcc_8492' is important as I'm getting usedStrain based on the number of underscores.
	usedStrain=$(echo $e | rev | cut -d "_" -f1-4 | rev |sed "s/\.fastq//")
	sampleID=$(echo $e | sed "s/_${usedStrain}//")
	
	# as alternative if file name is different from above, add the usedStrain manually in the form of B_uniformis_atcc_8492

	echo $usedStrain
	echo $sampleID

	# Transform genome into long format
	cat ../data/genomes/${usedStrain}.fna | grep -v "^>" |  tr -d "\n" | sed "s/./\0\n/g" > ../data/genomes_long_format/${usedStrain}_long.txt

	# run prokka
	#source activate prokka
	#mkdir -p ../results/prokka; prokka --outdir ../data/prokka/ --force --prefix ${usedStrain}  --cpus 8 ../data/genomes/${usedStrain}.fna
	#cat ../data/prokka/${usedStrain}.gff | grep "ID=" > tmpfile; mv tmpfile ../data/prokka/${usedStrain}.gff
	#conda deactivate

	# Filter reads and get barcodes. TODO: Remove low quality reads?!
	mkdir -p ../data/experiments/${experimentName}/filtered
	mkdir -p ../results/out/${experimentName}/
	mkdir -p ../results/bam/${experimentName}/
	python filter_reads_and_get_barcodes.py ../data/experiments/${experimentName}/raw/${sampleID}_${usedStrain}.fastq ../data/experiments/${experimentName}/filtered/${sampleID}_${usedStrain}_filtered.fastq ../results/out/${experimentName}/readname_barcode_${sampleID}_${usedStrain}.tab

	# map reads
	#source activate bowtie2
	bowtie2 -p 4 -x ../data/genomes/${usedStrain} -U ../data/experiments/${experimentName}/filtered/${sampleID}_${usedStrain}_filtered.fastq -S ../results/bam/${experimentName}/${sampleID}_${usedStrain}.bam
	#conda deactivate # Could also be conda deactivate?

	# Sort bam
	#source activate samtools
	samtools sort ../results/bam/${experimentName}/${sampleID}_${usedStrain}.bam > ../results/bam/${experimentName}/${sampleID}_${usedStrain}_sorted.bam

	# Clean up
	rm ../results/bam/${experimentName}/${sampleID}_${usedStrain}.bam

	# Get depth along each position. Cap the depth at 1000 for simplicity. Save as "depth" file
	samtools depth -a -d 1000 ../results/bam/${experimentName}/${sampleID}_${usedStrain}_sorted.bam > ../results/out/${experimentName}/depth_${sampleID}_${usedStrain}.tab

	# Get important columns and save as "insertionInfoRaw" file
	samtools view ../results/bam/${experimentName}/${sampleID}_${usedStrain}_sorted.bam | cut -f 1,2,4,5,6,10 > ../results/out/${experimentName}/insertionInfoRaw_${sampleID}_${usedStrain}.tab
	#conda deactivate

	# Get read length directly (parsing in R is quite painfully slow...). This creates a file with the important columns of the bam file and read length info
	paste <(cat ../results/out/${experimentName}/insertionInfoRaw_${sampleID}_${usedStrain}.tab) <(cat ../results/out/${experimentName}/insertionInfoRaw_${sampleID}_${usedStrain}.tab | cut -f5 | sed "s/M$//") > ${sampleID}_${usedStrain}.tmp; mv ${sampleID}_${usedStrain}.tmp ../results/out/${experimentName}/insertionInfoRaw_${sampleID}_${usedStrain}.tab
done

# Generate plots
# Rscript analyse_and_visualize.r
