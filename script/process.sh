#!/bin/bash -e

#for sampleID in 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample1_sequence  000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample4_sequence 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample5_sequence 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample6_sequence 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample2_sequence 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample3_sequence CV001I10_all_forward_reads_concatenated_140000_reads_per_sample

experimentName="CV001N"

for sampleID in $(ls "../data/experiments/${experimentName}/" | sed "s/\.fastq//")
do

# Filter reads. TODO: Remove low quality reads?!
echo "Filtering reads and getting barcodes..."
python filter_reads_and_get_barcodes.py ../data/experiments/${experimentName}/${sampleID}.fastq ../data/experiments/${experimentName}/${sampleID}_filtered.fastq ../data/experiments/${experimentName}/${sampleID}_readname_barcode.tab

echo "mapping reads..."
source activate bowtie2
if [[ $sampleID == "000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample2_sequence" ]] || [[ $sampleID == "000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample3_sequence"  ]]
then
	echo "we need to map aginst the B. vulgatus genome in this case."
	bowtie2 -p 6 -x ../data/genomes/GCF_000154205.1_ASM15420v1_genomic -U ../data/reads/${experimentName}/${sampleID}_filtered.fastq -S ../results/bam/${sampleID}.bam
else
	bowtie2 -p 6 -x ../data/genomes/atcc_8492_concatenated -U ../data/reads/${experimentName}/${sampleID}_filtered.fastq -S ../results/bam/${sampleID}.bam
fi
source deactivate # Could also be conda deactivate?

# Sort bam
source activate samtools
samtools sort ../results/bam/${sampleID}.bam > ../results/bam/${sampleID}_sorted.bam

# Clean up
rm ${sampleID}.bam

# Get depth along each position. Cap the depth at 1000 for simplicity.
samtools depth -a -d 1000 ../results/bam/${sampleID}_sorted.bam > ../results/out/depth_${sampleID}.tab

# Get important columns
samtools view ../results/bam/${sampleID}_sorted.bam | cut -f 1,2,4,5,6,10 > ../results/out/insertionInfoRaw_${sampleID}.tab
source deactivate

# Get read length directly (parsing in R is quite painfully slow...)
paste <(cat ../results/out/insertionInfoRaw_${sampleID}.tab) <(cat ../results/out/insertionInfoRaw_${sampleID}.tab | cut -f5 | sed "s/M$//") > ${sampleID}.tmp; mv ${sampleID}.tmp ../results/out/insertionInfoRaw_${sampleID}.tab
done
