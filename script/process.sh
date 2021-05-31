#!/bin/bash -e

#for sampleID in 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample1_sequence  000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample4_sequence 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample5_sequence 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample6_sequence 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample2_sequence 000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample3_sequence
for sampleID in CV001I10_all_forward_reads_concatenated_140000_reads_per_sample
do
#sampleID="HMW5KAFX2_CV001L3_21s001573-1-1_Voogdt_lane1Sample2_sequence"

# Filter reads. TODO: Remove low quality reads?!
echo "Filtering reads and getting barcodes..."
python filter_reads_and_get_barcodes.py ${sampleID}.fastq ${sampleID}_filtered.fastq results/${sampleID}_readname_barcode.tab
#python filter_reads.py ${sampleID}.fastq ${sampleID}_filtered_stringent.fastq


# Transform reads to fasta
#paste <(cat ${sampleID}.fastq | paste - - - - | cut -f1) <(cat ${sampleID}.fastq | paste - - - - | cut -f2) | tr "\t" "\n" > data/${sampleID}.fasta

# Map reads to B. uniformis genome
# OLD
#bowtie2 -p 6 -x ~/for_carlos_3/seq/atcc_8492_concatenated --trim5 110 -U ${sampleID}_filtered.fastq -S ${sampleID}.bam
# New (with internal read splitting in filter_reads.py
source activate bowtie2
if [[ $sampleID == "000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample2_sequence" ]] || [[ $sampleID == "000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample3_sequence"  ]]
then
	echo "we need to map aaginst the B. vulgatus genome in this case."
	bowtie2 -p 6 -x ~/for_carlos_3/seq/GCF_000012825.1_ASM1282v1_genomic_atcc_8482 -U ${sampleID}_filtered.fastq -S ${sampleID}.bam
else
	bowtie2 -p 6 -x ~/for_carlos_3/seq/atcc_8492_concatenated -U ${sampleID}_filtered.fastq -S ${sampleID}.bam
fi
source deactivate # Could also be conda deactivate?

# Sort bam
source activate samtools
samtools sort ${sampleID}.bam > ${sampleID}_sorted.bam

# Clean up
rm ${sampleID}.bam

# Get depth along each position. Cap the depth at 1000 for simplicity.
samtools depth -a -d 1000 ${sampleID}_sorted.bam > results/depth_${sampleID}.tab

# Get the SAM flag, position, CIGAR flag and sequence of all reads.
# SAM flag is important because reads mapping on the antisense strand have to be 'reversed' in order to find the position of where the transposon actually inserted.
# Example: The following is the extracted columns for a reverse-mapping read (SAM flag 16):
# 16	69	51M	CTCCAGAAAGGAGGTGTTCCAGCCGCACCTTCCGGTACGGCTACCTTGTTA
# You acn see that the read ends in TA, so the actual insertion site is 69 + 51. This isn't necessary for forward reads
samtools view ${sampleID}_sorted.bam | cut -f 1,2,4,5,6,10 > results/insertionInfoRaw_${sampleID}.tab
source deactivate

# Get read length directly (parsing in R is quite painfully slow...)
paste <(cat results/insertionInfoRaw_${sampleID}.tab) <(cat results/insertionInfoRaw_${sampleID}.tab | cut -f5 | sed "s/M$//") > ${sampleID}.tmp; mv ${sampleID}.tmp results/insertionInfoRaw_${sampleID}.tab
done
