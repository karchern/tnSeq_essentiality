from Bio import SeqIO
from Bio.Seq import Seq
from sys import argv
from regex import search

fastq = SeqIO.parse(str(argv[1]), "fastq")
fastqBaseName = str(argv[1]).replace(".fastq", "")
def filterRecords(rec):
    # Select reads with inverted repeat
    #if "TACGAAGACCGGGGACTTATCATCCAACCTGT" in str(rec.seq):
    # TODO: Allow for one mismatch
    reg =  search("(TACGAAGACCGGGGACTTATCATCCAACCTGT){s<=1}", str(rec.seq))
    if reg:
        if not "TAGAAAGGCGATGAAAGAGATAAAACGAAA" in str(rec.seq):
            # For mapping, I am splitting the reads at the transposon inverted repeat. This makes sure that each read maps exactly at it's insertion site.
            # It also means I don't have to do --trim5 in bowtie2 anymore.
            #for qualAtPos in rec.letter_annotations['phred_quality']:
            #    if qualAtPos < 30:
            #        return(None)
            #return(rec[(str(rec.seq).find("TACGAAGACCGGGGACTTATCATCCAACCTGT") + len("TACGAAGACCGGGGACTTATCATCCAACCTGT")):])
            return(rec[reg.span()[1]:])

#fastqOut = (record for record in fastq if filterRecords(record))
fastqOut = (filterRecords(record) for record in fastq)
fastqOut = (record for record in fastqOut if record is not None)
SeqIO.write(fastqOut, str(argv[2]), "fastq")

# Also get barcodes...
fastq = SeqIO.parse(str(argv[1]), "fastq")

def getBarcodes(rec):
    # Same filtering step as in filterRecords. Not very elegant,  but I couldn't think of another way to solve this without loading the whole fastq into memory.
    # In any case, The filtering steps are the same so we look for the Barcode in a subset of filtered_reads.
    reg =  search("(TACGAAGACCGGGGACTTATCATCCAACCTGT){s<=1}", str(rec.seq))
    if reg:
        if not "TAGAAAGGCGATGAAAGAGATAAAACGAAA" in str(rec.seq):
            bridgingSite = search("CAGAATTGGGAGTC", str(rec.seq))
            if bridgingSite:
                # TODO: Include phred score check!
                indices = ((bridgingSite.span()[0] - 25), bridgingSite.span()[0])
                barcode_phred_scores = rec.letter_annotations['phred_quality'][indices[0]:indices[1]]
                if all([True if phredscore >= 30 else False for phredscore in barcode_phred_scores]):
                    barCode = rec.seq[indices[0]:indices[1]]
                else:
                    barCode = "NA"
            else:
                barCode = "NA"
            return((str(rec.id), str(barCode)))

infoOut = (getBarcodes(record) for record in fastq)
infoOut = (record for record in infoOut if record is not None)
with open(str(argv[3]), "w") as f:
    for e in infoOut:
        f.write(e[0] + '\t' + e[1] + '\n')

