## modules
import gzip
import pysam
from Bio import SeqIO
from Bio.Seq import Seq


def fasta_length(infile, outfile):
    fasta = (
        SeqIO.parse(infile, "fasta")
        if infile.endswith((".fa", ".fasta"))
        else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    )
    for seq in fasta:
        outfile.write("{}\t{}\n".format(seq.id, len(seq.seq)))


def fastq_length(infile, outfile):
    seqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  ## header
            seq_id = j.strip().replace("@", "")
        elif k == 1:  ## sequence
            seq_len = len(j)
        elif k == 2:
            continue  ## plus
        elif k == 3:  ## quality
            outfile.write("{}\t{}\n".format(seq_id, seq_len))


def seq_length(infile, outfile):
    if infile.endswith((".fa", ".fa.gz", ".fasta")):
        fasta_length(fasta_length)
    elif infile.endswith((".fq", ".fq.gz", ".fastq")):
        fastq_length(infile, outfile)
