# modules
import gzip
from Bio import SeqIO


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
        if k == 0:  # header
            seq_id = j.strip()
            seq_id = seq_id if not isinstance(seq_id, bytes) else seq_id.decode("utf-8")
            seq_id = seq_id.replace("@", "")
            zmw = seq_id.split("/")[1]
        elif k == 1:  # sequence
            seq = j.strip()
            seq = seq if not isinstance(seq, bytes) else seq.decode("utf-8")
            seq_len = len(seq)
        elif k == 2:
            continue  # plus
        elif k == 3:  # quality
            outfile.write("{}\t{}\n".format(zmw, seq_len))


def seq_length(infile, outfile):
    if infile.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
        fasta_length(infile, outfile)
    elif infile.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
        fastq_length(infile, outfile)
