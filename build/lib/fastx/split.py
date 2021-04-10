## modules
import os
import gzip
from Bio import SeqIO
from Bio.Seq import Seq

def chunkstring(string, string_length):
    chunks = [string[i:i+string_length] for i in range(0,len(string),string_length)]
    return(chunks)


def fasta_split(infile, outdir):
    fasta = (
        SeqIO.parse(infile, "fasta")
        if infile.endswith((".fa", ".fasta"))
        else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    )
    for seq in fasta:
        outfile = open(os.path.join(outdir, "{}.fa".format(seq.id.replace("/", "_"))), "w")
        outfile.write(">{}\n".format(seq.id))
        for chunk in chunkstring(seq.seq, 50):
            outfile.write("{}\n".format(chunk))
        outfile.close()

def fastq_split(infile, outdir):
    seqfile = (
        open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    )
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  ## header
            seq_id = j.strip()
            outfile = open(os.path.join(outdir, "{}.fq".format(seq_id.replace("@", "").replace("/", "_"))), "w")
        elif k == 1:  ## sequence
            seq = j.strip()
            seq_len = len(seq)
        elif k == 2:
            continue  ## plus
        elif k == 3:  ## quality
            seq_bq = j.strip()
            outfile.write("{}\n{}\n+\n{}\n".format(seq_id, seq, seq_bq))    
            outfile.close()

def seq_split(infile, outdir):
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir): os.mkdir(outdir)
    if infile.endswith((".fa", ".fa.gz", ".fasta")):
        fasta_split(infile, outdir)
    elif infile.endswith((".fq", ".fq.gz", ".fastq")):
        fastq_split(infile, outdir)
