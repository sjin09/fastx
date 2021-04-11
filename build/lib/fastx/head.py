## modules
import os
import gzip
import natsort
from Bio import SeqIO
from fastx.common import chunkstring

def fasta_head(infile, number, outfile):
    counter = 0
    fasta = (
        SeqIO.parse(infile, "fasta")
        if infile.endswith((".fa", ".fasta"))
        else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    )
    for seq in fasta:
        counter += 1
        outfile.write(">{}\n".format(seq.id))
        for chunk in chunkstring(seq.seq):
            outfile.write("{}\n".format(chunk))
        if counter == number:
            break
    

def fastq_head(infile, number, outfile):
    counter = 0
    seqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  ## header
            seq_id = j.strip()
        elif k == 1:  ## sequence
            seq = j.strip()
        elif k == 2: ## plus
            continue  
        elif k == 3:  ## quality
            seq_bq = j.strip()
            outfile.write("{}\n{}\n+\n{}\n".format(seq_id, seq, seq_bq))
            counter += 1
            if counter == number:
                break

def seq_head(infile, number, outfile):
    if infile.endswith((".fa", ".fa.gz", ".fasta")):
        fasta_head(infile, number, outfile)
    elif infile.endswith((".fq", ".fq.gz", ".fastq")):
        fastq_head(infile, number, outfile)