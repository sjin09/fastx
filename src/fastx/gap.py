## modules
import os
import gzip
from Bio import SeqIO


def seq_gaps(infile, outfile):
    fasta = (
        SeqIO.parse(infile, "fasta")
        if infile.endswith((".fa", ".fasta"))
        else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    )
    status = 0
    for seq in fasta:
        seq_len = len(seq.seq) - 1
        for i, j in enumerate(seq.seq):
            if status == 0:
                if j == "N":
                    start = i
                    status = 1
                else:
                    status = 0
            elif status == 1:
                if j == "N":
                    if i != seq_len:
                        status = 1
                    elif i == seq_len:
                        status = 0
                        end = i + 1
                        outfile.write("{}\t{}\t{}\n".format(seq.id, start, end))
                else:
                    end = i
                    status = 0
                    outfile.write("{}\t{}\t{}\n".format(seq.id, start, end))

