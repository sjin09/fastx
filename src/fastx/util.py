## modules
import gzip
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))


def chunkstring(string):
    chunks = [
        string[i : i + 60] for i in range(0, len(string), 60)
    ]
    return chunks


def fasta2fastq(infile, outfile):
    if infile.endswith((".fa", ".fa.gz", ".fasta")):
        fasta = (
            SeqIO.parse(infile, "fasta")
            if infile.endswith((".fa", ".fasta"))
            else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
        )
        for seq in fasta:
            outfile.write(
                "@{}\n{}\n+\n{}\n".format(seq.id, seq.seq, "!" * len(seq.seq))
            )
    else:
        print("Did you provide a FASTA file?")


def fastq2fasta(infile, outfile):

    if infile.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
        seqfile = (
            open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
        )
        for i, j in enumerate(seqfile):
            k = i % 4
            if k == 0:  ## header
                seq_id = j.strip()
                seq_id = seq_id if not isinstance(seq_id, bytes) else seq_id.decode("utf-8")
                seq_id.replace("@", ">")
            elif k == 1:  ## sequence
                seq = j.strip()
                seq = seq if not isinstance(seq, bytes) else seq.decode("utf-8")
            elif k == 2:
                continue  ## plus
            elif k == 3:  ## quality
                outfile.write("{}\n".format(seq_id))
                for chunk in chunkstring(seq):
                    outfile.write("{}\n".format(chunk))
    else:
        print("Did you provide a FASTQ file?")
