import os
import pyfastx


def build_index(infile):
    fxifile = infile + ".fxi"
    if os.path.exists(fxifile):
        print("fxi index is present")
    else:
        print("buliding fxi index for {}".format(infile))
        if infile.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            pyfastx.Fasta(infile)
        else:
            pyfastx.Fastq(infile)
        print("fxi index has been created for {}".format(infile))


def seq_index(infile):
    if infile.endswith(
        (".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fq", ".fq.gz", ".fastq", ".fastq.gz")
    ):
        build_index(infile)
    else:
        print("fastx index doesn't support the provided input file")
