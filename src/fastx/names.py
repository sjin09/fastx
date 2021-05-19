import os
import pyfastx

FASTA_SUFFIX = (".fa", ".fa.gz", ".fasta", ".fasta.gz")
FASTQ_SUFFIX = (".fq", ".fq.gz", ".fastq", ".fastq.gz")
INFILE_SUFFIX = (
    ".fa",
    ".fa.gz",
    ".fasta",
    ".fasta.gz",
    ".fq",
    ".fq.gz",
    ".fastq",
    ".fastq.gz",
)

def load_seqfile(infile):
    fxifile = infile + ".fxi"
    if os.path.exists(fxifile) and infile.endswith(FASTA_SUFFIX):
        seqfile = pyfastx.Fasta(infile, build_index=False)
    elif not os.path.exists(fxifile) and infile.endswith(FASTA_SUFFIX):
        seqfile = pyfastx.Fasta(infile, build_index=True)
    elif os.path.exists(fxifile) and infile.endswith(FASTQ_SUFFIX):
        seqfile = pyfastx.Fastq(infile, build_index=False)
    elif not os.path.exists(fxifile) and infile.endswith(FASTQ_SUFFIX):
        seqfile = pyfastx.Fastq(infile, build_index=True)
    return seqfile


def return_names(infile, outfile):
    seqfile = load_seqfile(infile)
    id_lst = seqfile.keys()
    for id in id_lst:
        outfile.write("{}\n".format(id))


def seq_names(infile, outfile):
    if infile.endswith(INFILE_SUFFIX):
        return_names(infile, outfile)
    else:
        print("fastx blacklist doesn't support the provided input file")
