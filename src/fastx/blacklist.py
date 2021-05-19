import os
import gzip
import pyfastx
import natsort

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


def load_blacklist(blacklist):
    blacklist_lst = [line.strip() for line in open(blacklist).readlines()]
    return blacklist_lst


def load_whitelist(zmw_lst, blacklist_set):
    zmw_set = set(zmw_lst)
    whitelist_lst = list(zmw_set.difference(blacklist_set))
    whitelist_lst = natsort.natsorted(whitelist_lst)
    return whitelist_lst


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


def return_whitelist(infile, blacklist, outfile):
    blacklist_set = set(load_blacklist(blacklist))
    seqfile = load_seqfile(infile)
    zmw_lst = seqfile.keys()
    whitelist_lst = load_whitelist(zmw_lst, blacklist_set)
    for zmw in whitelist_lst:
        read = seqfile[zmw]
        outfile.write("{}".format(read.raw))


def filter_blacklist(infile, blacklist, outfile):
    if infile.endswith(INFILE_SUFFIX):
        return_whitelist(infile, blacklist, outfile)
    else:
        print("fastx blacklist doesn't support the provided input file")
