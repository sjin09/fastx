## modules
import os
import gzip
import pyfastx
from Bio import SeqIO
from collections import defaultdict

# global
INFILE_SUFFIX = (
    ".fa",
    ".fq",
    ".fa.gz",
    ".fq.gz",
    ".fasta",
    ".fastq",
    ".fasta.gz",
    ".fastq.gz",
)
FASTA_SUFFIX = (".fa", ".fa.gz", ".fasta", ".fasta.gz")
FASTQ_SUFFIX = (".fa", ".fq.gz", ".fastq", ".fastq.gz")

class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))


def chunkstring(string):
    chunks = [string[i : i + 60] for i in range(0, len(string), 60)]
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
                seq_id = (
                    seq_id if not isinstance(seq_id, bytes) else seq_id.decode("utf-8")
                )
                seq_id = seq_id.replace("@", ">")
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
    return_whitelist(infile, blacklist, outfile)
