import gzip
import natsort
from Bio import SeqIO


def get_target(infile: str, outfile) -> None:
    contig_hsh = {}
    SEX_CHROMOSOME = ("W", "X", "Y", "Z")
    sequences = SeqIO.parse(infile, "fasta") if infile.endswith(".fasta") else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    for i in sequences:
        contig = i.id
        contig_length = len(i.seq)
        if (contig.startswith("Super") or contig.startswith("SUPER")) and not "unloc" in contig and not contig.endswith(SEX_CHROMOSOME) and contig_length > 1000000:
            contig_hsh[contig] = contig_length
    contig_lst = natsort.natsorted(list(contig_hsh.keys()))
    for contig in contig_lst:
        outfile.write("{}\n".format(contig))

