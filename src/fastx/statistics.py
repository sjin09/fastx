## modules
import gzip
from Bio import SeqIO


def length_list(infile):
    lst = []
    if infile.endswith((".fa", ".fa.gz", ".fasta")):
        fasta = (
            SeqIO.parse(infile, "fasta")
            if infile.endswith((".fa", ".fasta"))
            else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
        )
        for seq in fasta:
            lst.append(len(seq.seq))
    elif infile.endswith((".fq", ".fq.gz", ".fastq")):
        seqfile = (
            open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
        )
        for i, j in enumerate(seqfile):
            k = i % 4
            if k == 0:  ## header
                continue
            elif k == 1:  ## sequence
                lst.append(len(j.strip()))
            elif k == 2:
                continue  ## plus
            elif k == 3:  ## quality
                continue
    return lst


def calculate_n50(length):
    total = sum(length)
    cumulative_total = 0
    for l in length:
        cumulative_total += l
        if cumulative_total >= total / 2:
            return l


def seq_statistics(infile, outfile):
    lst = length_list(infile)
    total = sum(lst)
    count = len(lst)
    min_length = min(lst)
    max_length = max(lst)
    mean_length = total / count
    n50 = calculate_n50(sorted(lst, reverse=True))
    ## return
    outfile.write("infile: {}".format(infile))
    outfile.write("number_of_seuqences: {}".format(count))
    outfile.write("N50: {}".format(n50))
    outfile.write("min_length: {}".format(min_length))
    outfile.write("mean_length: {}".format(mean_length))
    outfile.write("max_length: {}".format(max_length))
    outfile.write("total (bp): {}".format(total))
