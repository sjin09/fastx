
## modules
import os
import gzip
import numpy as np
from fastx.common import NestedDefaultDict

## global
dna = list("ATGC")
positions = [format(x, ".2f") for x in np.arange(0.0, 1.01, 0.01)]

def fastqc(infile, prefix):
    ## init
    global dna
    global positions
    bq_count_hash = {i:0 for i in range(1,94)}
    bq_cumsum_hash = {}
    base_proportion_hash = NestedDefaultDict()
    bq_proportion_outfile = open("{}.bq_proportion.txt".format(prefix), "w")
    base_proportion_outfile = open("{}.base_proportion.txt".format(prefix), "w")
    for p in positions:
        for base in dna:
            base_proportion_hash[p][base] = 0

    seq_sum = 0 
    seqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  ## header
            seq_id = j.strip()
        elif k == 1:  ## sequence
            seq = j.strip()
            seq_len = len(seq)
            seq_lst = list(seq)
            seq_sum += seq_len
            seq_pos = np.asarray(range(0, seq_len)) / seq_len
            seq_pos = [format(m, ".2f") for m in seq_pos]
            ## append
            for m, n in zip(seq_pos, seq_lst):
                base_proportion_hash[m][n] += 1
        elif k == 2:
            continue  ## plus
        elif k == 3:  ## quality
            bq_ascii = j.strip()
            bq_int_lst = [ord(bq) - 33 for bq in list(bq_ascii)]
            for bq in bq_int_lst:
                bq_count_hash[bq] += 1

    ## return
    base_proportion_outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("POS", "A", "T", "G", "C", "COUNTS"))
    for pos in base_proportion_hash:
        # print(pos, base_proportion_hash[pos], sum(base_proportion_hash[pos].values()))
        base_count = base_proportion_hash[pos].values()
        total = sum(base_count)
        base_proportion = [format(count/total, ".3f") for count in base_count]
        base_count = [":".join([str(_count) for _count in base_count])]
        base_proportion_outfile.write("{}\n".format("\t".join([pos] + base_proportion + base_count)))
    base_proportion_outfile.close()

    ## return
    for upper_bq in list(range(10, 100, 10)):
        count = 0
        for _bq in range(1, upper_bq):
            count += bq_count_hash[_bq]
        bq_cumsum_hash[upper_bq] = count
    bq_cumsum_hash[93] = bq_count_hash[93]
    for _bq in bq_cumsum_hash:
        bq_proportion_outfile.write("Q{}: {}/{} = {:.2f}%\n".format(_bq, bq_cumsum_hash[_bq], seq_sum, (bq_cumsum_hash[_bq]/seq_sum) * 100))
    


def seq_fastqc(infile, prefix):
    if infile.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
        fastqc(infile, prefix)
    else:
        print("fastqc does not support the provided input file")
