
## modules
import os
import gzip
from Bio import SeqIO
from fastx.shared import NestedDefaultDict

def fastqc(infile, prefix):
    ## init
    bq_proportion_outfile = open("{}.bq_proportion.txt".format(prefix), "w")
    base_proportion_outfile = open("{}.base_proportion.txt".format(prefix), "w")
    seqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  ## header
            seq_id = j.strip()
        elif k == 1:  ## sequence
            seq = j.strip()
            seq_len = len(seq)
        elif k == 2:
            continue  ## plus
        elif k == 3:  ## quality
            seq_bq = j.strip()
            # outfile.write("{}\n{}\n+\n{}\n".format(seq_id, seq, seq_bq))
            # outfile.close()


    global dna
    global positions
    foo_hash = NestedDefaultDict()
    for k in positions:
        for v in dna:
            foo_hash[k][v] = 0

    for i, j in enumerate(open(foo)):
        k = i % 4
        if k == 0: ## header
            sequence_id = j.strip().replace("@" ,"")
        elif k == 1: ## sequence
            sequence = j.strip()
            sequence_list = list(sequence)
            sequence_length = len(sequence)
            sequence_position = np.asarray(range(0, sequence_length)) / sequence_length
            sequence_position = [format(m, ".2f") for m in sequence_position]

            ## append
            for m,n in zip(sequence_position, sequence_list):
                foo_hash[m][n] += 1

        elif k == 2: ## plus
            continue
        elif k == 3: ## quality scores
            continue

def fastqc(infile, prefix):
    if infile.endswith((".fq", ".fq.gz", ".fastq")):
        fastqc(infile, prefix)
