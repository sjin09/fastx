# modules
import gzip
import natsort
from Bio import SeqIO
from collections import defaultdict


def chunkstring(string):
    chunks = [
        string[i : i + 60] for i in range(0, len(string), 60)
    ]
    return chunks


def fasta_sort(infile, outfile):
    counter = 0
    seq_hash = defaultdict()
    fasta = (
        SeqIO.parse(infile, "fasta")
        if infile.endswith((".fa", ".fasta"))
        else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    )
    for seq in fasta:
        seq_hash[seq.id] = seq.seq  # assumes
        counter += 1
    seq_hash_count = len(seq_hash.keys())

    if seq_hash_count == counter:
        unsorted_keys = seq_hash.keys()
        sorted_keys = natsort.natsorted(seq_hash.keys())
        k = all(sorted_keys[i] == j for i, j in enumerate(unsorted_keys))
        if k:
            print("Sequences are already sorted")
        else:
            for seq_id in sorted_keys:
                outfile.write(">{}\n".format(seq_id))
                for chunk in chunkstring(seq_hash[seq_id]):
                    outfile.write("{}\n".format(chunk))
    else:
        print("Sequences do not have unique IDs")


def fastq_sort(infile, outfile):
    counter = 0
    seq_hash = defaultdict()
    seqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  # header
            seq_id = j.strip()
        elif k == 1:  # sequence
            seq = j.strip()
        elif k == 2:
            continue  # plus
        elif k == 3:  # quality
            counter += 1
            seq_bq = j.strip()
            seq_hash[seq_id] = "{}\n{}\n+\n{}\n".format(seq_id, seq, seq_bq)
    # return
    seq_hash_count = len(seq_hash.keys())
    if seq_hash_count == counter:
        unsorted_keys = seq_hash.keys()
        sorted_keys = natsort.natsorted(seq_hash.keys())
        if all(sorted_keys[i] == j for i, j in enumerate(unsorted_keys)):
            print("Sequences are already sorted")
        else:
            for seq_id in sorted_keys:
                outfile.write("{}".format(seq_hash[seq_id]))
    else:
        print("Sequences do not have unique IDs")


def seq_sort(infile, outfile):
    if infile.endswith((".fa", ".fa.gz", ".fasta")):
        fasta_sort(infile, outfile)
    elif infile.endswith((".fq", ".fq.gz", ".fastq")):
        fastq_sort(infile, outfile)
