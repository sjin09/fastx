
## modules
import os
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

## global
dna = list("ATGC")
positions = [format(x, ".2f") for x in np.arange(0.0, 1.01, 0.01)]

class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))

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

    bar = open(foo.replace(".fastq",".per_position_sequence_content"), "w")
    bar.write("{}\t{}\t{}\t{}\t{}\n".format("position", "A", "T", "G", "C"))
    for k in foo_hash:
        base_count = np.asarray(list(foo_hash[k].values()))
        base_sum = sum(foo_hash[k].values())
        base_proportion = [format(v, ".3f") for v in base_count/base_sum]
        bar.write("{}\n".format("\t".join([k] + base_proportion)))
    bar.close()


def fasta_split(infile, outdir):
    fasta = (
        SeqIO.parse(infile, "fasta")
        if infile.endswith((".fa", ".fasta"))
        else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    )
    for seq in fasta:
        outfile = open(
            os.path.join(outdir, "{}.fa".format(seq.id.replace("/", "_"))), "w"
        )
        outfile.write(">{}\n".format(seq.id))
        for chunk in chunkstring(seq.seq, 50):
            outfile.write("{}\n".format(chunk))
        outfile.close()


def fastq_split(infile, outdir):
    seqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  ## header
            seq_id = j.strip()
            outfile = open(
                os.path.join(
                    outdir, "{}.fq".format(seq_id.replace("@", "").replace("/", "_"))
                ),
                "w",
            )
        elif k == 1:  ## sequence
            seq = j.strip()
            seq_len = len(seq)
        elif k == 2:
            continue  ## plus
        elif k == 3:  ## quality
            seq_bq = j.strip()
            outfile.write("{}\n{}\n+\n{}\n".format(seq_id, seq, seq_bq))
            outfile.close()


def seq_split(infile, outdir):
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if infile.endswith((".fa", ".fa.gz", ".fasta")):
        fasta_split(infile, outdir)
    elif infile.endswith((".fq", ".fq.gz", ".fastq")):
        fastq_split(infile, outdir)
