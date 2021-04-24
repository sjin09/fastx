import gzip
from Bio import SeqIO
from fastx.mutsig import purine, purine2pyrimidine, sig_lst, tri_lst


def fa_tricounts(infile: str, threshold: int, outfile):
    # init
    tri2count = {_tri: 0 for _tri in tri_lst}
    seqfile = SeqIO.parse(infile, "fasta") if infile.endswith(".fasta") else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    
    for i in seqfile:
        seq = str(i.seq)
        slen = len(seq)
        if slen < threshold: continue
        seq_tri_lst = [seq[sdex:sdex+3] for sdex, _base in enumerate(seq)]
        # seq_tri_lst = [_ for _ in seq_tri_lst if len(_) == 3 and _.count("N") == 0]
        for _tri in seq_tri_lst:
            if _tri.count ("N") != 0:
                print(_tri)
            elif len(_tri) != 3:
                print(_tri)
        # for tri in tri_list:
        # ref = tri[1]
        # if ref in purine:
        #     tri = "".join([purine2pyrimidine[base] for base in tri[::-1]])
        #     tri2count[tri] += 1
        # else:
        #     tri2count[tri] += 1

        
def tricounts(infile, threshold, outfile):
    if infile.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
        fa_tricounts(infile, threshold, outfile)
    # elif infile.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
    #     fq_tricounts(infile, outfile)
    # else:
    #     print("tri_counts doesn't support the provided input")


        ## main
        # for j in sequences:
        #     seq, slen = str(j.seq), len(j.seq)
        #     if slen < threshold: continue
        #     tri_list = [seq[index:index+3] for index, _ in enumerate(seq)]
        #     tri_list = [_ for _ in tri_list if len(_) == 3 and _.count("N") == 0]
        #     for tri in tri_list:
        #         ref = tri[1]
        #         if ref in purine:
        #             tri = "".join([purine2pyrimidine[base] for base in tri[::-1]])
        #             tri2count[tri] += 1
        #         else:
        #             tri2count[tri] += 1
        # for _tri in tri_lst:
        #     outfile.write("{}\t{}\n".format(_tri, tri2count[_tri]))
