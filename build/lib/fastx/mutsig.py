import itertools
import natsort

# init
purine = ["A", "G"]
purine2pyrimidine = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

# mutational signatures
dna = list("atgc")
sig_lst = []
sig2sub = {}
sig2tri = {}
sub_lst = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
tri_lst = []
for sub in sub_lst:
    for upstream, downstream in itertools.product(dna, repeat=2):
        ref, alt = sub.split(">")
        sig = upstream + sub + downstream
        tri = upstream.upper() + ref + downstream.upper()
        sig2sub[sig] = sub
        sig2tri[sig] = tri
        sig_lst.append(sig)
        tri_lst.append(tri)
sig_lst = natsort.natsorted(sig_lst)
tri_lst = natsort.natsorted(list(set(tri_lst)))
