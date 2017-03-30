import os
import sys
from kmer_cluster import KmerCluster

sf_high_freq=sys.argv[1]
sf_all_kmer=sys.argv[2]
k=21
edit_distance=1
n_round=int(sys.argv[3])
sf_select_high_freq=sys.argv[4]

kmcluster=KmerCluster(sf_high_freq, sf_all_kmer, k, edit_distance, n_round)
kmcluster.select_high_freq_kmers(True, sf_select_high_freq)
kmcluster.cluster_kmers()
