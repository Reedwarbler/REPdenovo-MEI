import os
import sys

def is_clipped(cigar):
    cnt=0
    l=len(cigar)
    if cigar[l-1]=="S" or cigar[l-1]=="H":#right clipped
        cnt=2

    for i in range(l):
        if cigar[i]>="0" and cigar[i]<="9":
            continue
        else:
            if cigar[i]=="S" or cigar[i]=="H": ##left clipped
                cnt=cnt+1
            break
    return cnt


def is_hit(m_high_freq, seq, k):
    cnt=0
    for i in range(len(seq)-k+1):
        kmer=seq[i:i+k]
        if m_high_freq.has_key(kmer):
            cnt=cnt+1
    return cnt


def check_reads_in_kmer_lib(sf_high_freq_kmers, sf_reads, k):
    m_high_freq={}
    with open(sf_high_freq_kmers) as fin_high_freq:
        for line in fin_high_freq:
            if line[0]==">":
                continue
            m_high_freq[line.rstrip()]=1

    with open(sf_reads) as fin_algn:
        for line in fin_algn:
            fields=line.split()
            cigar=fields[5]
            seq=fields[9]
            if is_clipped(cigar)!=0:
                cnt=is_hit(m_high_freq, seq, k)
                if cnt>0:
                    print fields[0], cigar, cnt

if __name__ == "__main__":
    sf_high_freq_kmers=sys.argv[1]
    k=22
    sf_reads=sys.argv[2]
    check_reads_in_kmer_lib(sf_high_freq_kmers, sf_reads, k)
