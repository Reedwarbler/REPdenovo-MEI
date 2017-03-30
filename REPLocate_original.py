import sys
import os
import collections
from Bio import SeqIO
import operator

insert_size=300
derivation=50

dist1=insert_size-3*derivation
dist2=insert_size+3*derivation

def alpha_belt_hash():
    m_alpha={}
    m_alpha["a"]="0"
    m_alpha["A"]="0"
    m_alpha["c"]="1"
    m_alpha["C"]="1"
    m_alpha["g"]="2"
    m_alpha["G"]="2"
    m_alpha["t"]="3"
    m_alpha["T"]="3"
    m_alpha["n"]="4"
    m_alpha["N"]="4"
    return m_alpha

def preprocess_seq(m_alpha,seq):
    pseq=""
    for i in range(len(seq)):
        pseq=pseq+m_alpha[seq[i]]
    return pseq

def gnrt_reverse_complementary(s):
    lth=len(s)
    s_rc=""
    for i in range(lth-1,-1,-1):
        if s[i]=="0":
            s_rc=s_rc+"3"
        elif s[i]=="3":
            s_rc=s_rc+"0"
        elif s[i]=="1":
            s_rc=s_rc+"2"
        elif s[i]=="2":
            s_rc=s_rc+"1"
        else:
            s_rc=s_rc+"4"
    return s_rc

def is_qualified_clipped(cigar, cutoff_len):
    l=len(cigar)
    signal=[]
    lenth=[]
    temp=""
    for i in range(l):
        if cigar[i]>="0" and cigar[i]<="9":
            temp=temp+cigar[i]
        else:
            signal.append(cigar[i])
            try:
                lenth.append(int(temp))
            except ValueError:
                print "Error: ", cigar, temp
            temp=""

    b_qualified=False
    cnt_m=0
    clip_flag=0
    left_clip_len=0
    right_clip_len=0

    for i in range(len(signal)):
        if signal[i]=="M":
            cnt_m=cnt_m+lenth[i]

    if (signal[0]=="S" or signal[0]=="H") and lenth[0]>=cutoff_len:#left-clip
        clip_flag=1
        if signal[0]=="S":
            left_clip_len=lenth[0]
        b_qualified=True
    if (signal[len(signal)-1]=="S" or signal[len(signal)-1]=="H") and lenth[len(signal)-1]>=cutoff_len: #right-clip
        clip_flag=clip_flag+2 #if this is 3, then both side clipped
        if signal[len(signal)-1]=="S":
            right_clip_len=lenth[len(signal)-1]
        b_qualified=True

    return b_qualified, cnt_m, clip_flag, left_clip_len, right_clip_len

#count how many kmers hit the high frequence kmer library
def cnt_hit_kmer_lib(m_high_freq, seq, k):
    cnt=0
    for i in range(len(seq)-k+1):
        kmer=seq[i:i+k]
        if m_high_freq.has_key(kmer):
            cnt=cnt+1
            m_high_freq[kmer]=m_high_freq[kmer]+1 ########################################for temp usage
    return cnt


def cluster_sites(candidate_sites, n_extend):
    merged_sites={}
    marked={}
    for pos in candidate_sites:
        if marked.has_key(pos)==True:
            continue
        max_pos=pos
        sum_hit=candidate_sites[pos]
        for i in range(pos-n_extend, pos+n_extend+1):
            if candidate_sites.has_key(i) and i!=pos and marked.has_key(i)==False:
                sum_hit=sum_hit+candidate_sites[i]
                marked[i]=1
                if candidate_sites[i] > candidate_sites[max_pos]:
                    max_pos=i
        merged_sites[max_pos]=sum_hit

    return merged_sites


def cnt_discordant_support(potential_discordant_pe_pos, clip_pos):
    cnt=0
    for i in range(clip_pos-dist2, clip_pos+dist2+1):
        if potential_discordant_pe_pos.has_key(i):
            cnt=cnt+1
    return cnt


#first using clipped reads to get candidate sites
#then, check whether there are enough discordant pe support
def locate_rep_sites_clip(sf_high_freq_kmers, sf_simple_rep, anchor_mapq, pe_cutoff,
                          min_clip_lenth, k, hit_ratio, n_extend, clip_cutoff):
    m_alpha=alpha_belt_hash()
    m_high_freq={}
    #m_high_freq_test={} ############################################################################# temp usage
    with open(sf_high_freq_kmers) as fin_high_freq:
        s_cnt=0#####################################################################################
        for line in fin_high_freq:
            if line[0]==">":
                lfields=line.split("_")
                s_cnt=lfields[1]
                continue
            #pseq=preprocess_seq(m_alpha, line.rstrip())
            pseq=line.rstrip()#####################################################################################################################
            m_high_freq[pseq]=0
            #m_high_freq_test[pseq]=s_cnt #############################################################
            rc_seq=gnrt_reverse_complementary(pseq)
            m_high_freq[rc_seq]=0
            #m_high_freq_test[rc_seq]=s_cnt #############################################################

    m_sim_rep_kmer={}
    for record in SeqIO.parse(sf_simple_rep, "fasta"):
        sim_rep=str(record.seq)
        p_sim_rep=preprocess_seq(m_alpha, sim_rep)
        for i in range(len(p_sim_rep)-k+1):
            kmer=p_sim_rep[i:i+k]
            if m_sim_rep_kmer.has_key(kmer)==False:
                m_sim_rep_kmer[kmer]=1
            rc_kmer=gnrt_reverse_complementary(kmer)
            if m_sim_rep_kmer.has_key(rc_kmer)==False:
                m_sim_rep_kmer[rc_kmer]=1

    candidate_sites_left={}
    candidate_sites_right={}
    candidate_simple_left={}
    candidate_simple_right={}
    potential_discordant_pe_pos={}

    for sam_record in sys.stdin:
        sam_fields=sam_record.split()

        cigar=sam_fields[5]
        if cigar=="*":#unmapped reads
            continue

        map_pos=int(sam_fields[3])
        map_quality=int(sam_fields[4])
        ref=sam_fields[2]
        mate_ref=sam_fields[6]
        mate_pos=int(sam_fields[7])
        seq=sam_fields[9]

        ##save discordant reads information
        if mate_ref!="=" and mate_ref!="*" and map_quality >= anchor_mapq:
            if potential_discordant_pe_pos.has_key(map_pos)==False:
                potential_discordant_pe_pos[map_pos]=0
            potential_discordant_pe_pos[map_pos]=potential_discordant_pe_pos[map_pos]+1


        #check clipped reads
        bqualified_clip,cnt_m, clip_flag, left_clip_len, right_clip_len=is_qualified_clipped(cigar, min_clip_lenth)

        if left_clip_len>0:
            left_clip_seq= preprocess_seq(m_alpha, seq[:left_clip_len])
            n_left_hit=cnt_hit_kmer_lib(m_high_freq, left_clip_seq, k)
            n_left_all=left_clip_len-k+1
            if n_left_hit >= (n_left_all*hit_ratio):
                if candidate_sites_left.has_key(map_pos)==False:
                    candidate_sites_left[map_pos]=0
                candidate_sites_left[map_pos]=candidate_sites_left[map_pos]+1

            n_left_hit_simple=cnt_hit_kmer_lib(m_sim_rep_kmer, left_clip_seq,k)
            if n_left_hit_simple >= (n_left_all*hit_ratio):
                if candidate_simple_left.has_key(map_pos)==False:
                    candidate_simple_left[map_pos]=1
                candidate_simple_left[map_pos]=candidate_simple_left[map_pos]+1

        if right_clip_len>0:
            rd_len=len(seq)
            right_clip_seq=preprocess_seq(m_alpha, seq[rd_len-right_clip_len:])
            n_right_hit=cnt_hit_kmer_lib(m_high_freq, right_clip_seq, k)
            n_right_all=right_clip_len-k+1
            if n_right_hit >= (n_right_all*hit_ratio):
                clip_pos=map_pos+cnt_m-1
                if candidate_sites_right.has_key(clip_pos)==False:
                    candidate_sites_right[clip_pos]=0
                candidate_sites_right[clip_pos]=candidate_sites_right[clip_pos]+1

            n_right_hit_simple=cnt_hit_kmer_lib(m_sim_rep_kmer, right_clip_seq,k)
            if n_right_hit_simple >= (n_right_all*hit_ratio):
                if candidate_simple_right.has_key(map_pos)==False:
                    candidate_simple_right[map_pos]=1
                candidate_simple_right[map_pos]=candidate_simple_right[map_pos]+1

    merged_left=cluster_sites(candidate_sites_left, n_extend)
    merged_right=cluster_sites(candidate_sites_right, n_extend)

    candidate_sites={}
    candidate_site_simple_rep={}
    for pos in merged_left:
        if candidate_sites.has_key(pos)==False and merged_left[pos]>=clip_cutoff:
            n_discordant=cnt_discordant_support(potential_discordant_pe_pos, pos)
            if n_discordant >= pe_cutoff:
                candidate_sites[pos]=merged_left[pos]
                if candidate_simple_left.has_key(pos)==False:
                    print pos, merged_left[pos], n_discordant
                else:
                    candidate_site_simple_rep[pos]=merged_left[pos]

    for pos in merged_right:
        if candidate_sites.has_key(pos)==False and merged_right[pos]>=clip_cutoff:
            n_discordant=cnt_discordant_support(potential_discordant_pe_pos, pos)
            if n_discordant >= pe_cutoff:
                if candidate_sites.has_key(pos):
                    candidate_sites[pos]=candidate_sites[pos]+merged_right[pos]
                else:
                    candidate_sites[pos]=merged_right[pos]
                if candidate_simple_right.has_key(pos)==False:
                    print pos, merged_right[pos], n_discordant
                else:
                    candidate_site_simple_rep[pos]=candidate_sites[pos]

    print "-------------------------------------------------------------"

    for pos in candidate_site_simple_rep: ##simple repeats sites
        print pos, candidate_site_simple_rep[pos]


    #############################################################################for temp usage
    # m_sorted=sorted(m_high_freq.items(), key=operator.itemgetter(1), reverse=True)
    # with open("hit_kmer_statistic.txt", "w") as fout_sta:
    #     for (kmer, freq) in m_sorted:
    #         fout_sta.write(kmer+"\t"+str(freq)+"\t"+m_high_freq_test[kmer]+"\n")

    #return candidate_sites, candidate_site_simple_rep


if __name__ == "__main__":
    #anchor_mapq=10
    #sf_out="candidate_sites.txt"
    #parse_by_chrom(anchor_mapq,  sf_out)
    sf_high_freq_kmers=sys.argv[1]
    sf_simple_rep=sys.argv[2]
    min_clip_lenth=21
    k=21
    hit_ratio=float(sys.argv[3])
    n_extend=20
    clip_cutoff=int(sys.argv[4])
    anchor_mapq=20
    pe_cutoff=int(sys.argv[5])
    locate_rep_sites_clip(sf_high_freq_kmers, sf_simple_rep, anchor_mapq, pe_cutoff,
                          min_clip_lenth, k, hit_ratio, n_extend, clip_cutoff)
