import os
import sys
import random
import itertools
from subprocess import *
from multiprocessing import Pool
from multiprocessing import Process, Manager

class KmerCluster:
    def __init__(self, sf_high_freq, sf_all_kmer, k, edit_distance, n_round):
        self.m_alpha=self.alpha_belt_hash()
        self.m_high_freq={}
        self.m_seed=[]
        self.edit_distance=edit_distance #allowed edit distance
        self.k=k #kmer length
        self.sf_high_freq=sf_high_freq
        self.random_colums=[]
        self.n_round=n_round
        self.sf_all_kmer=sf_all_kmer
        self.m_exclude_kmers=self.gnrt_exclusive_table()

    def read_in_kmer(self):
        m_kmer={}
        with open(self.sf_high_freq) as fin_high_freq:
            i_cnt=0
            for line in fin_high_freq:
                if line[0]==">":
                    lfields=line.split("_")
                    i_cnt=int(lfields[1].rstrip())
                    continue
                fields=line.split()
                pseq=self.preprocess_seq(self.m_alpha, fields[0].rstrip())
                m_kmer[pseq]=i_cnt
                rc_seq=self.gnrt_reverse_complementary(pseq)
                m_kmer[rc_seq]=i_cnt
        return m_kmer

    def alpha_belt_hash(self):
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
        m_alpha["0"]="A"
        m_alpha["1"]="C"
        m_alpha["2"]="G"
        m_alpha["3"]="T"
        return m_alpha

    def preprocess_seq(self, m_alpha,seq):
        pseq=""
        for i in range(len(seq)):
            pseq=pseq+m_alpha[seq[i]]
        return pseq

    def gnrt_reverse_complementary(self, s):
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

    def gnrt_random_column(self):
        random.seed()
        l_column=[]
        for i in range(self.k):
            l_column.append(-1)

        for i in range(self.edit_distance):
            #for j in range(i):
            while True:
                c=random.randint(0,self.k-1)
                if l_column[c]==-1:
                    l_column[c]=1
                    break

        l_sorted=[]
        for i in range(self.k):
            if l_column[i]==1:
                l_sorted.append(i)
        return l_sorted

    def gnrt_all_rounds_random_column(self, n_round):
        for i in range(n_round):
            l_column=self.gnrt_random_column()
            self.random_colums.append(l_column) ##save the column information

    #generate the exclusive kmer we do not want to involve in
    def gnrt_exclusive_table(self):
        m_exclusive={}
        l=self.k-self.edit_distance
        s1=""
        s2=""
        s3=""
        s4=""
        for i in range(l):
            s1=s1+"0"
            s2=s2+"3"
            if i%2==0:
                s3=s3+"0"
                s4=s4+"3"
            else:
                s3=s3+"3"
                s4=s4+"0"
        m_exclusive[s1]=1
        m_exclusive[s2]=1
        m_exclusive[s3]=1
        m_exclusive[s4]=1
        return m_exclusive


    def gnrt_seed(self, kmer, l_column):
        start=0
        seed=""
        for pos in l_column:
            seed=seed+kmer[start:pos]
            start=pos+1
        seed=seed+kmer[start:]
        return seed


    def one_n_round_select(self, m_kmers, i, m_exclude_kmers, b_save, m_high_freq):
        l_column=self.random_colums[i]
        m_temp={}
        for kmer in m_kmers:
            seed=self.gnrt_seed(kmer, l_column)
            if m_exclude_kmers.has_key(seed)==False:
                if b_save:
                    m_temp[seed]=True ###########################################
                if m_high_freq.has_key(kmer)==False:
                    m_high_freq[kmer]=m_kmers[kmer]
            else:
                m_high_freq[kmer]=-1
        self.m_seed.append(m_temp)

    def exclude_kmers_by_rounds(self, m_kmers, rounds, b_save):
        m_high_freq={}
        for i in range(rounds):
            self.one_n_round_select(m_kmers, i, self.m_exclude_kmers, b_save, m_high_freq)

        m_kmers.clear()#clear
        m_del={}
        for key in m_high_freq:
            if m_high_freq[key]<0:
                m_del[key]=True
        for key in m_del:
            del m_high_freq[key]
        return m_high_freq

    def dump_kmers(self, m_kmers, sf_out):
        with open(sf_out,"w") as fout_kmer:
            cnt=0
            for kmer in m_kmers:
                scontent=">{0}_{1}\n".format(cnt, m_kmers[kmer])
                fout_kmer.write(scontent+ kmer+"\n")
                cnt=cnt+1

    def select_high_freq_kmers(self, b_save):
        #given the high frequent kmers, remove AAAAAAAAAA, TTTTTT, ATATATAT..., TATATA...
        m_kmers=self.read_in_kmer()
        print "High frequent kmers are loaded in ..."
        self.gnrt_all_rounds_random_column(self.n_round)
        #print self.random_colums
        #return ########################################################################################################
        self.m_high_freq=self.exclude_kmers_by_rounds(m_kmers, self.n_round, b_save)
        #self.dump_kmers(self.m_high_freq, sf_selected)


    def run_kmer_cluster(self, sf_all_kmers, start, end, m_seed, sf_temp):
        with open(sf_all_kmers, "r") as fin_kmers:
            i_cnt=0
            with open(sf_temp,"w") as fout_temp:
                for line in itertools.islice(fin_kmers, start, end):
                    if line[0]==">":
                        lfields=line.split("_")
                        i_cnt=int(lfields[1].rstrip())
                        continue
                    kmer=self.preprocess_seq(self.m_alpha, line.rstrip())

                    for i in range(self.n_round):
                        l_column=self.random_colums[i]
                        seed=self.gnrt_seed(kmer, l_column)
                        if m_seed[i].has_key(seed):
                            fout_temp.write(kmer+" "+str(i_cnt)+"\n")
                            break

    def cluster_kmers_multi_thread(self, n_threads, sf_out):
        processes=[]
        n_line=0
        with open(self.sf_all_kmer) as f:#count number of lines of all the kmers
            for n_line, ltmp in enumerate(f):
                pass
        #print "Number of lines:{0}\n".format(n_line)
        fstart=0
        fend=0
        favrg=n_line/n_threads
        for i in range(n_threads):#split the lines
            fstart=fend
            fend=(i+1)*favrg
            if fend%2!=0:
                fend=fend+1
            if i==n_threads-1:
                fend=n_line+1
            sf_temp="tmp_thread_{0}_out.txt".format(i)
            p = Process(target=self.run_kmer_cluster, args=(self.sf_all_kmer, fstart, fend,
                                                            self.m_seed, sf_temp))
            p.start()
            processes.append(p)
        for p in processes:
            p.join()

        sf_01mer="tmp_01mer.txt"
        cmd="cat "
        for i in range(n_threads):
            cmd=cmd+"tmp_thread_{0}_out.txt ".format(i)
        cmd=cmd+" > {0}".format(sf_01mer)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        self.cvt_back_to_kmer(sf_01mer, sf_out)
        cmd="rm tmp_*"
        Popen(cmd, shell = True, stdout = PIPE).communicate()



    def cvt_01mer_kmer(self, s01):
        lenth=len(s01)
        kmer=""
        for i in range(lenth):
            kmer=kmer+self.m_alpha[s01[i]]
        return kmer

    def cvt_back_to_kmer(self, sf_01mer, sf_kmer):
        cnt=0
        with open(sf_kmer,"w") as fout_kmer:
            with open(sf_01mer) as fin_01mer:
                for line in fin_01mer:
                    fields=line.split()
                    s01=fields[0]
                    freq=fields[1]
                    kmer=self.cvt_01mer_kmer(s01)

                    fout_kmer.write(">0{0}_{1}\n".format(cnt,freq))
                    fout_kmer.write(kmer+"\n")
                    cnt=cnt+1


    def rm_tandem_repeat_kmer(self, sf_kmer):
        cmd="/data2/chongchu/tools/trf409.linux64 {0} 2 7 7 80 10 35 21".format(sf_kmer)

if __name__ == "__main__":
    sf_high_freq=sys.argv[1]
    sf_all_kmer=sys.argv[2]
    k=21
    edit_distance=1
    #n_round=30
    n_round=int(sys.argv[3])
    #sf_select_high_freq=sys.argv[4]
    n_threads=4
    sf_parsed=sys.argv[4]
    #sf_parsed="new_involved_kmers_of_low_freq.fa"

    kc=KmerCluster(sf_high_freq, sf_all_kmer, k, edit_distance, n_round)
    kc.select_high_freq_kmers(True)
    kc.cluster_kmers_multi_thread(n_threads,sf_parsed)

