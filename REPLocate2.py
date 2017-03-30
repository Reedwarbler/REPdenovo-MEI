import sys
import os
import collections
from Bio import SeqIO
import operator
import getopt


class REPLocator:
    def __init__(self, inzt, std):
        self.insert_size=inzt
        self.derivation=std
        self.dist2=self.insert_size+3*self.derivation

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
        return m_alpha

    def preprocess_seq(self, m_alpha,seq):
        pseq=""
        for i in range(len(seq)):
            pseq=pseq+m_alpha[seq[i]]
        return pseq

    def gnrt_reverse_complementary(self,s):
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

    def is_qualified_clipped(self, cigar, cutoff_len):
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

        b_hardclip=False
        cnt_m=0
        clip_flag=0
        left_clip_len=0
        right_clip_len=0

        for i in range(len(signal)):
            if signal[i]=="M":
                cnt_m=cnt_m+lenth[i]

        if (signal[0]=="S" or signal[0]=="H") and lenth[0]>=cutoff_len:#left-clip
            clip_flag=1
            left_clip_len=lenth[0]
            if signal[0]=="H":
                b_hardclip=True
        if (signal[len(signal)-1]=="S" or signal[len(signal)-1]=="H") and lenth[len(signal)-1]>=cutoff_len: #right-clip
            clip_flag=clip_flag+2 #if this is 3, then both side clipped
            right_clip_len=lenth[len(signal)-1]
            if signal[len(signal)-1]=="H":
                b_hardclip=True
        return b_hardclip, cnt_m, clip_flag, left_clip_len, right_clip_len

    #count how many kmers hit the high frequency kmer library
    def cnt_hit_kmer_lib(self, m_high_freq, seq, k):
        cnt=0
        for i in range(len(seq)-k+1):
            kmer=seq[i:i+k]
            if m_high_freq.has_key(kmer):
                cnt=cnt+1
                #m_high_freq[kmer]=m_high_freq[kmer]+1 ########################################for temp usage
        return cnt

    def cluster_sites(self, candidate_sites, n_extend):
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

    def cnt_clip_in_range(self, m_soft_hard_pos, pos, n_extend):
        n_clip=0
        for i in range(pos-n_extend, pos+n_extend+1):
            if m_soft_hard_pos.has_key(i)==True:
                n_clip=n_clip+m_soft_hard_pos[i]
        return n_clip

    def cnt_discordant_support(self, potential_discordant_pe_pos, clip_pos, dist):
        cnt=0
        for i in range(clip_pos-dist, clip_pos+dist+1):
            if potential_discordant_pe_pos.has_key(i):
                cnt=cnt+1
        return cnt

    def load_high_freq_kmers(self, sf_high_freq_kmers):
        m_alpha=self.alpha_belt_hash()
        m_high_freq={}

        print "Loading high frequent kmers into memory ..."
        with open(sf_high_freq_kmers) as fin_high_freq:
            #s_cnt=0#####################################################################################
            for line in fin_high_freq:
                if line[0]==">":
                    #lfields=line.split("_")
                    #s_cnt=lfields[1]
                    continue
                pseq=self.preprocess_seq(m_alpha, line.rstrip())
                m_high_freq[pseq]=0
                rc_seq=self.gnrt_reverse_complementary(pseq)
                m_high_freq[rc_seq]=0
        return m_high_freq, m_alpha
    
    #collect candidate clipped reads
    #collect candidate 01pair reads
    def collect_candidate_clipped_01pair_reads(self, sf_high_freq_kmers, sf_simple_rep, anchor_mapq,
                              min_clip_lenth, k, hit_ratio, rlenth, sf_tmp_clip_reads, sf_tmp_01_pair, sf_tmp_discordt):
        m_high_freq, m_alpha=self.load_high_freq_kmers(sf_high_freq_kmers)
        print "Processing alignments ..."
        with open(sf_tmp_clip_reads,"w") as fout_clip, open(sf_tmp_01_pair,"w") as fout_01pair, \
                open(sf_tmp_discordt,"w") as fout_discordt:
            for sam_record in sys.stdin:
                sam_fields=sam_record.split()
                cigar=sam_fields[5]
                map_quality=int(sam_fields[4])
                map_pos=int(sam_fields[3])
                qname=sam_fields[0]
                flag=int(sam_fields[1])
                bfirst=1#first in the pair
                if flag&64==0:
                    bfirst=0
                seq=sam_fields[9]
                rnext=sam_fields[6]
                pnext=int(sam_fields[7])

                if flag&4!=0 and flag&8==0 and rnext=="=":#read unmapped but mate mapped
                    read_seq=self.preprocess_seq(m_alpha, seq)#check the unmapped sequence, whether is high frequency sequence
                    n_hit=self.cnt_hit_kmer_lib(m_high_freq, read_seq, k)
                    n_all=len(seq)-k+1
                    if n_hit>=(n_all*hit_ratio):
                        fout_01pair.write(qname+" "+str(bfirst)+" "+str(pnext)+"\n")

                if cigar=="*":
                    continue

                #check clipped reads
                b_hardclip, cnt_m, clip_flag, left_clip_len, right_clip_len=self.is_qualified_clipped(cigar, min_clip_lenth)

                if left_clip_len>0:
                    if b_hardclip==True:
                        fout_clip.write(qname+" "+str(map_pos)+" "+str(map_quality)+" "+str(bfirst)+" 2"+"\n")
                        #print cigar, "#########" ######################################################################################33
                    else:
                        left_clip_seq= self.preprocess_seq(m_alpha, seq[:left_clip_len])
                        n_left_hit=self.cnt_hit_kmer_lib(m_high_freq, left_clip_seq, k)
                        n_left_all=left_clip_len-k+1
                        if n_left_hit >= (n_left_all*hit_ratio):
                            fout_clip.write(qname+" "+str(map_pos)+" "+str(map_quality)+" "+str(bfirst)+" 1"+"\n")#save read_id, clip_pos, mapq, first-or-second-in-pair
                        else:
                            fout_clip.write(qname+" "+str(map_pos)+" "+str(map_quality)+" "+str(bfirst)+" 0"+"\n")#save the unsupported ones

                if right_clip_len>0:
                    clip_pos=map_pos+cnt_m-1
                    if b_hardclip==True:
                        fout_clip.write(qname+" "+str(clip_pos)+" "+str(map_quality)+" "+str(bfirst)+" 2"+"\n")
                        #print cigar, "************" ######################################################################################33
                    else:
                        right_clip_seq=self.preprocess_seq(m_alpha, seq[rlenth-right_clip_len:])
                        n_right_hit=self.cnt_hit_kmer_lib(m_high_freq, right_clip_seq, k)
                        n_right_all=right_clip_len-k+1
                        if n_right_hit >= (n_right_all*hit_ratio):
                            fout_clip.write(qname+" "+str(clip_pos)+" "+str(map_quality)+" "+str(bfirst)+" 1"+"\n")
                        else:
                            fout_clip.write(qname+" "+str(clip_pos)+" "+str(map_quality)+" "+str(bfirst)+" 0"+"\n")

                #collect reads that map to other chromosomes
                if rnext!="=" and rnext!="*" and map_quality >= anchor_mapq:
                    fout_discordt.write(qname+" "+str(bfirst)+" "+str(map_pos)+" "+rnext+" "+str(pnext)+"\n") #


    #first using clipped reads to get candidate sites
    #then, check whether there are enough discordant pe support
    def locate_rep_sites(self, sf_high_freq_kmers, sf_clipped_reads, sf_01_pair, sf_discord, sf_gap_pos, chrom,
                         anchor_mapq, k, hit_ratio, n_extend, soft_clip_cutoff, all_clip_cutoff,
                         sf_discord_valid, sf_candidate_only_clip):

        m_high_freq, m_alpha=self.load_high_freq_kmers(sf_high_freq_kmers)

        m_clipped_reads={}
        m_unsupport_clip_pos={}
        m_hard_clip_reads={}

        with open(sf_clipped_reads) as fin_clipped_reads:
            for line in fin_clipped_reads:
                fields=line.split()
                read_id=fields[0]
                pos=int(fields[1])
                mapq=int(fields[2])
                bfirst=int(fields[3])
                bflag=fields[4]
                if bflag=="1":
                    m_clipped_reads[read_id]=(pos, mapq, bfirst)
                elif bflag=="2":
                    m_hard_clip_reads[read_id]=(pos, mapq, bfirst)
                else:
                    if m_unsupport_clip_pos.has_key(pos)==False:
                        m_unsupport_clip_pos[pos]=1
                    else:
                        m_unsupport_clip_pos[pos]=m_unsupport_clip_pos[pos]+1

        m_01pair={}
        with open(sf_01_pair) as fin_01pair:
            for line in fin_01pair:
                fields=line.split()
                qname=fields[0]
                bfirst=int(fields[1])
                pos=int(fields[2])
                m_01pair[qname]=(bfirst, pos)

        m_discord_need_validate={}
        with open(sf_discord) as fin_discord:
            for line in fin_discord:
                fields=line.split() #read_id, bfirst, rnext, pnext
                m_discord_need_validate[fields[0]]=(fields[1], fields[2], fields[3])

        potential_discordant_pe_pos={}
        potential_clip_pos={}
        potential_hard_clip_pos={}
        potential_map_mate_unmap_pos={}
        with open(sf_discord_valid, "w") as fout_valid_pe:
            for sam_record in sys.stdin:
                sam_fields=sam_record.split()

                cigar=sam_fields[5]
                map_quality=int(sam_fields[4])
                if cigar=="*":#unmapped reads
                    continue

                qname=sam_fields[0]
                flag=int(sam_fields[1])
                bfirst=1#first in the pair
                if flag&64==0:
                    bfirst=0
                read_seq=sam_fields[9]
                map_pos=int(sam_fields[3])
                mate_ref=sam_fields[6]

                #filter clipped reads by using the anchor's mapping quality
                if m_clipped_reads.has_key(qname) and bfirst!=m_clipped_reads[qname][2] and map_quality>=anchor_mapq:
                    clip_pos=int(m_clipped_reads[qname][0])
                    if potential_clip_pos.has_key(clip_pos)==False:
                        potential_clip_pos[clip_pos]=1
                    else:
                        potential_clip_pos[clip_pos]=potential_clip_pos[clip_pos]+1
                elif m_hard_clip_reads.has_key(qname) and bfirst!=m_hard_clip_reads[qname][2] and map_quality>=anchor_mapq:
                    clip_pos=int(m_hard_clip_reads[qname][0])
                    if potential_hard_clip_pos.has_key(clip_pos)==False:
                        potential_hard_clip_pos[clip_pos]=1
                    else:
                        potential_hard_clip_pos[clip_pos]=potential_hard_clip_pos[clip_pos]+1

                #filter 01-pair
                if m_01pair.has_key(qname) and bfirst!=m_01pair[qname][0] and map_quality>=anchor_mapq:
                    pos=m_01pair[qname][1]
                    if potential_map_mate_unmap_pos.has_key(pos)==False:
                        potential_map_mate_unmap_pos[pos]=1
                    else:
                        potential_map_mate_unmap_pos[pos]=potential_map_mate_unmap_pos[pos]+1

                # ##save discordant reads information
                # if mate_ref!="=" and mate_ref!="*" and map_quality >= anchor_mapq:
                #     if potential_discordant_pe_pos.has_key(map_pos)==False:
                #         potential_discordant_pe_pos[map_pos]=1
                #     else:
                #         potential_discordant_pe_pos[map_pos]=potential_discordant_pe_pos[map_pos]+1

                #filter discordant pair
                if m_discord_need_validate.has_key(qname) and bfirst!=int(m_discord_need_validate[qname][0]):
                    seq=self.preprocess_seq(m_alpha, read_seq)
                    n_hit=self.cnt_hit_kmer_lib(m_high_freq, seq, k)
                    n_all=len(seq)-k+1
                    valid_chrm=m_discord_need_validate[qname][1]
                    valid_pos=m_discord_need_validate[qname][2]
                    if n_hit >= (n_all*hit_ratio):
                        fout_valid_pe.write(valid_chrm+" "+valid_pos+" "+"1\n")
                    else:
                        fout_valid_pe.write(valid_chrm+" "+valid_pos+" "+"0\n")
            #insert size changed pairs (become shorter)

        #merge the neighbored clip positions
        merged_clip_pos=self.cluster_sites(potential_clip_pos, n_extend)
        #merged_hard_clip_pos=self.cluster_sites(potential_hard_clip_pos, n_extend)

        potential_soft_hard_clip_pos={}
        for pos in potential_clip_pos:
            n_soft=potential_clip_pos[pos]
            n_hard=0
            if potential_hard_clip_pos.has_key(pos):
                n_hard=potential_hard_clip_pos[pos]
            n_all=n_soft+n_hard
            potential_soft_hard_clip_pos[pos]=n_all

        for pos in potential_hard_clip_pos:
            if potential_soft_hard_clip_pos.has_key(pos)==False:
                potential_soft_hard_clip_pos[pos]=potential_hard_clip_pos[pos]

        #get all the gap positions
        m_gap_pos={}
        with open(sf_gap_pos) as fin_gap_pos:
            for line in fin_gap_pos:
                fields=line.split()
                chrm=fields[3]
                pstart=int(fields[0])
                pend=int(fields[1])
                if m_gap_pos.has_key(chrm)==False:
                    m_gap_pos[chrm]={}
                else:
                    m_gap_pos[chrm][pstart]=1
                    m_gap_pos[chrm][pend]=1

        #call the sites according to cutoff
        with open(sf_candidate_only_clip,"w") as fout_candidate:
            for pos in merged_clip_pos: #print pos, merged_left[pos], n_discordant, "left"
                if merged_clip_pos[pos]<soft_clip_cutoff: ##first satisfy the clip cutoff
                    continue
                n_all_clip=self.cnt_clip_in_range(potential_soft_hard_clip_pos, pos, n_extend)
                if n_all_clip<all_clip_cutoff:
                    continue
                #check whether is a gap, if so, then continue
                n_gap_support=self.cnt_discordant_support(m_gap_pos[chrom], pos, rlenth)
                if n_gap_support>0:
                    continue

                # #check discordant pe support
                # n_discordant=self.cnt_discordant_support(potential_discordant_pe_pos, pos, self.dist2)
                n_01pair=self.cnt_discordant_support(potential_map_mate_unmap_pos, pos, self.dist2)
                #
                # if (n_discordant+n_01pair) < pe_cutoff: ##next satisfy the pe cutoff
                #     continue

                n_unsupport=0
                n_hard_clip=0
                for ipos in range(pos-n_extend, pos+n_extend):
                    if m_unsupport_clip_pos.has_key(ipos):
                        n_unsupport=n_unsupport+m_unsupport_clip_pos[ipos]
                    if potential_hard_clip_pos.has_key(ipos):
                        n_hard_clip=n_hard_clip+potential_hard_clip_pos[ipos]

                #print pos, merged_clip_pos[pos], n_unsupport, n_hard_clip, n_discordant, n_01pair
                #print pos, merged_clip_pos[pos], n_unsupport, n_hard_clip
                fout_candidate.write(str(pos)+" "+str(merged_clip_pos[pos])+" "+str(n_unsupport)
                                     +" "+str(n_hard_clip)+" "+str(n_01pair)+"\n")


    def call_final_sites(self, sf_fai, sf_out_folder, pe_cutoff):
        m_chrm_discord_pos={}
        m_chrm_discord_pos_unsupt={}
        chrm_list=[]
        with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                chrm=fields[0]
                chrm_list.append(chrm)
                sf_discord_valided="{0}{1}_pe_validated.tmp".format(sf_out_folder, chrm)#PEs (of all chroms) validated from chrm
                with open(sf_discord_valided) as fin_discord:
                    for record in fin_discord:
                        flds=record.split() ##chrm position
                        supt_chrm=flds[0]
                        supt_pos=int(flds[1])
                        b_valid=int(flds[2])

                        if b_valid==1:
                            if m_chrm_discord_pos.has_key(supt_chrm)==False:
                                m_chrm_discord_pos[supt_chrm]={}
                            if m_chrm_discord_pos[supt_chrm].has_key(supt_pos)==False:
                                m_chrm_discord_pos[supt_chrm][supt_pos]=1
                            else:
                                m_chrm_discord_pos[supt_chrm][supt_pos]=m_chrm_discord_pos[supt_chrm][supt_pos]+1
                        else:
                            if m_chrm_discord_pos_unsupt.has_key(supt_chrm)==False:
                                m_chrm_discord_pos_unsupt[supt_chrm]={}
                            if m_chrm_discord_pos_unsupt[supt_chrm].has_key(supt_pos)==False:
                                m_chrm_discord_pos_unsupt[supt_chrm][supt_pos]=1
                            else:
                                m_chrm_discord_pos_unsupt[supt_chrm][supt_pos]=m_chrm_discord_pos_unsupt[supt_chrm][supt_pos]+1

        for chrm in chrm_list:
            sf_candidate_only_clip="{0}{1}_candidate_sites.tmp".format(sf_out_folder, chrm)
            sf_final_sites="{0}{1}_sites.txt".format(sf_out_folder, chrm)
            with open(sf_candidate_only_clip) as fin_candidate, open(sf_final_sites,"w") as fout_final:
                for line in fin_candidate:
                    fields=line.split() #pos, n_clip, n_unsupport_clip, n_hard_clip, n_01pair
                    pos=int(fields[0])
                    n_clip=fields[1]
                    n_unsupt_clip=fields[2]
                    n_hard_clip=fields[3]
                    n_01pair=int(fields[4])

                    n_discordant=self.cnt_discordant_support(m_chrm_discord_pos[chrm], pos, self.dist2)
                    if n_discordant<pe_cutoff:
                        continue

                    n_unsupt_discordant=self.cnt_discordant_support(m_chrm_discord_pos_unsupt[chrm], pos, self.dist2)

                    fout_final.write(str(pos)+" "+n_clip+" "+n_unsupt_clip+" "+n_hard_clip+" "+str(n_01pair)+
                                     " "+str(n_discordant)+" "+str(n_unsupt_discordant)+"\n")

def usage():
    return

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "clfk:h:s:r:p:a:q:d:e:t:z:v:g:m:o:x:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    b_call_clip=False
    b_locate=False
    b_final=False
    k=21
    min_clip_lenth=21
    sf_high_freq_kmers=""
    sf_simple_rep="simple.ref"
    hit_ratio=0.5
    soft_clip_cutoff=3
    all_clip_cutoff=10
    pe_cutoff=3
    n_extend=15
    anchor_mapq=30
    rlenth=100 #read length
    insert_size=300
    derivation=50
    sf_gap_pos="GRCh38_gap_pos.txt"
    chrm=""
    sf_out_folder="./"
    sf_fai=""
    for o, a in opts:
        if o == "-c":
            b_call_clip=True
        elif o == "-l":
            b_locate=True
        elif o=="-f":
            b_final=True
        elif o == "-k":
            k=int(a)
            min_clip_lenth=k
        elif o == "-h":
            sf_high_freq_kmers=a
        elif o == "-s":
            sf_simple_rep=a
        elif o == "-r":
            hit_ratio=float(a)
        elif o == "-p":
            soft_clip_cutoff=int(a)
        elif o == "-a":
            all_clip_cutoff=int(a)
        elif o == "-q":
            anchor_mapq=int(a)
        elif o == "-d":
            pe_cutoff=int(a)
        elif o == "-e":
            n_extend=int(a)
        elif o == "-t":
            rlenth=int(a)
        elif o == "-z":
            insert_size=int(a)
        elif o == "-v":
            derivation=int(a)
        elif o =="-g":
            sf_gap_pos=a
        elif o == "-m":
            chrm=a
        elif o == "-o":
            sf_out_folder=a
            if sf_out_folder[-1]!="/":
                sf_out_folder=sf_out_folder+"/"
        elif o=="-x":
            sf_fai=a
        else:
            assert False, "unhandled option"
    
    rl=REPLocator(insert_size, derivation)
    sf_tmp_clip_reads="{0}{1}_clipped_reads.tmp".format(sf_out_folder, chrm)
    sf_tmp_01_pair="{0}{1}_unmap_matemap_reads.tmp".format(sf_out_folder, chrm)
    sf_tmp_discord="{0}discordant_gnrted_from_{1}.tmp".format(sf_out_folder, chrm)
    sf_need_validate="{0}{1}_pe_need_validate.tmp".format(sf_out_folder, chrm) #gnrted from Run_REPLocate_multi_thread2.py
    sf_discord_valided="{0}{1}_pe_validated.tmp".format(sf_out_folder, chrm)#PEs (of all chroms) validated from chrm
    sf_candidate_only_clip="{0}{1}_candidate_sites.tmp".format(sf_out_folder, chrm)#candidate pos only with cutoff of clipped reads
    if b_call_clip==True:
        rl.collect_candidate_clipped_01pair_reads(sf_high_freq_kmers, sf_simple_rep, anchor_mapq,
                              min_clip_lenth, k, hit_ratio, rlenth, sf_tmp_clip_reads, sf_tmp_01_pair, sf_tmp_discord)
    elif b_locate==True:
        rl.locate_rep_sites(sf_high_freq_kmers, sf_tmp_clip_reads, sf_tmp_01_pair, sf_need_validate, sf_gap_pos, chrm, anchor_mapq,
                            k, hit_ratio, n_extend, soft_clip_cutoff, all_clip_cutoff, sf_discord_valided, sf_candidate_only_clip)
    elif b_final==True:
        rl.call_final_sites(sf_fai, sf_out_folder, pe_cutoff)
