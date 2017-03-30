import sys
import os
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from subprocess import *
from multiprocessing import Pool

CLIP_SPT_GNTP_0=0
CLIP_SPT_GNTP_1=1
CLIP_SPT_GNTP_2=2

#read in the repeatmasker output, and get the copy positions
def pick_copies_from_rmsk(sf_rmsk, len_cutoff):
    chrom_pos={}
    with open(sf_rmsk) as fin_rmsk:
        for line in fin_rmsk:
            fields=line.split()
            chrm=fields[4]
            chrm_start=int(fields[5])
            chrm_end=int(fields[6])

            if chrm_end-chrm_start < len_cutoff:
                continue

            if chrom_pos.has_key(chrm)==False:
                chrom_pos[chrm]=[]
            chrom_pos[chrm].append((chrm_start,chrm_end))
    return chrom_pos

def load_pos_hap_from_file(sf_cstt):
    chrom_pos={}
    with open(sf_cstt) as fin_rmsk:
        for line in fin_rmsk:
            fields=line.split()
            chrm=fields[0]
            chrm_start=int(fields[1])
            chrm_end=int(fields[2])

            if chrom_pos.has_key(chrm)==False:
                chrom_pos[chrm]=[]
            chrom_pos[chrm].append((chrm_start,chrm_end))
    return chrom_pos


def construct_ref(sf_ref, m_pos, sf_new_ref, sf_new_pos):
    l_new_seqs=[]
    m_new_pos={}
    for record in SeqIO.parse(sf_ref, "fasta"):
        sid=str(record.id)
        if m_pos.has_key(sid)==False:
            continue

        pntr=0
        seq=""
        if m_new_pos.has_key(sid)==False:
            m_new_pos[sid]=[]
        accumulation=0
        #itemp=-1
        for start,end in m_pos[sid]:
            # itemp=itemp+1
            # if mhap[sid][itemp]==0:##if 1 then remove the copy, otherwise keep no change
            #     continue

            seq=seq+str(record.seq[pntr:start])
            pntr=end

            m_new_pos[sid].append(start-accumulation)
            accumulation=accumulation+(end-start)

        seq=seq+str(record.seq[pntr:])
        new_record = SeqRecord(Seq(seq,generic_dna), sid, '', '')
        l_new_seqs.append(new_record)

    output_handle = open(sf_new_ref, "w")
    SeqIO.write(l_new_seqs, output_handle, "fasta")
    output_handle.close()

    with open(sf_new_pos,"w") as fout_new_pos:
        for chrm in m_new_pos:
            for pos in m_new_pos[chrm]:
                fout_new_pos.write(chrm+" "+str(pos)+"\n")


#randomly generate the genotype of each copy, according to the 1/2 percentage
def random_genotype(percent):
    r=random.randint(1, 100)
    if r <= percent:
        return 1
    else:
        return 2

#generate the two haplotypes
def gnrt_genotype(m_chrom_pos, percent, sf_gntp):
    mhap1={}
    mhap2={}

    with open(sf_gntp,"w") as fout_gntp:
        for chrm in m_chrom_pos:
            if mhap1.has_key(chrm)==False:
                mhap1[chrm]=[]
            if mhap2.has_key(chrm)==False:
                mhap2[chrm]=[]

            for start,end in m_chrom_pos[chrm]:
                g=random_genotype(percent)
                if g==1:
                    r=random.randint(1,2)
                    if r==1:
                        mhap1[chrm].append(1)
                        mhap2[chrm].append(0)
                        fout_gntp.write(chrm+":"+str(start)+"_"+str(end)+" 1/0\n")
                    else:
                        mhap1[chrm].append(0)
                        mhap2[chrm].append(1)
                        fout_gntp.write(chrm+":"+str(start)+"_"+str(end)+" 0/1\n")
                else:
                    mhap1[chrm].append(1)
                    mhap2[chrm].append(1)
                    fout_gntp.write(chrm+":"+str(start)+"_"+str(end)+" 1/1\n")
    return mhap1, mhap2


#load the genotype from file
def load_genotype_from_file(sf_gntp):
    mhap1={}
    mhap2={}
    with open(sf_gntp) as fin_gntp:
        for line in fin_gntp:
            fields=line.split()
            chrm_pos=fields[0]
            chrm_pos_fields=chrm_pos.split(":")
            chrm=chrm_pos_fields[0]
            # all_pos=chrm_pos_fields[1]
            # pos_fields=all_pos.split("_")

            gntp=fields[1]
            gntp_fields=gntp.split("/")
            hap1=int(gntp_fields[0])
            hap2=int(gntp_fields[1])

            if mhap1.has_key(chrm)==False:
                mhap1[chrm]=[]
            if mhap2.has_key(chrm)==False:
                mhap2[chrm]=[]

            mhap1[chrm].append(hap1)
            mhap2[chrm].append(hap2)

    return mhap1, mhap2


from itertools import izip
#generate child's genotype from parent's genotype
def gnrt_child_gntp(sf_p1, sf_p2, sf_c):
    with open(sf_p1) as fin_p1:
        with open(sf_p2) as fin_p2:
            with open(sf_c, "w") as fout_c:
                for line1, line2 in izip(fin_p1, fin_p2):
                    fields1=line1.split()
                    fields2=line2.split()
                    g1=fields1[1]
                    g2=fields2[1]

                    if fields1[0]!=fields2[0]:
                        print "ERROR: deletion positions are different for parents!!!!"
                        return

                    pos=fields1[0]
                    gntp_c=""
                    if g1=="1/1" and g2=="1/1":
                        gntp_c="1/1"
                        #fout_c.write(pos+" 1/1\n")
                    elif (g1=="1/1" and g2=="1/0") or (g1=="1/0" and g2=="1/1"):
                        r1=random.randint(1,2)
                        if r1==1:
                            gntp_c="1/1"
                        else:
                            r2=random.randint(1,2)
                            if r2==1:
                                gntp_c="0/1"
                            else:
                                gntp_c="1/0"

                    elif (g1=="1/0" and g2=="0/1") or (g1=="0/1" and g2=="1/0"):
                        r3=random.randint(1,2)
                        if r3==1:
                            r4=random.randint(1,2)
                            if r4==1:
                                gntp_c="0/0"
                            else:
                                gntp_c="1/1"
                        else:
                            r5=random.randint(1,2)
                            if r5==1:
                                gntp_c="0/1"
                            else:
                                gntp_c="1/0"
                    fout_c.write(pos+" "+gntp_c+"\n")

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

        b_hardclip=False
        cnt_m=0
        clip_flag=0
        left_clip_len=0
        right_clip_len=0

        if len(signal)==1 and signal[0]=="M":
            return 0

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
        return clip_flag

#CLIP_SPT_GNTP_0
def cnt_clipped_reads_at_sites(sf_sam, site_pos):
    cnt_clip=0
    cnt_fullmap=0
    with open(sf_sam) as fin_sam:
        for line in fin_sam:
            fields=line.split()
            pos=int(fields[3])
            if abs(pos-site_pos) <100:
                cigar=fields[5]
                if cigar=="*":
                    continue
                cutoff_len=10
                rtn=is_qualified_clipped(cigar, cutoff_len)
                if rtn>0:
                    cnt_clip=cnt_clip+1
                elif (site_pos-pos)>12:
                    cnt_fullmap=cnt_fullmap+1

    print cnt_clip, cnt_fullmap################################################################################

    cnt_all=cnt_fullmap+cnt_clip
    if cnt_all==0:
        return CLIP_SPT_GNTP_0
    ratio_fullmap=float(cnt_fullmap)/float(cnt_fullmap+cnt_clip)

    if cnt_clip<2 and cnt_all>=30:
        return CLIP_SPT_GNTP_2
    elif ratio_fullmap>0.35 and ratio_fullmap<0.65 and cnt_all>=30:
        return CLIP_SPT_GNTP_1
    else:
        return CLIP_SPT_GNTP_0


def check_one_pos(sf_bam, chrm, istart):
    sf_out="tmp/{0}_{1}.sam".format(chrm, istart)
    if os.path.exists(sf_out)==False:
        return CLIP_SPT_GNTP_0

    print sf_out ##################################################################################################

    #cmd="samtools view {0} {1}:{2}-{3} > {4}".format(sf_bam, chrm, istart-1000, istart+1000, sf_out)
    #Popen(cmd, shell = True, stdout = PIPE).communicate()
    b_left=cnt_clipped_reads_at_sites(sf_out, istart)

    #cmd="rm {0}".format(sf_out)
    #Popen(cmd, shell = True, stdout = PIPE).communicate()
    return b_left

def run_cmd(cmd):
    Popen(cmd, shell = True, stdout = PIPE).communicate()


from random import randint
def check_individual_has_rep(m_chrom_pos, sf_bam, sf_rslt):
    l_sites=[]
    for chrm in m_chrom_pos:
        for start, end in m_chrom_pos[chrm]:
            # vlu=randint(0,9)
            # if vlu%10!=0:
            #     continue
            istart=int(start)
            iend=int(end)
            l_sites.append((chrm, istart, iend))

    cmd_rmsk_list=[]

    for (chrm, istart, iend) in l_sites:
        sf_out="tmp/{0}_{1}.sam".format(chrm, istart)
        cmd="samtools view {0} {1}:{2}-{3} > {4}".format(sf_bam, chrm, istart-1000, istart+1000, sf_out)
        sf_out1="tmp/{0}_{1}.sam".format(chrm, iend)
        cmd1="samtools view {0} {1}:{2}-{3} > {4}".format(sf_bam, chrm, iend-1000, iend+1000, sf_out1)
        cmd_rmsk_list.append(cmd)
        cmd_rmsk_list.append(cmd1)

    n_threads=15
    bunch_size=1 #every time just run one command
    pool = Pool(n_threads)
    pool.map(run_cmd, cmd_rmsk_list, bunch_size)
    pool.close()
    pool.join()

    with open(sf_rslt, "w") as fout_rslt:
        for (chrm, istart, iend) in l_sites:
            b_left=check_one_pos(sf_bam, chrm, istart)
            if b_left==CLIP_SPT_GNTP_0:
                continue

            b_right=check_one_pos(sf_bam, chrm, iend)
            if b_right==b_left:
                if b_right==CLIP_SPT_GNTP_1:
                    #print chrm, start, end, 1
                    fout_rslt.write(chrm+" "+str(istart)+" "+str(iend)+" 1\n")
                else:
                    #print chrm, start, end, 2
                    fout_rslt.write(chrm+" "+str(istart)+" "+str(iend)+" 2\n")


def parse(sf_all, sf_pos):
    m_all={}
    with open(sf_all) as fin_all:
        for line in fin_all:
            fields=line.split()
            id="{0}_{1}_{2}".format(fields[4], fields[5],fields[6])
            m_all[id]=line.rstrip()
    with open(sf_pos) as fin_pos:
        for line in fin_pos:
            fields=line.split()
            id="{0}_{1}_{2}".format(fields[0], fields[1],fields[2])
            if m_all.has_key(id):
                print m_all[id]

def check_conflict(gf,gm,gc):
    if gm==2 and gf==2 and gc!=2:
        return True
    elif ((gm==2 and gf==0) or (gm==0 and gf==2) and gc!=1):
        return True
    elif ((gm==2 and gf==1) or (gm==1 and gf==2) and gc==0):
        return True
    elif ((gm==0 and gf==1) or (gm==1 and gf==0) and gc==2):
        return True
    elif gm==0 and gf==0 and gc!=0:
        return True
    return False


def check_consistent(sf_f, sf_m, sf_c):
    m_f={}
    m_m={}
    m_c={}
    with open(sf_f) as fin_f, open(sf_m) as fin_m, open(sf_c) as fin_c:
        for line in fin_f:
            fields=line.split()
            id=fields[0]+" "+fields[1]+" "+fields[2]
            m_f[id]=fields[3]

        for line in fin_m:
            fields=line.split()
            id=fields[0]+" "+fields[1]+" "+fields[2]
            m_m[id]=fields[3]

        for line in fin_c:
            fields=line.split()
            id=fields[0]+" "+fields[1]+" "+fields[2]
            m_c[id]=fields[3]

    m_conflict={}
    for sid in m_f:
        gf=int(m_f[sid])
        gm=0
        gc=0
        if m_m.has_key(sid):
            gm=int(m_m[sid])
        if m_c.has_key(sid):
            gc=int(m_c[sid])
        if check_conflict(gf,gm,gc):
            m_conflict[sid]=1
    print "# of conflict ones: ", len(m_conflict)

    for sid in m_m:
        if m_f.has_key(sid):
            continue
        gf=0
        gm=int(m_m[sid])
        gc=0
        if m_c.has_key(sid):
            gc=int(m_c[sid])

        if check_conflict(gf,gm,gc):
            m_conflict[sid]=1

    print "# of conflict ones: ", len(m_conflict)

    for sid in m_c:
        if m_f.has_key(sid):
            continue

        gf=0
        gm=0
        if m_m.has_key(sid):
            gm=int(m_m[sid])
        gc=int(m_c[sid])

        if check_conflict(gf,gm,gc):
            m_conflict[sid]=1

    print "# of conflict ones: ", len(m_conflict)

    sf_f_f="consistent_"+sf_f
    sf_m_f="consistent_"+sf_m
    sf_c_f="consistent_"+sf_c
    with open(sf_f_f,"w") as fout_f, open(sf_m_f,"w") as fout_m, open(sf_c_f,"w") as fout_c:
        for sid in m_f:
            if m_conflict.has_key(sid):
                continue
            fout_f.write(sid+" "+m_f[sid]+"\n")

        for sid in m_m:
            if m_conflict.has_key(sid):
                continue
            fout_m.write(sid+" "+m_m[sid]+"\n")

        for sid in m_c:
            if m_conflict.has_key(sid):
                continue
            fout_c.write(sid+" "+m_c[sid]+"\n")


def test_trio(sf_f, sf_m, sf_c):
    m_f={}
    m_m={}
    m_c={}
    with open(sf_f) as fin_f, open(sf_m) as fin_m, open(sf_c) as fin_c:
        for line in fin_f:
            fields=line.split()
            id=fields[0]+" "+fields[1]+" "+fields[2]
            m_f[id]=fields[3]

        for line in fin_m:
            fields=line.split()
            id=fields[0]+" "+fields[1]+" "+fields[2]
            m_m[id]=fields[3]

        for line in fin_c:
            fields=line.split()
            id=fields[0]+" "+fields[1]+" "+fields[2]
            m_c[id]=fields[3]

    for sid in m_c:
        if m_f.has_key(sid)==False and m_m.has_key(sid)==False:
            print sid

if __name__ == "__main__":
    len_cutoff=50
    #len_cutoff=400
    percent=25
    #sf_gntp="copy_genotype.txt"

    # sf_bam=sys.argv[1]
    # sf_rmsk=sys.argv[2]
    # sf_out=sys.argv[3]
    # m_chrom_pos=pick_copies_from_rmsk(sf_rmsk, len_cutoff)
    # check_individual_has_rep(m_chrom_pos, sf_bam, sf_out)

    #check_consistent(sys.argv[1], sys.argv[2], sys.argv[3])
    #test_trio(sys.argv[1], sys.argv[2], sys.argv[3])

    sf_ref=sys.argv[1]
    sf_consistent=sys.argv[2]
    sf_new_ref=sys.argv[3]
    sf_new_pos=sys.argv[4]
    m_chrom_pos=load_pos_hap_from_file(sf_consistent)
    construct_ref(sf_ref, m_chrom_pos, sf_new_ref, sf_new_pos)
