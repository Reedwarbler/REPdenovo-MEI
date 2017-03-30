import sys
import os
from Bio import SeqIO
from subprocess import *

def read_pos(sf_f):
    m_temp={}
    with open(sf_f) as fin_f:
        for line in fin_f:
            fields=line.split()
            if len(fields)<2:
                continue
            pos=int(fields[0])
            m_temp[pos]=1
    return m_temp

def gnrt_flank(sf_pos, sf_ref, flank_lenth, sf_out):
    chrm_seq=""
    for record in SeqIO.parse(sf_ref, "fasta"):
        if str(record.id)=="chr1":
            chrm_seq=str(record.seq)
            break

    with open(sf_out,"w") as fout_flank:
        cnt=0
        with open(sf_pos) as fin_pos:
            for line in fin_pos:
                fields=line.split()
                pos=int(fields[0])

                fout_flank.write(">"+str(cnt)+"_left\n")
                fout_flank.write(chrm_seq[pos-flank_lenth-5:pos-5]+"\n")

                fout_flank.write(">"+str(cnt)+"_right\n")
                fout_flank.write(chrm_seq[pos+5:pos+flank_lenth+5]+"\n")
                cnt=cnt+1


def check_alignment(sf_algnmt, dist):
    m_algn={}
    with open(sf_algnmt) as fin_algn:
        for line in fin_algn:
            fields=line.split()
            qname=fields[0]
            ref=fields[5] ####for mini ##map to which long read
            pos=int(fields[7]) ##map position on the long read

            qfields=qname.split("_")
            id=qfields[0] #start from 0
            dir=qfields[1]

            if m_algn.has_key(id)==False:
                m_algn[id]={}
            if m_algn[id].has_key(dir)==False:
                m_algn[id][dir]={}
            m_algn[id][dir][ref]=pos

    #f_one_flank_hit=open("one_flank_hit_id.txt","w")

    m_hit={}
    for id in m_algn: ##for each candidate site
        if len(m_algn[id])<2:
            # if len(m_algn[id]==1):
            #     f_one_flank_hit.write(str(id)+"\n")
            continue
        for rid in m_algn[id]["left"]:
            if m_algn[id]["right"].has_key(rid):
                lpos=m_algn[id]["left"][rid]
                rpos=m_algn[id]["right"][rid]

                if m_hit.has_key(id)==False:
                    m_hit[id]={}
                if abs(lpos-rpos)>=dist: ##here have potential bug that, if one flank map on more than 1 positions on read
                    m_hit[id][rid]=(lpos,rpos,1) ## 1 indicate the insertion exist
                else:
                    m_hit[id][rid]=(lpos,rpos,0) ## 0 indicate the insertion absent
    return m_hit

def gnrt_reverse_complementary(s):
    lth=len(s)
    s_rc=""
    for i in range(lth-1,-1,-1):
        if s[i]=="A" or s[i]=="a":
            s_rc=s_rc+"T"
        elif s[i]=="T" or s[i]=="t":
            s_rc=s_rc+"A"
        elif s[i]=="C" or s[i]=="c":
            s_rc=s_rc+"G"
        elif s[i]=="G" or s[i]=="g":
            s_rc=s_rc+"C"
        else:
            s_rc=s_rc+s[i]
    return s_rc

##generate the genotype
#also call out the repeat copies
def get_gntp_copy_seqs(m_hit, sf_reads, flank_lenth, sf_out_folder):
    m_read_id={}
    m_gntp={}
    for id in m_hit:
        g0=0
        g1=0
        for rid in m_hit[id]:
            if m_hit[id][rid][2]==1: # indicates the copy exists
                g1=g1+1
                if m_read_id.has_key(rid)==False:
                    m_read_id[rid]={}
                m_read_id[rid][id]=(m_hit[id][rid][0],m_hit[id][rid][1])
            else:
                g0=g0+1
        m_gntp[id]=(g0,g1)

    m_copies={}

    for record in SeqIO.parse(sf_reads, "fasta"):
        rid=str(record.id)
        if m_read_id.has_key(rid):
            for id in m_read_id[rid]:
                lstart=m_read_id[rid][id][0]
                rstart=m_read_id[rid][id][1]

                start=lstart+flank_lenth-50
                end=rstart+50
                b_rc=False
                if lstart>rstart:
                    start=rstart+flank_lenth-50
                    end=lstart+50
                    b_rc=True

                seq=str(record.seq[start:end])
                if b_rc==True:
                    seq=gnrt_reverse_complementary(seq)

                if m_copies.has_key(id)==False:
                    m_copies[id]={}
                m_copies[id][rid]=seq


    if os.path.exists(sf_out_folder)==False:
        cmd="mkdir {0}".format(sf_out_folder)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    for id in m_copies:
        g0=m_gntp[id][0]
        g1=m_gntp[id][1]
        with open(sf_out_folder+"/"+str(id)+"_"+str(g0)+"_"+str(g1)+".fa","w") as fout_copy:
            for rid in m_copies[id]:
                fout_copy.write(">"+rid+"\n")
                fout_copy.write(m_copies[id][rid]+"\n")


def align_to_long_reads(sf_flank, sf_reads_index, sf_out):
    cmd="/data2/chongchu/tools/minimap/minimap -l {0} {1} > {2}".format(sf_reads_index, sf_flank, sf_out)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


if __name__ == "__main__":
    flank_lenth=300

    sf_pos=sys.argv[1]
    sf_ref=sys.argv[2]
    sf_out_flank=sys.argv[3]
    gnrt_flank(sf_pos, sf_ref, flank_lenth, sf_out_flank)

    #sf_reads_prefix="/data2/chongchu/pacbio/human.polished"
    sf_reads_prefix="/data2/chongchu/1000G_pacbio_high_cov/NA19239_pacbio/NA19239_unique_id"
    sf_reads_index=sf_reads_prefix+".mmi"
    sf_out_algnmt=sys.argv[4]
    #align_to_long_reads(sf_out_flank, sf_reads_index, sf_out_algnmt)


    dist=int(sys.argv[5])
    m_hit=check_alignment(sf_out_algnmt, dist)

    print len(m_hit)

    sf_reads=sf_reads_prefix+".fasta"
    sf_out_folder=sys.argv[6]
    get_gntp_copy_seqs(m_hit, sf_reads, flank_lenth, sf_out_folder)
