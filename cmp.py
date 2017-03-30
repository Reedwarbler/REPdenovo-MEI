import sys
import os

def read_in_rmsk(sf_rmsk):
    m_mask={}
    with open(sf_rmsk) as fin_rmsk:
        for line in fin_rmsk:
            fields=line.split()
            if len(fields)>2 and fields[0].isdigit()==True:
                qid=fields[4]
                start=int(fields[5])
                end=int(fields[6])
                type=fields[10]

                qfields=qid.split("_")
                pos_id=qfields[-1]

                if m_mask.has_key(pos_id)==False:
                    m_mask[pos_id]={}
                if m_mask[pos_id].has_key(type)==False:
                    m_mask[pos_id][type]=0
                m_mask[pos_id][type]=m_mask[pos_id][type]+(end-start)
    return m_mask


def check_diff(sf_repdenovo_rmsk, sf_repdenovo_pos, sf_melt_rmsk, sf_melt_pos):
    m_rep=read_in_rmsk(sf_repdenovo_rmsk)
    m_melt=read_in_rmsk(sf_melt_rmsk)

    l_rep_id_real_pos=[]
    m_rep_pos={}
    with open(sf_repdenovo_pos) as fin_rep:
        for line in fin_rep:
            fields=line.split()
            pos=int(fields[0])
            l_rep_id_real_pos.append(pos)
            m_rep_pos[pos]=1

    l_melt_id_real_pos=[]
    m_melt_pos={}
    with open(sf_melt_pos) as fin_melt:
        for line in fin_melt:
            fields=line.split()
            pos=int(fields[0])
            l_melt_id_real_pos.append(pos)
            m_melt_pos[pos]=1

    for index in m_rep:
        rep_pos=l_rep_id_real_pos[int(index)]
        b_hit=False
        for i in range(rep_pos-100,rep_pos+100):
            if m_melt_pos.has_key(i):
                b_hit=True
                print "1\t", index+"\t", m_rep[index]
                break

        if b_hit==False:
            print "0\t", index+"\t", m_rep[index]


def pick_rmsk(sf_index, sf_rmsk):
    m_index={}
    with open(sf_index) as fin_index:
        for line in fin_index:
            field=line.split()
            if field[0]=="0":
                m_index[field[1]]=1
    with open("0_rmst.txt","w") as fout_0:
        with open("1_rmsk.txt","w") as fout_1:
            with open(sf_rmsk) as fin_rmsk:
                for line in fin_rmsk:
                    fields=line.split()
                    if len(fields)>2 and fields[0].isdigit()==True:
                        qid=fields[4]
                        qfields=qid.split("_")
                        pos_id=qfields[-1]
                        if m_index.has_key(pos_id):
                            fout_0.write(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"
                                         +fields[5]+"\t"+fields[6]+"\t"+fields[10]+"\n")
                        else:
                            fout_1.write(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"
                                         +fields[5]+"\t"+fields[6]+"\t"+fields[10]+"\n")

def check_overlap_by_chrom(sf_index):
    cnt_overlap=0
    with open(sf_index) as fin_index:
        for line in fin_index:
            fields=line.split()
            chrm=fields[0]

            m_rslt1={}
            sf_rslt2="/data2/chongchu/Call_repeats_from_short_reads/run_melt/NA19239_high_cov/{0}_sites_gntp.txt".format(chrm)
            sf_rslt1="/data2/chongchu/Call_repeats_from_short_reads/run_repdenovo-L/NA19239_high_cov/{0}_sites_gntp.txt".format(chrm)
            with open(sf_rslt1) as fin_rslt1, open(sf_rslt2) as fin_rslt2:
                for line in fin_rslt1:
                    fields=line.split()
                    m_rslt1[int(fields[0])]=1
                for line in fin_rslt2:
                    field2=line.split()
                    pos=int(field2[0])
                    for i in range(-20,21):
                        if m_rslt1.has_key(i+pos):
                            cnt_overlap=cnt_overlap+1
    print cnt_overlap

def check_overlap_by_all():
    sf_final1="/data2/chongchu/Call_repeats_from_short_reads/run_melt/NA19239_high_cov/final_validated.vcf"
    sf_final2="/data2/chongchu/Call_repeats_from_short_reads/run_repdenovo-L/NA19239_high_cov/final_validated.vcf"
    cnt_overlap=0
    mpos1={}
    mpos2={}
    with open(sf_final1) as fin_1, open(sf_final2) as fin_2:
        for line in fin_1:
            fields=line.split()
            chrm=fields[0]
            pos=int(fields[1])
            if mpos1.has_key(chrm)==False:
                mpos1[chrm]={}
            mpos1[chrm][pos]=1
        for line in fin_2:
            fields2=line.split()
            chrm2=fields2[0]
            pos2=int(fields2[1])

            if mpos1.has_key(chrm2)==False:
                continue
            for i in range(-20,21):
                if mpos1[chrm2].has_key(i+pos2):
                    cnt_overlap=cnt_overlap+1
    print cnt_overlap


if __name__ == "__main__":
    # sf_repdenovo_rmsk=sys.argv[1]
    # sf_repdenovo_pos=sys.argv[2]
    # sf_melt_rmsk=sys.argv[3]
    # sf_melt_pos=sys.argv[4]

    # sf_repdenovo_rmsk=sys.argv[3]
    # sf_repdenovo_pos=sys.argv[4]
    # sf_melt_rmsk=sys.argv[1]
    # sf_melt_pos=sys.argv[2]
    # check_diff(sf_repdenovo_rmsk, sf_repdenovo_pos, sf_melt_rmsk, sf_melt_pos)

    #sf_index=sys.argv[1]
    # sf_rmsk=sys.argv[2]
    # pick_rmsk(sf_index, sf_rmsk)
    #check_overlap(sf_index)
    check_overlap_by_all()