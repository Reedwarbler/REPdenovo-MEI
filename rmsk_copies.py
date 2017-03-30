import sys
import os
from subprocess import *
from multiprocessing import Pool
import copy_reg
import types

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

class CopyRMSK:
    def __init__(self, sf_vcf, rmsk_path, sf_copy_folder, sf_out_folder):
        self.sf_final_vcf=sf_vcf
        self.rmsk_path=rmsk_path
        self.sf_out_folder=sf_out_folder
        self.sf_copy_folder=sf_copy_folder

    def run_mask(self, sf_copy_name):
        sf_copy=self.sf_copy_folder+sf_copy_name
        if os.path.exists(sf_copy):
            cmd="{0} {1} -dir {2}".format(self.rmsk_path, sf_copy, self.sf_out_folder)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

    def read_in_list(self):
        l_copies=[]
        with open(self.sf_final_vcf) as fin_vcf:
            for line in fin_vcf:
                if line[0]=="#":
                    continue
                fields=line.split()
                chrm=fields[0]
                pos=fields[1]
                g0=fields[-2]
                g1=fields[-1]

                fname="{0}_{1}_{2}_{3}.fa".format(chrm, pos, g0, g1)
                l_copies.append(fname)
        return l_copies


    def rmsk_copies(self, n_threads):
        l_copies=self.read_in_list()
        bunch_size=1 #every time just run one command
        pool = Pool(n_threads)
        pool.map(self.run_mask, l_copies, bunch_size) #[{}]
        pool.close()
        pool.join()


    def parse_rmsk_out(self):
        l_copies=self.read_in_list()
        for sf_copy_name in l_copies:
            sf_mask=self.sf_out_folder+sf_copy_name+".out"
            with open(sf_mask) as fin_mask:
                m_read_family={}
                b_out=True
                for line in fin_mask:
                    fields=line.split()
                    if len(fields)<=1 or fields[0]=="SW" or fields[0]=="score":
                        continue
                    elif fields[0]=="There":
                        print sf_copy_name, "null"
                        b_out=False
                        continue
                    else:
                        read_id=fields[4]
                        qstart=int(fields[5])
                        qend=int(fields[6])
                        div=float(fields[1])
                        match_lenth=qend-qstart
                        n_mismatch=int(float(match_lenth)*div)

                        family=fields[10]

                        if m_read_family.has_key(read_id)==False:
                            m_read_family[read_id]={}
                        if m_read_family[read_id].has_key(family)==False:
                            m_read_family[read_id][family]=[]
                        m_read_family[read_id][family].append((match_lenth, n_mismatch))

                        # if m_family.has_key(family)==False:
                        #     m_family[family]=[]
                        # m_family[family].append(div)

                if b_out==True:
                    print sf_copy_name,

                    #first check whether all the reads are consistent
                    m_read_domnt={}
                    m_all_family_div={}
                    for read_id in m_read_family:
                        domnt_family=""
                        max_len=0
                        for family in m_read_family[read_id]:
                            for lenth, n_mismatch in m_read_family[read_id][family]:
                                if max_len<lenth:
                                    domnt_family=family
                                    max_len=lenth

                                if m_all_family_div.has_key(family)==False:
                                    m_all_family_div[family]=[]
                                m_all_family_div[family].append((lenth, n_mismatch))

                        m_read_domnt[domnt_family]=0


                    if len(m_read_domnt)>1: #########reads doesn't support same repeat family############################################
                        print "conflicts"
                        continue

                    for family in m_all_family_div:
                        total_lenth=0
                        total_mismatch=0
                        for lenth, n_mismatch in m_all_family_div[family]:
                            total_lenth=total_lenth+lenth
                            total_mismatch=total_mismatch+n_mismatch
                        print family, float(total_mismatch)/float(total_lenth), total_lenth,
                    print ""##

if __name__ == "__main__":
    sf_vcf=sys.argv[1]
    rmsk_path=sys.argv[2]
    sf_copy_folder=sys.argv[3]
    sf_out_folder=sys.argv[4]

    if sf_copy_folder[-1]!="/":
        sf_copy_folder=sf_copy_folder+"/"
    if sf_out_folder[-1]!="/":
        sf_out_folder=sf_out_folder+"/"

    n_threads=int(sys.argv[5])

    cm=CopyRMSK(sf_vcf, rmsk_path, sf_copy_folder, sf_out_folder)
    cm.rmsk_copies(n_threads)
    cm.parse_rmsk_out()
