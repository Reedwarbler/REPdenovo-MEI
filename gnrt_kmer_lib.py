import sys
import os
import getopt
from subprocess import *

class KmerLib:
    # def __init__(self):
    #     # self.kmc_path=self.check_path(kmc_path)
    #     # self.kmer_lenth=k
    #     # self.sf_reads_list=reads_file
    #     # self.out_folder=self.check_path(out_folder)
    #     # self.gnrt_res="{0}kmer.res".format(self.out_folder)

    def check_path(self, path):
        if path[-1]!="/":
            path=path+"/"
        return path

    def count_kmer(self, kmc_path, k, reads_file, out_folder):
        kmc_path=self.check_path(kmc_path)
        kmer_lenth=k
        sf_reads_list=reads_file
        out_folder=self.check_path(out_folder)
        gnrt_res="{0}kmer.res".format(out_folder)

        tmp="{0}KMC_temp/".format(out_folder)
        if os.path.exists(tmp)==False:
            cmd="mkdir {0}".format(tmp)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="{0}kmc -k{1} -m32 -cs1000000000 @{2} {3} {4}"\
            .format(kmc_path, kmer_lenth, sf_reads_list, gnrt_res, tmp)
        print cmd
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    def dump_kmer(self, kmc_path, k, out_folder, min_freq, max_freq, gnrt_res):
        sf_kmc_kmer="{0}tmp_kmc_kmer.txt".format(out_folder)
        cmd="{0}kmc_dump -ci{1} -cx{2} {3} {4}".format(kmc_path, min_freq, max_freq, gnrt_res, sf_kmc_kmer)
        print cmd
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        sf_kmer="{0}{1}mers_freq{2}_{3}.fa".format(out_folder, k, min_freq, max_freq)
        cmd="awk -f cvtKMC_2_Fa.awk {0} > {1}".format(sf_kmc_kmer, sf_kmer)
        print cmd
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        #clean temp files
        cmd="rm {0}tmp_*".format(out_folder)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    def filter_tr(self, sf_kmer, sf_masked):
        cmd="/data2/chongchu/tools/trf409.linux64 {0} 2 7 7 80 10 35 21 -m -h -l 21".format(sf_kmer)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        sf_ori_mask="{0}.2.7.7.80.10.35.21.mask".format(sf_kmer)
        cmd="mv ./*.mask {0}".format(sf_ori_mask)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        with open(sf_ori_mask) as fin_m:
            with open(sf_masked,"w") as fout_m:
                for line in fin_m:
                    if line.rstrip()=="":
                        continue
                    fout_m.write(line)

    def gnrt_kmer_lib(self, sf_high, sf_low, n_round, sf_out):
        tmp="{0}_tmp".format(sf_out)
        cmd="python kmer_cluster_multi_thread.py {0} {1} {2} {3}".format(sf_high, sf_low, n_round, tmp)
        print cmd
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="cat {0} {1} > {2}".format(sf_high, tmp, sf_out)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

def usage():
    return

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "cdmgr:o:n:x:i:f:h", ["help", "output="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    b_count=False
    b_mask=False
    b_gnrt=False
    sout=""
    reads_file=""
    sinput=""
    freq=0
    min_cutoff=0
    max_cutoff=0
    for o, a in opts:
        if o == "-c":
            b_count = True
        elif o=="-d":
            b_count=False
        elif o=="-m":
            b_mask=True
        elif o=="-g":
            b_gnrt=True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o == "-r":
            reads_file = a
        elif o=="-i":
            sinput=a
        elif o in ("-o", "--output"):
            sout = a
        elif o == "-n":
            min_cutoff = int(a)
        elif o == "-x":
            max_cutoff = int(a)
        elif o == "-f":
            freq = int(a)
        else:
            assert False, "unhandled option"

    kmc_path="/data2/chongchu/tools/kmc2.3/"
    k=21

    kl=KmerLib()

    if b_count==True:
        kl.count_kmer(kmc_path, k, reads_file, sout)
    elif b_mask==True:
        kl.filter_tr(sinput, sout)
    elif b_gnrt==True:
        kl.gnrt_kmer_lib(reads_file, sinput ,freq, sout) #sf_high, sf_low, n_round, sf_out
    else:
        kl.dump_kmer(kmc_path, k, sout, min_cutoff, max_cutoff, sinput) #here sout is output folder
