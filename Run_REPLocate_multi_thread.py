import os
import sys
from subprocess import *
from multiprocessing import Pool

def run_cmd(cmd):
    Popen(cmd, shell = True, stdout = PIPE).communicate()

def locate_sites_mt(sf_fai, sf_sam, sf_kmer_lib, sf_sim, cover_ratio, nclip, npe, n_threads,sf_out_folder):
    cmd_list=[]
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            chrm=fields[0].rstrip()
            if sf_out_folder[-1]!="/":
                sf_out_folder=sf_out_folder+"/"
            cmd="samtools view {0} \"{1}\" | python REPLocate.py {2} {3} {4} {5} {6} - > {7}{8}_sites.txt"\
                .format(sf_sam, chrm, sf_kmer_lib, sf_sim, cover_ratio, nclip, npe, sf_out_folder, chrm)
            print cmd #############################################################################################
            cmd_list.append(cmd)

    bunch_size=1 #every time just run one command
    pool = Pool(n_threads)
    pool.map(run_cmd, cmd_list, bunch_size)
    pool.close()
    pool.join()

def merge_output(sf_fai, sf_out_folder):
    if sf_out_folder[-1]!="/":
        sf_out_folder=sf_out_folder+"/"
    sf_out_merged=sf_out_folder+"All_sites.txt"

    with open(sf_out_merged,"w") as fout_merged:
        with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                chrm=fields[0].rstrip()
                sf_sites="{0}{1}_sites.txt".format(sf_out_folder,chrm)
                with open(sf_sites) as fin_sites:
                    for record in fin_sites:
                        fout_merged.write(chrm+" "+record)

if __name__ == "__main__":
    sf_fai=sys.argv[1]
    sf_sam=sys.argv[2]
    sf_kmer_lib=sys.argv[3]
    sf_sim=sys.argv[4]
    cover_ratio=float(sys.argv[5])
    nclip=int(sys.argv[6])
    npe=int(sys.argv[7])
    n_threads=int(sys.argv[8])
    sf_out_folder=sys.argv[9]
    locate_sites_mt(sf_fai, sf_sam, sf_kmer_lib, sf_sim, cover_ratio, nclip, npe, n_threads,sf_out_folder)
    merge_output(sf_fai, sf_out_folder)
