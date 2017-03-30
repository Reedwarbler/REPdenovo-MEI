import os
import sys
from subprocess import *
from multiprocessing import Pool

def run_cmd(cmd):
    Popen(cmd, shell = True, stdout = PIPE).communicate()

def locate_sites_mt(sf_fai, sf_sam, sf_kmer_lib, sf_sim, anchor_mapq, k, cover_ratio, rlength, n_softclip, n_allclip,
                    npe, n_extend, n_threads,sf_out_folder):
    #First round
    cmd_list=[]
    cmd2_list=[]
    chrm_list=[]
#
# samtools view {0} chr1 | python REPLocate2.py -c -m chr1 -o ./NA19239_high_cov/ -h low_cov_kmer_lib.fa
# samtools view {0} chr1 | python REPLocate2.py -l -m chr1 -o ./NA19239_high_cov/ - > ./NA19239_high_cov/chr1_sites.txt
#
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            chrm=fields[0].rstrip()
            chrm_list.append(chrm)
            if sf_out_folder[-1]!="/":
                sf_out_folder=sf_out_folder+"/"
            cmd="samtools view {0} \"{1}\" | python REPLocate2.py -c " \
                "-m {2} -h {3} -q {4} -k {5} -r {6} -t {7} -o {8} - "\
                .format(sf_sam, chrm, chrm, sf_kmer_lib, anchor_mapq, k, cover_ratio, rlength, sf_out_folder)
            print cmd #############################################################################################
            cmd_list.append(cmd)

            cmd2="samtools view {0} \"{1}\" | python REPLocate2.py -l " \
                "-m {2} -q {3} -d {4} -e {5} -p {6} -o {7} -h {8} -k {9} -r {10} -a {11} - > {12}{13}_sites.txt "\
                .format(sf_sam, chrm, chrm, anchor_mapq, npe, n_extend, n_softclip, sf_out_folder,
                        sf_kmer_lib, k, cover_ratio, n_allclip, sf_out_folder, chrm)
            print cmd2#############################################################################################
            cmd2_list.append(cmd2)

    bunch_size=1 #every time just run one command
    pool = Pool(n_threads)
    pool.map(run_cmd, cmd_list, bunch_size)
    pool.close()
    pool.join()

    ##here need to merge the read ids according to chromosomes
    merge_readid_by_chrom(chrm_list, sf_out_folder)

    #Second round
    bunch_size=1 #every time just run one command
    pool = Pool(n_threads)
    pool.map(run_cmd, cmd2_list, bunch_size)
    pool.close()
    pool.join()

    cmd="python REPLocate2.py -f -x {0} -d {1} -o {2}".format(sf_fai, npe, sf_out_folder)
    run_cmd(cmd)


def merge_readid_by_chrom(chrm_list, sf_out_folder):
    m_chrm_read={}
    for chrm in chrm_list:
        sf_tmp_discord="{0}discordant_gnrted_from_{1}.tmp".format(sf_out_folder, chrm) ##gnrted from REPLocate2.py
        if os.path.exists(sf_tmp_discord)==False:
            continue
        with open(sf_tmp_discord) as fin_discord:
            for line in fin_discord:
                fields=line.split() #read_id, bfirst, map_pos, rnext, pnext
                read_id=fields[0]
                bfirst=fields[1]
                map_pos=int(fields[2])
                rnext=fields[3]
                if m_chrm_read.has_key(rnext)==False:
                    m_chrm_read[rnext]={}
                m_chrm_read[rnext][read_id]=(bfirst, chrm, map_pos)

    for chrm in chrm_list:
        sf_need_validate="{0}{1}_pe_need_validate.tmp".format(sf_out_folder, chrm) #Will be used in REPLocate2.py
        with open(sf_need_validate, "w") as fout_valdt:
            for read_id in m_chrm_read[chrm]:
                fout_valdt.write(read_id+" "+m_chrm_read[chrm][read_id][0]+" "+m_chrm_read[chrm][read_id][1]+
                                 " "+str(m_chrm_read[chrm][read_id][2])+"\n")


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
    n_softclip=int(sys.argv[6])
    n_allclip=int(sys.argv[7])
    npe=int(sys.argv[8])
    n_threads=int(sys.argv[9])
    anchor_mapq=30
    k=21
    cover_ratio=0.5
    rlength=126
    n_extend=15
    sf_out_folder=sys.argv[10]
    if sf_out_folder[-1]!="/":
        sf_out_folder=sf_out_folder+"/"
    locate_sites_mt(sf_fai, sf_sam, sf_kmer_lib, sf_sim, anchor_mapq, k, cover_ratio, rlength,
                    n_softclip, n_allclip, npe, n_extend, n_threads, sf_out_folder)
