import os
import sys
import collections

#0. remove the confilicts, null, simple, Low_complexity repeats
def parse_out_MEIs(sf_sites):
    m_sites={}
    with open(sf_sites) as fin_sites:
        for line in fin_sites:
            fields=line.split()
            if fields[1]=="null" or fields[1]=="conflicts":
                continue
            n_types=(len(fields)-1)/3
            max_len=0
            max_div=0.0
            max_type=""
            for i in range(n_types):
                type=fields[i*3+1]
                divgt=float(fields[i*3+2])
                lenth=int(fields[i*3+3])

                if lenth>max_len:
                    max_len=lenth
                    max_type=type
                    max_div=divgt
            id=fields[0]
            id_fields=id.split(".")
            id_detail_fields=id_fields[0].split("_")
            chrm_pos="{0}_{1}".format(id_detail_fields[0], id_detail_fields[1])
            lr0=int(id_detail_fields[2])
            lr1=int(id_detail_fields[3])
            #if max_type=="Simple_repeat" or max_type=="Low_complexity": #Simple_repeat #Low_complexity
            if max_type!="Simple_repeat" and max_type!="Low_complexity" and max_type!="Satellite": #Simple_repeat #Low_complexity
                print line.rstrip()
                m_sites[chrm_pos]=(max_type, max_div, max_len, lr0, lr1)
    return m_sites

#1. get out the validated ones what's their features
def get_validate_sites_features(sf_fai, sf_folder, sf_sites):
    cnt_hybird=0
    m_sites=parse_out_MEIs(sf_sites)
    sf_vcf=sf_folder+"final_validated.vcf"
    m_hybird={}
    with open(sf_vcf) as fin_pacbio:
        for line in fin_pacbio:
            fields=line.split()
            chrm=fields[0]
            pos=fields[1]
            check_id="{0}_{1}".format(chrm,pos)
            m_hybird[check_id]=0

    #         if m_sites.has_key(check_id):
    #             cnt_hybird=cnt_hybird+1
    #

    cnt_validated=0
    sf_out=sf_sites+"_MEI.txt"
    with open(sf_fai) as fin_fai, open(sf_out,"w") as fout_validated:
        for line in fin_fai:
            chrm_fields=line.split()
            chrm=chrm_fields[0]

            sf_ram_sites="{0}{1}_sites.txt".format(sf_folder, chrm)
            #sf_ram_sites="{0}{1}_sites_filtered.txt".format(sf_folder, chrm)
            with open(sf_ram_sites) as fin_raw:
                for raw_line in fin_raw:
                    raw_fields=raw_line.split()
                    raw_pos=raw_fields[0]
                    check_id="{0}_{1}".format(chrm,raw_pos)
                    if m_sites.has_key(check_id):
                        cnt_validated=cnt_validated+1
                        out_content="{0} {1} {2} {3} {4} {5} {6}\n".format(chrm, raw_line.rstrip(),
                                                                            m_sites[check_id][3], m_sites[check_id][4],
                                                                            m_sites[check_id][0],m_sites[check_id][1],
                                                                            m_sites[check_id][2])
                        fout_validated.write(out_content)
                        # print chrm, raw_line.rstrip(), m_sites[check_id][3], m_sites[check_id][4], m_sites[check_id][0], \
                        #     m_sites[check_id][1], m_sites[check_id][2]
                    if m_hybird.has_key(check_id):
                        cnt_hybird=cnt_hybird+1

    print "Hybird: ", cnt_hybird
    print "Validated: ", cnt_validated

#.  filter based on unsupported reads
def filter_according_to_unsupported_reads(sf_fai, sf_folder):
    cnt_raw=0
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            chrm_fields=line.split()
            chrm=chrm_fields[0]
            sf_raw_sites="{0}{1}_sites.txt".format(sf_folder, chrm)
            sf_new_raw_sites="{0}{1}_sites_filtered.txt".format(sf_folder, chrm)

            with open(sf_raw_sites) as fin_raw, open(sf_new_raw_sites,"w") as fout_new_raw:
                m_new_raw_pos={}
                for raw_line in fin_raw:
                    raw_fields=raw_line.split()
                    raw_pos=int(raw_fields[0])
                    clip_support=int(raw_fields[1])
                    clip_unsupport=int(raw_fields[2])
                    hard_clip=int(raw_fields[3])

                    pe_support=int(raw_fields[5])
                    pe_unsupport=int(raw_fields[6])

                    #if float(clip_unsupport)/float(clip_support+hard_clip)>=1.5  or float(pe_unsupport)/float(pe_support) > 1.5 :
                    #if float(clip_unsupport)/float(clip_support)>2.0:##filter
                    if clip_unsupport/clip_support>=2 or clip_unsupport>100:##filter
                        continue
                    m_new_raw_pos[raw_pos]=(raw_line,clip_support)
                #merge nearby sites
                m_new_raw_pos_sorted=collections.OrderedDict(sorted(m_new_raw_pos.items()))
                pre_pos=-1
                pre_clip=-1
                m_save={}
                for ipos, content in m_new_raw_pos_sorted.iteritems():
                    cur_clip=content[1]
                    if (ipos-pre_pos)>150 and pre_pos!=-1:
                        if m_save.has_key(pre_pos)==False:
                            fout_new_raw.write(m_new_raw_pos_sorted[pre_pos][0])
                            cnt_raw=cnt_raw+1
                            m_save[pre_pos]=1
                        pre_pos=ipos
                        pre_clip=cur_clip

                    else:
                        if pre_clip>cur_clip:
                            if m_save.has_key(pre_pos)==False:
                                fout_new_raw.write(m_new_raw_pos_sorted[pre_pos][0])
                                cnt_raw=cnt_raw+1
                                m_save[pre_pos]=1
                        else:
                            pre_pos=ipos
                            pre_clip=cur_clip
                if m_save.has_key(pre_pos)==False:
                    fout_new_raw.write(m_new_raw_pos_sorted[pre_pos][0])
                cnt_raw=cnt_raw+1
    print "Raw called out:", cnt_raw


#2. merge sites close to each other
def merge_nearby_sites(sf_sites):
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            chrm=fields[0]

            sf_round1="{0}{1}_sites.txt".format(sf_folder,chrm)


if __name__ == "__main__":
    sf_sites=sys.argv[1]
    m=parse_out_MEIs(sf_sites)
    #print len(m)


    sf_fai=sys.argv[1]
    sf_folder=sys.argv[2]
    sf_sites=sys.argv[3]
    #parse_out_MEIs(sf_sites)
    if sf_folder[-1]!="/":
        sf_folder=sf_folder+"/"

    #filter_according_to_unsupported_reads(sf_fai, sf_folder)
    get_validate_sites_features(sf_fai, sf_folder, sf_sites)
