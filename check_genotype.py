import os
import sys

def check_trio_consistent_from_pos(sf_p1, sf_p2, sf_c):
    m_p1=read_pos(sf_p1)
    m_p2=read_pos(sf_p2)
    m_c=read_pos(sf_c)

    cnt=0
    for pos in m_c:
        for i in range(pos-10,pos+11):
            if m_p1.has_key(i) or m_p2.has_key(i):
                cnt=cnt+1
                break
    print cnt


if __name__ == "__main__":
    # sf_p1=sys.argv[1]
    # sf_p2=sys.argv[2]
    # sf_c=sys.argv[3]
    # check_trio_consistent(sf_p1, sf_p2, sf_c)

