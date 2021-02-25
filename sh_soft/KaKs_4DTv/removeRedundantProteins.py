#!/usr/bin/env python3
#此脚本用于从基因组的蛋白质文件中，找出每个基因最长的转录本
import sys,getopt
def usage():
    print('usage:python3 removeRedundantProteins.py -i <in_fasta> -o <out_fasta> <-h>')
    return
def removeRedundant(in_file,out_file):
    gene_dic = {}
    flag = ''
    with open (in_file) as in_fasta:
        for line in in_fasta:
            if '>' in line:
                line1 = line.strip('>\n')
                line2 = line1.split('.')
                li = line2[0]
                flag = li
                try:
                    gene_dic[li]
                except KeyError:
                    gene_dic[li] = [line]
                else:
                    gene_dic[li].append(line)
            else:
                gene_dic[flag][-1] += line
    with open (out_file,'w') as out_fasta:
        for k,v in gene_dic.items():
            if len(v) == 1:
                out_fasta.write(gene_dic[k][0])
            else:
                trans_max = ''
                for trans in gene_dic[k]:
                    a = len(list(trans))
                    b = len(list(trans_max))
                    if a > b:
                        trans_max = trans
                out_fasta.write(trans_max)
def main(argv):
    try:
        opts, args = getopt.getopt(argv,'hi:o:')
    except getopt.GetoptError:
        usage()
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            in_fasta_name = arg
        elif opt == '-o':
            outfile_name = arg
    try:
        removeRedundant(in_fasta_name,outfile_name)
    except UnboundLocalError:
        usage()
    return
if __name__ == '__main__':
    main(sys.argv[1:])
