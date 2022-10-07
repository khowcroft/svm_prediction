# usage example:

# python table.py "<infile> <outfile> <model_lists>"
#Input format:
#SeqID	Seq1	seq2
#rsXXX	CCCGTTTCCATGGCAACCAGA	CCCGTTTCCACGGCAACCAGA

#Script: 
#Make predictions using all the models
#Output format:
#seqID Model1_pred Model2_pred
#rsXXX -24 -16 ...

import os
import sys
import csv

os.chdir("/mnt/d/table")
def invert_char(str):
    str=str.upper()
    if str=="A":
        return "T"
    elif str=="T":
        return "A"
    elif str=="C":
        return "G"
    elif str=="G":
        return "C"
    else:
        return str
    

def invert(str):
    nstr=""
    for c in str:
        nstr=invert_char(c)+nstr
    return nstr

def model_score(model, seq1, seq2):
    wtfile = model
    line=seq1+"\t"+seq2
    infile=open(wtfile,"r")
    wt={}
    for line in infile:
        f=line.strip().split()
        k=len(f[0])
        wt[f[0]]=f[1]
        wt[invert(f[0])]=f[1]
    infile.close()

    if len(seq1)<2*k-1 or len(seq2)<2*k-1 :
        print("seqs should be at least {0}bp to accurately score with {1}-mer weights".format(2*k-1,k))
        exit()

    s1=0
    s2=0
    for i in range(0,len(seq1)-k+1):
        a=seq1[i:i+k]
        s1=s1+float(wt[a])
        b=seq2[i:i+k]
        s2=s2+float(wt[b])
    return "{0:6.2f}".format(s2-s1) #we could add more decimals


def main(argv = sys.argv):
    if len(argv) != 4:
        print("Usage: python", argv[0],"<infile> <outfile> <model_lists>")
        exit(0)

    out_rows = []
    model_file = open(argv[3], "r")
    header = ["seqID"]
    model_list = []
    for model in model_file:
        model = model.strip()
        model_list.append(model)
        model = model.replace("/mnt/d/table/", "")
        header.append(model)
    out_rows.append(header)
    infile=open(argv[1],"r")
    #if theres a header or something
    for line in infile:
        line =line.strip().split(",")
        seqID = line[0]
        seq1 = line[1]
        seq2 = line[2]

        row = [seqID]
        for model in model_list:
            # model ex: TF_E3_530_hg38_300_top10k_vs_neg1x_avg_weights.out
            score = model_score(model, seq1, seq2)
            row.append(score)
        out_rows.append(row)

    out = open(argv[2], 'w')
    csvwriter = csv.writer(out)
    csvwriter.writerows(out_rows)
    out.close        

    


    

if __name__ == '__main__': main()