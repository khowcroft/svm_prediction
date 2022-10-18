# usage example:
#Requires uppercase

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
import cProfile
import string
#cProfile.run('foo()')

#python table.py test_input.csv test_output.csv test_model_list.csv
#python table.py test_input.csv test_output.csv table_files.csv
os.chdir("/mnt/d/table")

def invert(str, char_dict):
    str = str.translate(char_dict)
    return str[::-1]
    #nstr = ""
    #for c in str:
        #nstr = char_dict.get(c) +nstr
    #return str

def model_score(model, seq_list, out_rows, char_dict):
    wtfile = model
    infile=open(wtfile,"r")
    wt={}
    for line in infile:
        f=line.strip().split()
        wt[f[0]]=f[1]
        wt[invert(f[0], char_dict)]=f[1]
    k=len(f[0])
    infile.close()

    for j in range(len(seq_list)):
        seq1 = seq_list[j][1]
        seq2 = seq_list[j][2]
        if len(seq1)<2*k-1 or len(seq2)<2*k-1 :
            print("seqs should be at least {0}bp to accurately score with {1}-mer weights".format(2*k-1,k))
            exit()
        

        s1=0
        s2=0
        for i in range(0,len(seq1)-k+1):
            a=seq1[i:i+k]
            s1 += float(wt[a])
            b=seq2[i:i+k]
            s2 += float(wt[b])
        out_rows[j].append( "{0:6.2f}".format(s2-s1)) #we could add more decimals


def main(argv = sys.argv):
    char_dict = {"A" : "T", "T": "A", "C": "G", "G": "C"}
    char_dict = "a".maketrans(char_dict) #why does this work (what does the string do?)

    if len(argv) != 4:
        print("Usage: python", argv[0],"<infile> <outfile> <model_lists>")
        exit(0)

    
    model_file = open(argv[3], "r")
    header = ["seqID"]
    model_list = []
    for model in model_file:
        model = model.strip()
        model_list.append(model)
        model = model.replace("/mnt/d/table/", "")
        header.append(model)
    infile=open(argv[1],"r")
    #if theres a header or something
    seq_list = []
    out_rows = []
    for line in infile:
        line =line.strip().split(",")
        seqID = line[0]
        seq1 = line[1]
        seq2 = line[2]
        seq_list.append([seqID, seq1, seq2])
        out_rows.append([seqID])

    infile.close()
    
    for model in model_list:
        # model ex: TF_E3_530_hg38_300_top10k_vs_neg1x_avg_weights.out
        model_score(model, seq_list, out_rows, char_dict)

    model_file.close()
    out_rows.insert(0,header)
    out = open(argv[2], 'w')
    csvwriter = csv.writer(out)
    csvwriter.writerows(out_rows)
    out.close        
    
#seqID,DHS_E3_100_300_noproms_nc30_hg38_top10k_vs_neg1x_avg_weights.out,TF_E3_530_hg38_300_top10k_vs_neg1x_avg_weights.out
#1,-16.05,-24.94
#2,-16.05,-24.94

    

if __name__ == '__main__': cProfile.run('main()')