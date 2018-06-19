'''
Created on 26 Apr 2018

@author: jingwenwang
'''
from sys import argv
import pysam
import collections

input_file=pysam.AlignmentFile(argv[1],"rb")
end5=open(argv[2],"w")
end3=open(argv[3],"w")

TSS=collections.defaultdict(lambda:collections.defaultdict())
TES=collections.defaultdict(lambda:collections.defaultdict())
last_chr=None

for read in input_file.fetch(until_eof=True):
    if read.is_reverse==False:
        ch=read.reference_name
        start=read.reference_start+1
        end=read.reference_start+read.template_length
        if ch != last_chr:
            for key in sorted(TSS.keys()):
                for strand in TSS[key].keys():
                    end5.write(last_chr+"\t"+str(key)+"\t"+strand+"\t"+str(TSS[key][strand])+"\n")
            for key in sorted(TES.keys()):
                for strand in TES[key].keys():
                    end3.write(last_chr+"\t"+str(key)+"\t"+strand+"\t"+str(TES[key][strand])+"\n")
            last_chr=ch
            TSS.clear()
            TES.clear()
        
        if read.is_read1:
            try: TSS[start]["+"] += 1
            except KeyError: TSS[start]["+"] = 1
            try: TES[end]["+"] += 1 
            except KeyError: TES[end]["+"] = 1
        else:
            try: TSS[end]["-"] += 1
            except KeyError: TSS[end]["-"] = 1
            try: TES[start]["-"] += 1
            except KeyError: TES[start]["-"] = 1

        

input_file.close()
end5.close()
end3.close()
