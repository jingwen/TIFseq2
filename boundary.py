'''
Created on 26 Apr 2018

@author: jingwenwang
'''
from sys import argv
import pysam
import collections

input_file=pysam.AlignmentFile(argv[1],"rb")
end5_file=argv[1].replace(".bam","_5end.ctss")
end3_file=argv[1].replace(".bam","_3end.ctss")
boundary_file=argv[1].replace(".bam","_boundary.txt")
end5=open(end5_file,"w")
end3=open(end3_file,"w")
border=open(boundary_file,"w")
border.write("seqname\tend5\tend3\tstrand\tcount\n")

TSS=collections.defaultdict(lambda:collections.defaultdict())
TES=collections.defaultdict(lambda:collections.defaultdict())
boundary=collections.defaultdict(lambda:collections.defaultdict())
last_chr=None

for read in input_file.fetch(until_eof=True):
    if read.is_reverse==False and "_" not in read.reference_name:
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
            for key in sorted(boundary.keys()):
                for strand in boundary[key].keys():
                    border.write(last_chr+"\t"+str(key[0])+"\t"+str(key[1])+"\t"+strand+"\t"+str(boundary[key][strand])+"\n")
            last_chr=ch
            TSS.clear()
            TES.clear()
            boundary.clear()
        
        if read.is_read1:
            try: TSS[start]["+"] += 1
            except KeyError: TSS[start]["+"] = 1
            try: TES[end]["+"] += 1 
            except KeyError: TES[end]["+"] = 1
            try: boundary[(start,end)]["+"] += 1
            except KeyError: boundary[(start,end)]["+"] = 1
        else:
            try: TSS[end]["-"] += 1
            except KeyError: TSS[end]["-"] = 1
            try: TES[start]["-"] += 1
            except KeyError: TES[start]["-"] = 1
            try: boundary[(end,start)]["-"] += 1
            except KeyError: boundary[(end,start)]["-"] = 1

input_file.close()
end5.close()
end3.close()
