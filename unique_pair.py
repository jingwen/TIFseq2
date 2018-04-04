from sys import argv
import pysam

hits_tag = "NH"
all_read=0

class Read(object):
    def __init__(self, name):
        self.name = name
        self.end5_count = 0
        self.end3_count = 0
        self.end5_reads = dict()
        self.end3_reads = dict()
        self.chro=None
        self.pair_count=0
        self.pairs=[]
        
    def input_5end(self, read):
        self.end5_count = read.get_tag(hits_tag)
        strand=read.is_reverse
        chro=read.reference_name
        self.fill_dict(chro,strand,self.end5_reads,read)
    
    def input_3end(self, read):
        self.end3_count = read.get_tag(hits_tag)
        strand=read.is_reverse
        chro=read.reference_name
        self.fill_dict(chro,strand,self.end3_reads,read)
    
    def fill_dict(self,chro,strand,reads,read):
        if chro in reads:
            try:
                reads[chro][strand].append(read)
            except KeyError:
                reads[chro][strand]=[read]
        else:
            reads[chro]={strand:[read]}
    
    def classify(self):
        intersect=[value for value in self.end5_reads if value in self.end3_reads]
        if len(intersect)==1:
            self.chro=intersect[0]
            chro_end5=self.end5_reads[self.chro]
            chro_end3=self.end3_reads[self.chro]
            if False in chro_end5 and True in chro_end3:
                read5=chro_end5[False]
                read3=chro_end3[True]
                for r1 in read5:
                    for r2 in read3:
                        if r1.reference_start <= r2.reference_start:
                            if r1.flag<256:
                                r1.flag=r1.flag+64
                            else:
                                r1.flag=r1.flag+64-256
                            if r2.flag<256:
                                r2.flag=r2.flag+128
                            else:
                                r2.flag=r2.flag+128-256
                            self.pair_count+=1
                            self.pairs.append([r1,r2])
            if True in chro_end5 and False in chro_end3:
                read5=chro_end5[True]
                read3=chro_end3[False]
                for r1 in read5:
                    for r2 in read3:
                        if r1.reference_start >= r2.reference_start:
                            if r1.flag<256:
                                r1.flag=r1.flag+64
                            else:
                                r1.flag=r1.flag+64-256
                            if r2.flag<256:
                                r2.flag=r2.flag+128
                            else:
                                r2.flag=r2.flag+128-256
                            self.pair_count+=1
                            self.pairs.append([r2,r1])
    
    def print_readpairs(self,unique_file,multi_file):
        if self.pair_count==1:
            unique_file.write(self.pairs[0][0])
            unique_file.write(self.pairs[0][1])
        else:
            for read_pairs in self.pairs:
                multi_file.write(read_pairs[0])
                multi_file.write(read_pairs[1])

input_file=pysam.AlignmentFile(argv[1],"rb")
unique_file=pysam.AlignmentFile(argv[2],"wb",template=input_file)
multi_file=pysam.AlignmentFile(argv[3],"wb",template=input_file)
current_read=Read(None)
for read in input_file.fetch(until_eof=True):
    name = read.qname
    if name!=current_read.name:
        try: 
            current_read.classify() 
        except KeyError: 
            pass
        if current_read.pair_count>0:
            current_read.print_readpairs(unique_file,multi_file)
        current_read=Read(name)
    if "5end" in read.get_tag("RG"):
        current_read.input_5end(read)
    elif "3end" in read.get_tag("RG"):
        current_read.input_3end(read)
unique_file.close()
multi_file.close()
input_file.close()

fix_mate=argv[2].replace("_unique","_fixmate")
sort_bam=fix_mate.replace("_fixmate","_sorted")
pysam.fixmate(argv[2],fix_mate)
pysam.sort("-o",sort_bam,fix_mate)
pysam.index(sort_bam)
