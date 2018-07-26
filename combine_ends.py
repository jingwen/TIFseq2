from sys import argv
import pysam
import collections

hits_tag = "NH"
all_read=0

class Read(object):
    #initiate a Read class to store reads with same ID
    def __init__(self, name):
        self.name = name
        self.end5_count = 0
        self.end3_count = 0
        self.end5_reads = collections.defaultdict(dict)
        self.end3_reads = collections.defaultdict(dict)
        self.chro=None
        self.pair_count=0
        self.pairs=[]
    #record reads at 5'ends    
    def input_5end(self, read):
        self.end5_count = read.get_tag(hits_tag)
        strand=read.is_reverse
        chro=read.reference_name
        self.fill_dict(chro,strand,self.end5_reads,read)
    #record reads at 3'ends
    def input_3end(self, read):
        self.end3_count = read.get_tag(hits_tag)
        strand=read.is_reverse
        chro=read.reference_name
        self.fill_dict(chro,strand,self.end3_reads,read)
    #update reads dict
    def fill_dict(self,chro,strand,reads,read):
        try:
            reads[chro][strand].append(read)
        except KeyError:
            reads[chro][strand]=[read]
    
    def test(self,r1,r2):
        if r1.flag<256:
            r1.flag=r1.flag+64
        else:
            r1.flag=r1.flag+64-256
        if r2.flag<256:
            r2.flag=r2.flag+128
        else:
            r2.flag=r2.flag+128-256
        self.pair_count+=1
    
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
                            if r1.reference_end<=r2.reference_end:
				self.test(r1,r2)
                            	self.pairs.append([r1,r2])
            if True in chro_end5 and False in chro_end3:
                read5=chro_end5[True]
                read3=chro_end3[False]
                for r1 in read5:
                    for r2 in read3:
                        if r1.reference_start >= r2.reference_start:
                            if r1.reference_end>=r2.reference_end:
				self.test(r1,r2)
                            	self.pairs.append([r2,r1])
    
    def print_readpairs(self,unique_file,multi_file):
        if self.pair_count==1:
            unique_file.write(self.pairs[0][0])
            unique_file.write(self.pairs[0][1])
        else:
            for read_pairs in self.pairs:
                multi_file.write(read_pairs[0])
                multi_file.write(read_pairs[1])
    
def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

unique_file_name=argv[1].replace(".name","_unique")
multi_file_name=argv[1].replace(".name","_multi")
input_file=pysam.AlignmentFile(argv[1],"rb")
unique_file=pysam.AlignmentFile(unique_file_name,"wb",template=input_file)
multi_file=pysam.AlignmentFile(multi_file_name,"wb",template=input_file)
current_read=Read(None)
for read in input_file.fetch(until_eof=True):
    umi=read.qname.split("_")[-1]
    if hamming_distance("GGGGGGGG",umi)>1:
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

fix_mate=unique_file_name.replace("unique","fixuniq")
sort_bam=fix_mate.replace("fixuniq","uniq_sorted")
pysam.fixmate(unique_file_name,fix_mate)
pysam.sort("-o",sort_bam,fix_mate)
pysam.index(sort_bam)
