#!/bin/awk -f
 
!/^#/ {s[NR]=$1;id5[NR]=$2;index5[NR]=$3;id3[NR]=$4;index3[NR]=$5}
/^#Investigator/ {split($0,p,":");PR=p[2]}
/^#Experiment/ {split($0,p,":");PJ=p[2]}
/^#Date/ {split($0,d,":");DT=d[2]}
/^#Read1Length/ {split($0,l,":");RL1=l[2]}
/^#Read2Length/ {split($0,l,":");RL2=l[2]}
END{
c=length(s)
index5[c+2]="GGGGGG";index3[c+2]="GGGGGG";id5[c+2]="unknown";id3[c+2]="unknown";s[c+2]="empty";
printf("[Header]\nInvestigator Name,%s\nExperiment Name,%s\nDate,%s\nWorkflow,GenerateFASTQ\nApplication,NextSeq FASTQ Only\n\n[Reads]\n%i\n%i\n",PR,PJ,DT,RL1,RL2)
printf("\n[Settings]\n\n")
printf("\n[Data]\nSample_ID,Sample_Name,I7_Index_ID,index,I5_Index_ID,index2\n")
for (i=2;i<c+2;i++){
	printf("%s_5_%s_3,%s_5_%s_3,%s,%s,%s,%s\n",s[i],s[i],s[i],s[i],id3[i],index3[i],id5[i],index5[i])
	printf("%s_3_%s_5,%s_3_%s_5,%s,%s,%s,%s\n",s[i],s[i],s[i],s[i],id5[i],index5[i],id3[i],index3[i])
#	printf("%s_5_empty_3,%s_5_empty_3,empty,%s,%s,%s\n",s[i],s[i],index3[c+2],id5[i],index5[i])
#	printf("empty_5_%s_3,empty_5_%s_3,%s,%s,empty,%s\n",s[i],s[i],id3[i],index3[i],index5[c+2])
#	printf("%s_3_empty_5,%s_3_empty_5,empty,%s,%s,%s\n",s[i],s[i],index5[c+2],id3[i],index3[i])
#	printf("empty_3_%s_5,empty_3_%s_5,%s,%s,empty,%s\n",s[i],s[i],id5[i],index5[i],index3[c+2])
	}
}

