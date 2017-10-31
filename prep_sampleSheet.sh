#!/bin/awk -f
 
!/^#/ {s[NR]=$1;id5[NR]=$2;index5[NR]=$3;id3[NR]=$4;index3[NR]=$5;SP[NR]=$6}
/^#Investigator/ {split($0,p,":");PR=p[2]}
/^#Experiment/ {split($0,p,":");PJ=p[2]}
/^#Date/ {split($0,d,":");DT=d[2]}
/^#Read1Length/ {split($0,l,":");RL1=l[2]}
/^#Read2Length/ {split($0,l,":");RL2=l[2]}
END{
c=length(s)
index5[c+2]="GGGGGG";index3[c+2]="GGGGGG";id5[c+2]="unknown";id3[c+2]="unknown";s[c+2]="empty";
printf("[Header]\nInvestigator Name,%s\nExperiment Name,%s\nDate,%s\nWorkflow,GenerateFASTQ\nApplication,NextSeq FASTQ Only\n\n[Reads]\n%i\n%i\n\n[Setting]\n",PR,PJ,DT,RL1,RL2)
printf("\n[Data]\nSample_ID,Sample_Name,I7_Index_ID,index,I5_Index_ID,index2,Sample_project\n")
for (i=2;i<=c+2;i++){
	for (j=2;j<=c+2;j++){
		if (i==c+2 && j==c+2) {break}
		else{
			if (id3[j]!="unknown" || id5[j]!="unknown"){
				printf("%s_5_%s_3,%s_5_%s_3,%s,%s,%s,%s\n",s[i],s[j],s[i],s[j],id3[j],index3[j],id5[i],index5[i])
				printf("%s_3_%s_5,%s_3_%s_5,%s,%s,%s,%s\n",s[i],s[j],s[j],s[i],id5[j],index5[j],id3[i],index3[i])
				}
			if (id3[i]!="unknown" || id5[i]!="unknown"){
				printf("%s_5_%s_5,%s_5_%s_5,%s,%s,%s,%s\n",s[i],s[j],s[i],s[j],id5[j],index5[j],id5[i],index5[i])
				printf("%s_3_%s_3,%s_3_%s_3,%s,%s,%s,%s\n",s[i],s[j],s[j],s[i],id3[j],index3[j],id3[i],index3[i])
				}			
			}
		}
	}
}

