#!/bin/awk -f
BEGIN { FS="[<>]" }
/Project name/{
	split($2,n,"\"");
	printf("Project:%s\nSample,Barcode,Lane_1,Lane_2,Lane_3,Lane_4,Barcode_all,Lane_1,Lane_2,Lane_3,Lane_4,Total\n",n[2]);
}
/Sample name/{
	split($2,n,"\"");
	if(n[2]=="all") printf ("%s,NA,NA,NA,NA,NA,",n[2])
	else printf("%s,",n[2])
}
/Barcode name/{
	split($2,n,"\"");
	printf("%s,",n[2])
	if(n[2]=="all") count=0
}
/Lane number/{
	getline;
	printf("%s,",$3);
	count+=$3
}
/\/Sample/{
	printf("%s\n",count);
	count=0
}
