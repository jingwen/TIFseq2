#!/bin/awk -f
BEGIN { FS="[<>]" }
/Project name/{
	split($2,n,"\"");
	printf("Project Name:%s\nSample Name\tBarcode Name\tLane_1\tLane_2\tLane_3\tLane_4\tBarcode_all\tLane_1\tLane_2\tLane_3\tLane_4\tTotal\n",n[2]);
}
/Sample name/{
	split($2,n,"\"");
	if(n[2]=="all") printf ("%s\t\t\t\t\t",n[2]);
	else printf("%s\t",n[2])
}
/Barcode name/{
	split($2,n,"\"");
	printf("%s\t",n[2])
}
/Lane number/{
	getline;
	printf("%s\t",$3);
	count+=$3
}
/\/Sample/{
	printf("%s\n",count);
	count=0
}
