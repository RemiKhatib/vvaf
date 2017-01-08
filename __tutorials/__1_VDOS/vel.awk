#The aim of this program is to produce a vel.xyz file thanks to a derivation of the pos-file
BEGIN{
    nb_line=2
    test=0
    i=1
}

{
    #Header recording
    if(NR%nb_line==1){
	line1_next=$0
	nb_line=$1+2
    }
    else if(NR%nb_line==2){
	line2_next=$0
    }




    else{
	#Body recording
	for(j=1;j<=4;j++){
	    tab1[j,i]=tab2[j,i]
	    tab2[j,i]=$j
	}
	i++

	#Printing the derivate
	if(NR%nb_line==0){
	    if(test==0){test=1}
	    else{
		print line1
		print line2
		for(i=1;i<=nb_line-2;i++){
		    print tab1[1,i]" "tab2[2,i]-tab1[2,i]" "tab2[3,i]-tab1[3,i]" "tab2[4,i]-tab1[4,i]
		}
	    }
	
	    i=1
	    line1=line1_next
	    line2=line2_next
	}

    }
}
