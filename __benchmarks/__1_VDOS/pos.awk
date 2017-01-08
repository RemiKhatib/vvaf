#The goal of this program is to simulate normal modes
#The modes are those obtained by a vibrational analysis
#The center of mass is (0,0,0)

BEGIN{
    nb_step=20000
    pi=atan2(1,1)*4
    f1=0.07
    f2=0.02
    dt=0.5 #fs.step-1
    c=0.00003#cm.fs-1
    
    #Frequencies in rad.step-1
    fB =1600*2*pi*c*dt
    fSS=3200*2*pi*c*dt
    fAS=3600*2*pi*c*dt
    
    
    #Everything is in the plane
    yO=0
    y1=0
    y2=0

    #Equilibrium positions
    xOe= 0.000 ; yOe=0.000 ; zOe=-0.066
    x1e=-0.768 ; y1e=0.000 ; z1e= 0.530
    x2e= 0.768 ; y2e=0.000 ; z2e= 0.530
    
    for(i=0;i<nb_step;i++){
	wB =fB*i
	wSS=fSS*i
	wAS=fAS*i
	
	#Header
	print "18"
	print "i = "i", time = "i*dt
	
	#AS)
	xO=xOe+0.07020*f1*cos(wAS) ; zO=zOe
	x1=x1e-0.55714*f1*cos(wAS) ; z1=z1e+0.43259*f1*cos(wAS)
	x2=x2e-0.55714*f1*cos(wAS) ; z2=z1e-0.43259*f1*cos(wAS)

	printf("O %f %f %f\n",xO,yO,zO+3)
	printf("H %f %f %f\n",x1,y1,z1+3)
	printf("H %f %f %f\n",x2,y2,z2+3)
	printf("O %f %f %f\n",-xO,-yO,-zO-3)
	printf("H %f %f %f\n",-x1,-y1,-z1-3)
	printf("H %f %f %f\n",-x2,-y2,-z2-3)
	
	#SS)
	xO=xOe                     ; zO=zOe+0.04926*f1*cos(wSS)
	x1=x1e+0.58814*f1*cos(wSS) ; z1=z1e-0.39099*f1*cos(wSS)
	x2=x2e-0.58814*f1*cos(wSS) ; z2=z1e-0.39099*f1*cos(wSS)

	printf("O %f %f %f\n",xO,yO,zO+6)
	printf("H %f %f %f\n",x1,y1,z1+6)
	printf("H %f %f %f\n",x2,y2,z2+6)
	printf("O %f %f %f\n",-xO,-yO,-zO-6)
	printf("H %f %f %f\n",-x1,-y1,-z1-6)
	printf("H %f %f %f\n",-x2,-y2,-z2-6)

	#B)
	xO=xOe                    ; zO=zOe+0.07114*f2*cos(wB)
	x1=x1e-0.42267*f2*cos(wB) ; z1=z1e-0.56464*f2*cos(wB)
	x2=x2e+0.42267*f2*cos(wB) ; z2=z1e-0.56464*f2*cos(wB)

	printf("O %f %f %f\n",xO,yO,zO+9)
	printf("H %f %f %f\n",x1,y1,z1+9)
	printf("H %f %f %f\n",x2,y2,z2+9)
	printf("O %f %f %f\n",-xO,-yO,-zO-9)
	printf("H %f %f %f\n",-x1,-y1,-z1-9)
	printf("H %f %f %f\n",-x2,-y2,-z2-9)
	
    }
}
