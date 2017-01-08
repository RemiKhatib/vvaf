!There are 2 subroutines:
!read_step1: which reads only the first step (to get information from the files and not from the user)
!reading: reads the full files (vel and pos)



!======================
!Getting information:
!-natoms
!-starting step (istep)
!-filling fmass
!-filling atom_type
!======================
!The first step is readed to get this information
subroutine read_step1(natoms,istep,kindmax,fmass,atom_type,ptr_info)
  use mymod, only : info_atoms !The defined type info_atoms is re-used
  implicit none

  !From main
  integer natoms,istep,kindmax
  double precision, dimension(:), allocatable :: fmass  
  character(len=3), dimension(:), allocatable :: atom_type
  type (info_atoms), dimension(:), pointer :: ptr_info 

  !Local
  integer i,j,iatom
  double precision amass_tot
  character char

  !Initialization
  i=0 ; j=0 ; iatom=0
  amass_tot = 0.d0
  char=''

  
  !Reading the header
  read(11,*)natoms
  read(11,*)char,char,istep
  
  allocate(fmass(natoms),atom_type(natoms))
  fmass=0.d0
  !Reading the body and filling fmass
  do iatom = 1,natoms
     read(11,*)atom_type(iatom)
     do i=1,kindmax
        if(atom_type(iatom).eq.ptr_info(i)%name)then 
           amass_tot = amass_tot + ptr_info(i)%mw
           fmass(iatom) = ptr_info(i)%mw
           ptr_info(i)%amount=ptr_info(i)%amount+1
           exit
        endif

        !Procedure if an non-referenced atom is present
        if(i==kindmax)then
           write(*,*)"The atom ",atom_type(iatom)," is not defined." 
           write(*,*)"&MAIN"
           write(*,*)"  &ATOM_WEIGHT"
           write(*,*)"    ",atom_type(iatom),"  'molar weight'"
           write(*,*)"  &END ATOM_WEIGHT"
           write(*,*)"&END MAIN"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        endif
     enddo
  enddo
  fmass(:)=fmass(:)/amass_tot

  do i=1,kindmax
     if(ptr_info(i)%amount>0)write(*,*)'M(',ptr_info(i)%name,')=',ptr_info(i)%mw,' g.mol-1'
  enddo
  write(*,*)'There are:'
  do i=1,kindmax
     if(ptr_info(i)%amount>0)write(*,*)ptr_info(i)%amount,ptr_info(i)%name,&
          ' %w=',100.d0*dble(ptr_info(i)%amount)*ptr_info(i)%mw/amass_tot
  enddo
  write(*,*)'Total weight',amass_tot
  write(*,*)

  !We come back to the beginning and we jump the 2 header-lines
  rewind(11)
  read(11,*)
  read(11,*)

end subroutine read_step1


!=====================================
!Files reading and velocity adjustment
!=====================================
subroutine reading(filepos,natoms,nmovie,istep,fmass,a,b,c,step_min,step_max,&
     vx,vy,vz,x,y,z)
  implicit none

  !From main
  integer natoms,nmovie,istep,step_min,step_max
  double precision a,b,c
  double precision, dimension(:), allocatable :: fmass  
  double precision, dimension(natoms,nmovie) :: vx,vy,vz,x,y,z !Use a 3D table v(3,nmax,nmovie_max) instead of vx,vy and vz is not efficient
  character(len=100) filepos
  
  !Local
  integer i,j,iatom,imovie,end_file
  double precision, dimension(:), allocatable :: vcmx,vcmy,vcmz !Use a 2D table vcm(3,nmovie_max) instead of vcmx, vcmy, vcmz does not change the calculation time
  character char
  character(len=3) symb

  !Initialization
  i=istep  ; j=istep ; iatom=0
  imovie=1 ; end_file=0
  char=''  ; symb=''


  write(*,*)'End of INPUT reading'
  write(*,*)
  write(*,*)'Number of atoms in molecule:',natoms
  write(*,*)'Maximal number of movie steps:',nmovie!The actual value of nmovie is >= to the true value

  do while ((i<=step_max).or.(imovie<=nmovie))!If the perfect case (good values of the user and no propblem in the vel-file), these 2 conditions should be true / false at the same time.
     if(mod(imovie,1000)==0)write(*,*)'imovie',imovie

     if((i<=j).and.(imovie>1))then
        write(*,*)"!!!The step",i," is redundant. Please check you velocity file!!! VVAF is still running."!Problems already encountered in some files
     endif
     if(imovie>1)j=i


     if(i>=step_min)then
        !Data record
        do iatom = 1,natoms
           read(11,*)symb,vx(iatom,imovie),vy(iatom,imovie),vz(iatom,imovie)
           if(filepos/='')then
              read(12,*)symb,x(iatom,imovie),y(iatom,imovie),z(iatom,imovie)
           endif
        enddo
        imovie=imovie+1

     else 
        do iatom = 1,natoms
           read(11,*)
           if(filepos/='')then
              read(12,*)
           endif
        enddo
     endif

     read(11,*,iostat=end_file)
     if(end_file/=0)exit
     read(11,*)char,char,i
     if(filepos/='')then
        read(12,*)
        read(12,*)
     endif
     
  enddo

  !======================================!
  !             BE CAREFUL:              !
  !      THE VALUE USED FOR NMOVIE       !
  !IS NOT THOSE WRITTEN IN THE INPUT FILE!
  !======================================!
  nmovie=imovie-1!imovie-1 is the exact value of steps really recorded
  allocate(vcmx(nmovie),vcmy(nmovie),vcmz(nmovie))
  vcmx=0.d0 ; vcmy=0.d0 ; vcmz=0.d0
  close(11)
  if(filepos/='')close(12)
  write(*,*)nmovie,"steps were recorded and the final step is ",i!Just to give informations in the case of redundant steps.
  write(*,*)'End of TRAJECTORY reading'
  write(*,*)


  !=====================================================================
  !Remove ensemble translation from velocities of atoms 
  !Problems could come from the standard deviation which is pretty high:
  !- vx,vy,vz = 1E-4 ua
  !- drift average = 1E-6 ua
  !- drift std. dev = 1E-5 ua
  !=====================================================================
  do imovie = 1,nmovie
     vcmx(imovie) = 0.d0
     vcmy(imovie) = 0.d0
     vcmz(imovie) = 0.d0
     do iatom = 1,natoms
        vcmx(imovie) = vcmx(imovie) + fmass(iatom)*vx(iatom,imovie)
        vcmy(imovie) = vcmy(imovie) + fmass(iatom)*vy(iatom,imovie)
        vcmz(imovie) = vcmz(imovie) + fmass(iatom)*vz(iatom,imovie)
     enddo

     do iatom=1,natoms
        !Removing the drift of the cell
        vx(iatom,imovie) = vx(iatom,imovie) - vcmx(imovie)
        vy(iatom,imovie) = vy(iatom,imovie) - vcmy(imovie)
        vz(iatom,imovie) = vz(iatom,imovie) - vcmz(imovie)

        !Centering the atoms according to pbc
        if(filepos/='')then
           x(iatom,imovie) = x(iatom,imovie) - dble(floor(x(iatom,imovie)/a))*a
           y(iatom,imovie) = y(iatom,imovie) - dble(floor(y(iatom,imovie)/b))*b
           z(iatom,imovie) = z(iatom,imovie) - dble(floor(z(iatom,imovie)/c))*c
        endif
     enddo
  enddo

  deallocate(vcmx,vcmy,vcmz)

end subroutine reading
