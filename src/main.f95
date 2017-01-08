!The aim of this program is to calculate velocity-velocity auto-/hetero-correlation
!Some Extra treatment can be done like:
!-Projection of the velocity on x/y/z
!-Using a prefactor depending on the bond orientation
!-Selection of a specific region of the box
!-Selection of the stretching mode

!=============================================================================================
!Main
!=============================================================================================
program main
  use mymod
  implicit none

  integer :: nb_read,nmovie,natoms,kindmax,istep,step_min,step_max,Di,shift
  double precision :: a,b,c,box_size,center,surface
  double precision, dimension(:), allocatable :: fmass
  double precision, dimension(:,:), allocatable :: vx,vy,vz,x,y,z !Use a 3D table v(3,nmax,nmovie_max) instead of vx,vy and vz is not efficient
  character(len=100) :: chain,filetr,filepos,slab_orient
  character(len=3), dimension(:), allocatable :: atom_type
  logical file_exist
  type (info_atoms), pointer, dimension(:) :: ptr_info!A pointer is used instead of a dynamic table because this table is "reallocated" inside a subroutine

  !==============
  !Initialization
  !==============
  nb_read=0;nmovie=0;natoms=0;kindmax=0;istep=0;step_min=0;step_max=0;Di=1000;shift=0
  a=0.d0;b=0.d0;c=0.d0;box_size=0.d0;center=0.d0;surface=0.d0
  filetr='';filepos='';slab_orient=''
  ptr_info => null()


  !======================
  !Reading the input file
  !======================
  call getarg(1,chain)
  if(chain=='')then
     write(*,*)"You have to specify the name of the input file"
     write(*,*)
     write(*,*)"!!! END OF PROGRAM !!!"
     call exit()
  else
     inquire(file=chain, exist=file_exist)
     if(.not. file_exist)then 
        write(*,*)'!!!File ', trim(chain),' does not exist!!!'
        stop
     endif
     open(10,file=chain,action='read')
  endif


  !==============================
  !Information about common atoms
  !==============================
  call weight_info(ptr_info,kindmax)


  !===================================================
  !Reading of the parameters from standard input
  !Question about the calculation which has to be done
  !===================================================
  !nb_read=0: Only the MAIN block will be read
  call read_input(nb_read,filetr,filepos,nmovie,natoms,istep,kindmax,Di,shift,&
       ptr_info,atom_type,fmass,a,b,c,box_size,center,surface,slab_orient,step_min,step_max,&
       vx,vy,vz,x,y,z)


  !nb_read=1: The VVAF block will be read and tested to know if everything is correctly written
  call read_input(nb_read,filetr,filepos,nmovie,natoms,istep,kindmax,Di,shift,&
       ptr_info,atom_type,fmass,a,b,c,box_size,center,surface,slab_orient,step_min,step_max,&
       vx,vy,vz,x,y,z)
  
  !=====================================
  !Files reading and velocity adjustment
  !This is the most time consumming part
  !=====================================
  allocate(vx(natoms,nmovie),vy(natoms,nmovie),vz(natoms,nmovie))
  allocate(x(natoms,nmovie),y(natoms,nmovie),z(natoms,nmovie))!Postion array are allocated even if there is no position file
  vx=0.d0 ; vy=0.d0 ; vz=0.d0
  x=0.d0 ; y=0.d0 ;z=0.d0 !These arrays are initialized at 0 (they will keep this value if there is no position file)
  call reading(filepos,natoms,nmovie,istep,fmass,a,b,c,step_min,step_max,&
       vx(1,1),vy(1,1),vz(1,1),x(1,1),y(1,1),z(1,1))

  !========================================================================================================
  !Calculation of the VVAF
  !The input file is reded one more time and each VVAF is calculated at the end of the VVAF block
  !The main information are already read (MAIN block), the integrity of the file as been checked previously
  !So no problem should occur after the the time consumming part
  !nb_read=2
  !========================================================================================================
  call read_input(nb_read,filetr,filepos,nmovie,natoms,istep,kindmax,Di,shift,&
       ptr_info,atom_type,fmass,a,b,c,box_size,center,surface,slab_orient,step_min,step_max,&
       vx,vy,vz,x,y,z)

  close(10)
  deallocate(ptr_info,atom_type,fmass)
  deallocate(vx,vy,vz,x,y,z)

  write(*,*)"Remark about intensity"
  write(*,*)"Divide by 2 if you have 2 surfaces"
  write(*,*)"Divide by 2 if you want the integral from 0 to +inf only"

end program main



!====================================
!Information about the weight of atom
!(for mass center)
!====================================
subroutine weight_info(ptr_info,kindmax)
  use mymod, only : info_atoms !The defined type info_atoms is re-used
  implicit none
  !From the main
  integer kindmax
  type (info_atoms), dimension(:), pointer :: ptr_info 

  kindmax=11
  allocate (ptr_info(kindmax))!The table is dynamic because some unknown elements can be added during the vvaf

  ptr_info(1)%name ='H'  ; ptr_info(1)%mw =1.0d0     ; ptr_info(1) %amount=0
  ptr_info(2)%name ='Li' ; ptr_info(2)%mw =7.0d0     ; ptr_info(2) %amount=0
  ptr_info(3)%name ='C'  ; ptr_info(3)%mw =12.0d0    ; ptr_info(3) %amount=0
  ptr_info(4)%name ='N'  ; ptr_info(4)%mw =14.0d0    ; ptr_info(4) %amount=0
  ptr_info(5)%name ='O'  ; ptr_info(5)%mw =16.0d0    ; ptr_info(5) %amount=0
  ptr_info(6)%name ='F'  ; ptr_info(6)%mw =19.0d0    ; ptr_info(6) %amount=0
  ptr_info(7)%name ='Na' ; ptr_info(7)%mw =22.99d0   ; ptr_info(7) %amount=0
  ptr_info(8)%name ='Cl' ; ptr_info(8)%mw =35.453d0  ; ptr_info(8) %amount=0
  ptr_info(9)%name ='Ca' ; ptr_info(9)%mw =40.1d0    ; ptr_info(9) %amount=0
  ptr_info(10)%name='I'  ; ptr_info(10)%mw=126.9d0   ; ptr_info(10)%amount=0
  ptr_info(11)%name='Pt' ; ptr_info(11)%mw=195.084d0 ; ptr_info(11)%amount=0

end subroutine weight_info




