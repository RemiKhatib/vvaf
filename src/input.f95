!================================
!List of subroutines in this file
!================================
!subroutine read_input reads the input file
!subroutine input_main reads the MAIN block
!subroutine input_vvaf reads the VVAF block
!subroutine input_vvaf_velocity reads the VVAF > VELOCITY block
!subroutine input_vvaf_velocity_atom reads the VVAF > VELOCITY >ATOM block
!subroutine input_vvaf_velocity_atom_bondw reads the VVAF > VELOCITY > ATOM > BOND_WITH block
!subroutine input_vvaf_velocity_atom_layer reads the VVAF > VELOCITY > ATOM > LAYER block
!subroutine input_vvaf_velocity_atom_annulus reads the VVAF > VELOCITY > ATOM > ANNULUS block
!subroutine wrong_keyword output when there is a wrong keyword
!subroutine read_word replaces read to read word by word without taking into account space, comments, etc...
!subroutine make_bond fills the table options%tab_at_bond
!subroutine symb2tab converts the desired symbols into a table
!subroutine recenter is just here to apply pbc on a value 



!===================================================
!Reading of the parameters from standard input
!Question about the calculation which has to be done
!===================================================
subroutine read_input(nb_read,filetr,filepos,nmovie,natoms,istep,kindmax,Di,shift,&
     ptr_info,atom_type,fmass,a,b,c,box_size,center,surface,slab_orient,step_min,step_max,&
     vx,vy,vz,x,y,z)

  use mymod, only : info_atoms,vvaf_options,input_main,input_vvaf,read_word !The defined type info_atoms and vvaf_options are re-used
  implicit none

  !From the main
  integer nb_read,nmovie,natoms,istep,kindmax,Di,shift,step_min,step_max
  double precision :: a,b,c,box_size,center,surface
  double precision, dimension(:), allocatable :: fmass
  double precision,dimension(:,:):: vx,vy,vz,x,y,z
  character slab_orient
  character(len=100):: filetr,filepos
  character(len=3), dimension(:), allocatable :: atom_type
  type (info_atoms), dimension(:), pointer :: ptr_info 

  !Local
  integer:: i,end_file
  character(len=100) :: chain
  logical next_new
  type (vvaf_options) :: options

  !Initialization
  i=0
  end_file=0
  chain=''
  next_new=.true.

  call read_word(chain,end_file,next_new)
  !Main loop which totally reads the input file
  do while (end_file<=0)!end_file=-2 if end of line, >0 if end of file or other problems, -1 if end of a "word", else 0 


     !=====================================================================================================
     !The input file is read 3 times
     !1)Only the MAIN block (can be at the end of the file) is read and we collect the main information
     !2)We control the VVAF blocks to know if everything is correct
     !3)If everything is correct, we reread the VVAF blocks and we calculate the VVAF 
     !(because of the step 2 we know we will not have any problem) 
     !=====================================================================================================
     if(next_new.and.(chain=='&MAIN'))then!If the first word of the line is associated with the MAIN block
        if(nb_read==0)then!The MAIN block is read only once
           call input_main(nb_read,nmovie,natoms,istep,kindmax,Di,shift,step_min,step_max,a,b,c,box_size,center,surface,&
                fmass,filetr,filepos,slab_orient,atom_type,ptr_info)
        else!We pass without reading the MAIN block
           do 
              call read_word(chain,end_file,next_new)
              if(chain=="&END")then
                 call read_word(chain,end_file,next_new)
                 if(chain=="MAIN")then
                    exit!End of the main block
                 endif
              endif
           end do
        endif
     elseif((nb_read>=1).and.(chain=='&VVAF'))then!If the information from the main file are correct, the VVAF block can be read

        !===========================================================
        !Initialization and allocation of the data for the next VVAF
        !===========================================================
        !These arrays are allocated at the beginning of each VVAF block
        !These arrays are deallocated at the end of each VVAF block
        options%filename=''
        options%filename_tmp=''
        options%cross=.true.
        options%inc=.false.
        options%inc_val=0.d0
        options%v0%vproj=0                  ; options%vt%vproj=0
        options%v0%scal_proj=.false.        ; options%vt%scal_proj=.false.
        options%v0%bond_diff=.false.        ; options%vt%bond_diff=.false.
        options%v0%P=0                      ; options%vt%P=0
        options%v0%Q=0                      ; options%vt%Q=0
        options%v0%R=0                      ; options%vt%R=0
        options%v0%dip_stretch(:)=0.d0      ; options%vt%dip_stretch(:)=0.d0
        options%v0%dip_bend=0.d0            ; options%vt%dip_bend=0.d0
        options%v0%dip_ss=0.d0              ; options%vt%dip_ss=0.d0
        options%v0%dip_as=0.d0              ; options%vt%dip_as=0.d0
        options%v0%pol_stretch(:,:)=0.d0    ; options%vt%pol_stretch(:,:)=0.d0
        options%v0%pol_bend(:)=0.d0         ; options%vt%pol_bend(:)=0.d0
        options%v0%pol_ss(:)=0.d0           ; options%vt%pol_ss(:)=0.d0
        options%v0%pol_as=0.d0              ; options%vt%pol_as=0.d0
        options%v0%nb_bond=0                ; options%vt%nb_bond=0
        options%v0%dmin2=0.d0               ; options%vt%dmin2=0.d0
        options%v0%dmax2=0.d0               ; options%vt%dmax2=0.d0
        options%v0%vol_meth='n'             ; options%vt%vol_meth='n'
        options%v0%v0_at=.false.            ; options%vt%v0_at=.false.
        options%v0%min1_r=0.d0              ; options%vt%min1_r=0.d0
        options%v0%max1_r=0.d0              ; options%vt%max1_r=0.d0
        options%v0%min2_r=0.d0              ; options%vt%min2_r=0.d0
        options%v0%max2_r=0.d0              ; options%vt%max2_r=0.d0
        options%v0%min1=0.d0                ; options%vt%min1=0.d0
        options%v0%max1=0.d0                ; options%vt%max1=0.d0
        options%v0%min2=0.d0                ; options%vt%min2=0.d0
        options%v0%max2=0.d0                ; options%vt%max2=0.d0
        options%v0%inc_param=.false.        ; options%vt%inc_param=.false.
        options%v0%stretch_mode=.false.     ; options%vt%stretch_mode=.false.
        options%v0%c2v_mode=0               ; options%vt%c2v_mode=0
        options%v0%c2v_deriv=.false.        ; options%vt%c2v_deriv=.false.
        options%v0%molf2q(:,:,:)=0.d0       ; options%vt%molf2q(:,:,:)=0.d0

        allocate(&
             options%v0%tab_select(natoms) ,options%vt%tab_select(natoms),&
             options%v0%tab_at_bond(natoms),options%vt%tab_at_bond(natoms),&
             options%v0%tab_center(natoms) ,options%vt%tab_center(natoms),&
             )

        options%v0%tab_select=.false.       ; options%vt%tab_select=.false.
        options%v0%tab_at_bond=.false.      ; options%vt%tab_at_bond=.false.
        options%v0%tab_center=.false.       ; options%vt%tab_center=.false.

        !Then the information in the VVAF block are recorded
        call read_word(options%filename,end_file,next_new)

        call input_vvaf(nb_read,filepos,nmovie,natoms,options,atom_type,fmass,a,b,c,box_size,center,shift,Di,slab_orient,&
             surface,vx,vy,vz,x,y,z)

     elseif(chain=='')then

        !During the first lecture only the MAIN block is read and the main informations about array allocation are extracted
        !While during the second lecture (and further) all the blocks (MAIN and VVAF) are read
        !So it is neccessary to check if everything is good at the second lecture (nb_read=1) and not before
     elseif(nb_read>=1)then
        write(*,*)"Error outside MAIN or VVAF blocks"
        write(*,*)"Problematic keyword:",trim(chain)
        write(*,*)
        write(*,*)"&MAIN"
        write(*,*)"&END MAIN"
        write(*,*)"&VVAF"
        write(*,*)"&END VVAF"
        write(*,*)
        write(*,*)"!!! END OF PROGRAM !!!"
        call exit()
     endif

     call read_word(chain,end_file,next_new)
  enddo!one more VVAF block read

  nb_read=nb_read+1!The input file has been fully read one more time
  rewind(10)!
end subroutine read_input



!======================
!Reading the MAIN block
!Contain basic info
!About the system
!======================
subroutine input_main(nb_read,nmovie,natoms,istep,kindmax,Di,shift,step_min,step_max,a,b,c,box_size,center,surface,&
     fmass,filetr,filepos,slab_orient,atom_type,ptr_info)
  use mymod, only : info_atoms,read_step1,wrong_keyword,read_word,input_main_atweight
  implicit none

  !From the upper subroutine
  integer nb_read,nmovie,natoms,istep,kindmax,Di,shift,step_min,step_max
  double precision a,b,c,box_size,center,surface
  double precision, dimension(:), allocatable :: fmass
  character slab_orient
  character(len=100) filetr,filepos
  character(len=3), dimension(:), allocatable :: atom_type
  logical:: file_exist=.false.
  type (info_atoms), dimension(:), pointer :: ptr_info 

  !Local
  integer :: i,end_file
  character(len=100) :: chain,block_name
  logical next_new

  !Initialization
  i=0
  end_file=0
  chain=''
  block_name='MAIN'
  next_new=.false.


  do while (end_file<=0)
     call read_word(chain,end_file,next_new)

     select case(chain)
     case("VEL_FILE")!Velocity file
        call read_word(filetr,end_file,next_new)

        !If it is the first time that the input file is read, we extract some iformation
        !on the contrary, we do not do anything
        if(nb_read==0)then
           inquire(file=filetr, exist=file_exist)
           if(.not. file_exist)then 
              write(*,*)'!!!File ', trim(filetr),' does not exist!!!'
              stop
           endif
           open(11,file=filetr,action='read')
        endif

     case("CELL_PARAM" )!Cell parameters in Ang
        call read_word(chain,end_file,next_new)
        read(chain,*)a
        call read_word(chain,end_file,next_new)
        read(chain,*)b
        call read_word(chain,end_file,next_new)
        read(chain,*)c

     case("NORMAL_AXIS")!x,y,z: axis perpendicular to the surface
        call read_word(chain,end_file,next_new)
        if((chain=='X').or.(chain=='Y').or.(chain=='Z'))slab_orient=chain(1:1)!Control of the value inside slab_orient at the end of the MAIN block

     case("FIRST_STEP")!First step of the dynamic to take into account
        call read_word(chain,end_file,next_new)
        if(nb_read==0)read(chain,*)step_min!The value of step_min may have changed during the first lecture

     case("FINAL_STEP")!Last step of the dynamic to take into account 
        call read_word(chain,end_file,next_new)
        read(chain,*)step_max

     case("SAMPLE_WINDOW")!The VVAF will be defined with t in [t0-SAMPLE_WINDOW;t0+SAMPLE_WINDOW]
        call read_word(chain,end_file,next_new)
        read(chain,*)Di

     case("SAMPLE_SHIFT")!Shift between 2 t0
        call read_word(chain,end_file,next_new)
        read(chain,*)shift

     case("POS_FILE")!Name of the position file
        call read_word(filepos,end_file,next_new)
        !If it is the first time that the input file is read, we pass the first two lines
        !on the contrary, we do not do anything
        if(nb_read==0)then
           inquire(file=filepos, exist=file_exist)
           if(.not. file_exist)then 
              write(*,*)'!!!File ', trim(filepos),' does not exist!!!'
              stop
           endif
           open(12,file=filepos,action='read')
           read(12,*)!First header line
           read(12,*)!Second header line
        endif

     case("SLAB_CENTER")!Slab center
        call read_word(chain,end_file,next_new)
        read(chain,*)center
        !The centering is done at the end

     case("&ATOM_WEIGHT")!Atoms without specific molecular weight
        call input_main_atweight(kindmax,ptr_info)

     case('&END')
        call read_word(chain,end_file,next_new)!&END as to be followed by the name of the block

        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        elseif(nb_read==0)then!We control the MAIN block only once

           !Checking if mandatory keywords are used
           if(filetr=='')then
              write(*,*)"Velocity file not defined"
              write(*,*)"&MAIN"
              write(*,*)"  VEL_FILE 'filename'"
              write(*,*)"&END MAIN"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           endif

           if(&
                (abs(a)<=epsilon(0d0)).or.&
                (abs(b)<=epsilon(0d0)).or.&
                (abs(c)<=epsilon(0d0)))then
              write(*,*)"Cell parameters not defined"
              write(*,*)"&MAIN"
              write(*,*)"  CELL_PARAM 'a' 'b' 'c'"
              write(*,*)"&END MAIN"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           elseif(slab_orient=='X')then
              box_size=a
              surface=b*c
              center=center-dble(floor(center/a))*a
           elseif(slab_orient=='Y')then
              box_size=b
              surface=a*c
              center=center-dble(floor(center/b))*b
           elseif(slab_orient=='Z')then
              box_size=c
              surface=a*b
              center=center-dble(floor(center/c))*c
           else
              write(*,*)"Normal axis not defined"
              write(*,*)"&MAIN"
              write(*,*)"  NORMAL_AXIS 'X'/'Y'/'Z'"
              write(*,*)"&END MAIN"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           end if

           !Sample shift
           if(shift==0)then
              shift=2*Di
           elseif(shift<0)then
              write(*,*)"Negative value for the shift of t0"
              write(*,*)"&MAIN"
              write(*,*)"  SAMPLE_SHIFT 'number' (2*SAMPLE_WINDOW by default)"
              write(*,*)"&END MAIN"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           endif


           !======================
           !Getting information:
           !-natoms
           !-starting step (istep)
           !-filling fmass
           !-filling atom_type
           !======================
           call read_step1(natoms,istep,kindmax,fmass,atom_type,ptr_info)

           !=================================================================!
           !                          BE CAREFUL:                            !
           !The value of step_max right here is just a preliminary estimation!
           !This value is update when the velocity file is read              !
           !The final value is the exact number of steps really recorded     !
           !Therefore here, nmovie is calculated only when we read the main  !
           !block for the first time (nb_read=0)                             !
           !=================================================================!
           if((step_max>=step_min).and.(step_max>=istep))then
              if(istep>step_min)then!If the first step of the dynamic is higher than the first step defined by the user
                 nmovie=step_max-istep+1
                 step_min=istep
              else
                 nmovie=step_max-step_min+1
              endif

              !If the amount of dat to store is big
              if(nmovie*natoms>48000000)write(*,*)'Be carefull, a dynamic table with a size over 3GB will be created!'
           else
              nmovie=floor(16000000.d0/dble(natoms))
              step_max=step_min+nmovie-1
              !nmovie=floor(85000000./natoms) !A dynamic table of 5.7GB (approximatively 180,000 steps with 500 atoms) can be created
              write(*,*)"In MAIN block"
              write(*,*)'FINAL_STEP value not precised (or lower than FIRST STEP).'
              write(*,*)'So, FINAL_STEP value fixed to ',step_max,'.'
              write(*,*)'A dynamic table of 1GB will be created.'
           endif
        endif
        exit

     case('')

     case default
        call wrong_keyword(chain,block_name)

     end select
  enddo

end subroutine input_main


!==============================================================
!Subroutine which changes / adds the atomic weigth of the atoms
!==============================================================
subroutine input_main_atweight(kindmax,ptr_info)
  use mymod, only : info_atoms,read_word,wrong_keyword
  implicit none

  !From upper subroutine
  integer::kindmax
  type (info_atoms),dimension(:), pointer :: ptr_info

  !Local variables
  integer i,j,end_file
  logical next_new
  character(len=100) chain,block_name
  type (info_atoms), dimension(:), allocatable :: tab_tmp

  !Initialization
  i=0 ; j=0 ; end_file=0
  next_new=.false.
  chain='' ; block_name='ATOM_WEIGHT'

  do while (end_file<=0)
     call read_word(chain,end_file,next_new)

     select case(chain)
     case('&END')
        call read_word(chain,end_file,next_new)
        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        else
           exit
        endif

     case default!It is expected to have a Symbol and a molar weight
        do i=1,kindmax+1

           !Add a new value to the array ptr_info if the atom is unkown
           if(i==kindmax+1)then
              allocate (tab_tmp(kindmax))
              tab_tmp=ptr_info
              deallocate (ptr_info)
              kindmax=kindmax+1
              allocate (ptr_info(kindmax))
              do j=1,kindmax-1
                 ptr_info(j)=tab_tmp(j)
              enddo
              deallocate(tab_tmp)

              !Initialization of the new type
              ptr_info(kindmax)%name=chain(1:3)
              ptr_info(kindmax)%amount=0
              call read_word(chain,end_file,next_new)
              read(chain,*)ptr_info(kindmax)%mw
              call read_word(chain,end_file,next_new)
              exit

           !Change the values if atomic weight already known
           elseif(chain==ptr_info(i)%name)then
              call read_word(chain,end_file,next_new)
              read(chain,*)ptr_info(i)%mw
              call read_word(chain,end_file,next_new)
              exit

           endif
           
        enddo

     end select
  end do

end subroutine input_main_atweight



!=====================================
!Subroutine which reads the block VVAF
!This block contains the information
!About the correlation to do between
!The two velocities
!=====================================
subroutine input_vvaf(nb_read,filepos,nmovie,natoms,options,atom_type,fmass,a,b,c,box_size,center,shift,Di,slab_orient,&
     surface,vx,vy,vz,x,y,z)
  use mymod, only : vvaf_options,vvaf,wrong_keyword,read_word,&
       input_vvaf_velocity,update_layer,update_filename,error
  implicit none

  !From upper subroutine
  integer :: nb_read,nmovie,natoms,Di,shift
  double precision a,b,c,box_size,center,surface
  double precision,dimension(:,:):: vx,vy,vz,x,y,z
  double precision, dimension(natoms) :: fmass
  character:: slab_orient
  character(len=3), dimension(natoms):: atom_type
  character(len=100) :: filepos
  type (vvaf_options) :: options  

  !Local
  integer:: i,end_file
  double precision param,diff1,diff2
  double precision, dimension(-Di:Di,9):: corr
  character(len=100):: chain,block_name
  logical next_new,end_vvaf_block
  !  double precision start, finish

  !corr(:,1) is the autocorrelation
  !corr(:,2) is the intramolecular correlation 
  !corr(:,3) is the intermolecular correlation 


  !Initialization
  i=0 ; end_file=0
  param=0.d0 ; diff1=0.d0 ; diff2=0.d0
  corr(:,:)=0.d0
  chain='' ; block_name='VVAF'
  next_new=.false. ;  end_vvaf_block=.false.


  do while (end_file<=0)
     call read_word(chain,end_file,next_new)

     select case(chain)
        !Velocity to take into account if auto-correlation
        !First velocity to take into account if cross-correlation
     case("&VELOCITY0")
        call input_vvaf_velocity(1,filepos,natoms,options%v0,options%inc,options%cross,atom_type,fmass,a,b,c)

     case("&VELOCITYT")
        call input_vvaf_velocity(2,filepos,natoms,options%vt,options%inc,options%cross,atom_type,fmass,a,b,c)

     case("LAYER_INC")
        !Value of the increment if there is a layering
        !(In absolute value to not have any mistake)
        call read_word(chain,end_file,next_new)
        read(chain,*)options%inc_val
        options%inc_val=abs(options%inc_val)

     case('&END')
        call read_word(chain,end_file,next_new)!&END as to be followed by the name of the block
        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        else

           !==========================================================================
           !Before to leave the VVAF block, some verifications have to be done
           !If autocorrelation VELOCITY0 and VELOCITYT some values have to be the same
           !This table has to be update when the derived type evolves
           !REMARK: Not all the elements of options%v0 have to be paste in options%vt
           !ARBEIT WORK TRAVAIL
           !MAYBE I SHOULD CHANGE THE WAY I DEFINE THE DERIVED TYPE
           !==========================================================================
           if(.not.options%cross)then
              options%vt%bond_diff     =options%v0%bond_diff
              options%vt%tab_select    =options%v0%tab_select
              options%vt%tab_at_bond   =options%v0%tab_at_bond
              options%vt%tab_at_bond   =options%v0%tab_at_bond
              options%vt%nb_bond       =options%v0%nb_bond
              options%vt%dmin2         =options%v0%dmin2
              options%vt%dmax2         =options%v0%dmax2
              options%vt%vol_meth      =options%v0%vol_meth
              options%vt%tab_center    =options%v0%tab_center
              options%vt%v0_at         =options%v0%v0_at
              options%vt%min1_r        =options%v0%min1_r
              options%vt%max1_r        =options%v0%max1_r
              options%vt%min2_r        =options%v0%min2_r
              options%vt%max2_r        =options%v0%max2_r
              options%vt%min1          =options%v0%min1
              options%vt%max1          =options%v0%max1
              options%vt%min2          =options%v0%min2
              options%vt%max2          =options%v0%max2
              options%vt%inc_param     =options%v0%inc_param
              options%vt%stretch_mode  =options%v0%stretch_mode
              options%vt%c2v_mode      =options%v0%c2v_mode
              options%vt%molf2q(:,:,:) =options%v0%molf2q(:,:,:)
              
           endif


           !It is necessary to have at least one bond (vel%tab_at_bond=.true.)
           !to have one of these options
           !The control is done now to not impose any order between VELOCITY0 and VELOCITYT
           !even with the AUTO flag for VELOCITYT
           if(  (options%v0%c2v_mode/=0).or.(options%v0%P/=0).or.(options%v0%Q/=0).or.(options%v0%R/=0).or.&
                (options%vt%c2v_mode/=0).or.(options%vt%P/=0).or.(options%vt%Q/=0).or.(options%vt%R/=0).or.&
                options%v0%stretch_mode.or.options%v0%bond_diff.or.options%v0%scal_proj.or.&
                options%vt%stretch_mode.or.options%vt%bond_diff.or.options%vt%scal_proj)then
              if(.not.any(options%v0%tab_at_bond))then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"You need to define some bonds. E.g."
                 write(*,*)"&VVAF"
                 write(*,*)"  &VELOCITY?"
                 write(*,*)"    &ATOM"
                 write(*,*)"      &BOND_WITH"
                 write(*,*)"        INDEX 1 2 3 / SYMBOL H O / ALL"
                 write(*,*)"      &END BOND_WITH"
                 write(*,*)"    &END ATOM"
                 write(*,*)"  &END VELOCITY?"
                 write(*,*)"&VVAF"
                 write(*,*)
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif
              
              !Control of the BOND_WITH block
              if((options%v0%nb_bond<=0).or.(options%vt%nb_bond<=0))then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"You need to define some bonds. E.g."
                 write(*,*)"&VVAF"
                 write(*,*)"  &VELOCITY?"
                 write(*,*)"    &ATOM"
                 write(*,*)"      &BOND_WITH"
                 write(*,*)"        NB_BOND 1/2/3/..."
                 write(*,*)"      &END BOND_WITH"
                 write(*,*)"    &END ATOM"
                 write(*,*)"  &END VELOCITY?"
                 write(*,*)"&VVAF"
                 write(*,*)
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif

              if(options%v0%dmin2<0)options%v0%dmin2=0.d0
              if(options%vt%dmin2<0)options%vt%dmin2=0.d0

              if((options%v0%dmax2<options%v0%dmin2).or.(options%vt%dmax2<options%vt%dmin2))then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"DMAX < DMIN"
                 write(*,*)"&VVAF"
                 write(*,*)"  &VELOCITY?"
                 write(*,*)"    &ATOM"
                 write(*,*)"      &BOND_WITH"
                 write(*,*)"        DMIN ???"
                 write(*,*)"        DMAX ???"
                 write(*,*)"      &END BOND_WITH"
                 write(*,*)"    &END ATOM"
                 write(*,*)"  &END VELOCITY?"
                 write(*,*)"&VVAF"
                 write(*,*)
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif

              if(  ((options%v0%c2v_mode>0).and.(options%v0%nb_bond/=2)).or.&
                   ((options%vt%c2v_mode>0).and.(options%vt%nb_bond/=2)))then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"Symmetric / antisymmetric stretching is implemented only with 2 bonds."
                 write(*,*)"&VVAF"
                 write(*,*)"  &VELOCITY?"
                 write(*,*)"    &ATOM"
                 write(*,*)"      SIMPLE_C2V / &HESSIAN XXX (eg. C2V_AX2) ... &END HESSIAN"
                 write(*,*)"      &BOND_WITH"
                 write(*,*)"        NB_BOND 2"
                 write(*,*)"      &END BOND_WITH"
                 write(*,*)"    &END ATOM"
                 write(*,*)"  &END VELOCITY?"
                 write(*,*)"&VVAF"
                 write(*,*)
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif
              
              
              !Verification to not do a scalar product with a scalar
              if(  (options%v0%scal_proj.and.(options%v0%vproj/=0)).or.&
                   (options%vt%scal_proj.and.(options%vt%vproj/=0)))then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"Be careful, the keywords V_BOND_REF and VEL_PROJ are enable."
                 write(*,*)"You cannot project a scalar on a vector."
                 write(*,*)
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif
              if(  (&
                ((options%v0%vproj/=0).or.options%v0%scal_proj.or.(options%v0%c2v_mode>0)).and..not.&
                ((options%vt%vproj/=0).or.options%vt%scal_proj.or.(options%vt%c2v_mode>0))&
                ).or.&
                ((options%vt%vproj/=0).or.options%vt%scal_proj.or.(options%vt%c2v_mode>0)).and..not.&
                ((options%v0%vproj/=0).or.options%v0%scal_proj.or.(options%v0%c2v_mode>0))&
                )then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"Be careful, one of the keywords V_BOND_REF / VEL_PROJ / SIMPLE_C2V / &HESSIAN"
                 write(*,*)"is enable for one velocity. And not for the other one."
                 write(*,*)"You cannot make a scalar product with a scalar and a vector."
                 write(*,*)
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif
           endif !Tests about the bonds etc...

           !Verification about the keywords associated with the incremental layering
           if((abs(options%inc_val)<=epsilon(0d0)).and.options%inc)then
              write(*,*)"Problem for the VVAF ",trim(options%filename)
              write(*,*)"You want an increment for the layering, but the value of the increment is not defined"
              write(*,*)"&VVAF"
              write(*,*)"  LAYER_INC 'value'"
              write(*,*)"&VVAF"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()

           elseif((abs(options%inc_val)>epsilon(0d0)).and..not.options%inc)then
              write(*,*)"Problem for the VVAF ",trim(options%filename)
              write(*,*)"You have specified a value for the increment, but no parameter has been specified."
              write(*,*)"&VVAF"
              write(*,*)"  &VELOCITY?"
              write(*,*)"    &ATOM"
              write(*,*)"      &LAYER"
              write(*,*)"        PARAM_INC 'M1-m2' or 'm1-M2'"
              write(*,*)"      &END LAYER"
              write(*,*)"    &END ATOM"
              write(*,*)"  &END VELOCITY?"
              write(*,*)"&VVAF"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           elseif((abs(options%inc_val)>epsilon(0d0)).and.options%inc)then
              !If there is a layering with an increment
              !It is necessary to check that max1-min1=max2-min2
              !In fact due to rounding reasons, we do not compare these numbers as
              !double precision but as real
              !It hast to be true for v0 and vt
              !Therefore the thickness of the different layer evolves concomitantly
              !If one slab has a thickness of 0 it it not taken into account

              !Control for t0
              if(options%v0%vol_meth=='X')then
                 param=a
              elseif(options%v0%vol_meth=='Y')then
                 param=b
              elseif(options%v0%vol_meth=='Z')then
                 param=c
              endif
              
              if(options%v0%max1_r>=options%v0%min1_r)then
                 diff1=options%v0%max1_r-options%v0%min1_r
              else!pbc
                 diff1=param+options%v0%max1_r-options%v0%min1_r
              endif
              if(options%v0%max1_r>=options%v0%min1_r)then
                 diff2=options%v0%max2_r-options%v0%min2_r
              else!pbc
                 diff2=param+options%v0%max2_r-options%v0%min2_r
              endif
              diff1=diff1-dble(floor(diff1/param))*param
              diff2=diff2-dble(floor(diff2/param))*param
              if(&
                   (abs(diff1-diff2)>epsilon(0.d0)).and.&!Because of rounding problems, real are used and not double precison 
                   (abs(diff1)>epsilon(0d0)).and.&
                   (abs(diff2)>epsilon(0d0)))then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"For incremental layering, the slab thickness has to be the same for all the SLABS"
                 write(*,*)""
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif

              !Control for t
              if(options%vt%vol_meth=='X')then
                 param=a
              elseif(options%vt%vol_meth=='Y')then
                 param=b
              elseif(options%vt%vol_meth=='Z')then
                 param=c
              endif

              if(options%vt%max1_r>=options%vt%min1_r)then
                 diff1=options%vt%max1_r-options%vt%min1_r
              else!pbc
                 diff1=param+options%vt%max1_r-options%vt%min1_r
              endif
              if(options%vt%max1_r>=options%vt%min1_r)then
                 diff2=options%vt%max2_r-options%vt%min2_r
              else!pbc
                 diff2=param+options%vt%max2_r-options%vt%min2_r
              endif
              diff1=diff1-dble(floor(diff1/param))*param
              diff2=diff2-dble(floor(diff2/param))*param
              if(&
                   (abs(diff1-diff2)>epsilon(0.d0)).and.&!Because of rounding problems, real are used and not double precison
                   (abs(diff1)>epsilon(0d0)).and.&
                   (abs(diff2)>epsilon(0d0)))then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"For incremental layering, the slab thickness has to be the same for all the SLABS"
                 write(*,*)""
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif

              !if SLAB1 at t0 and at t have two different definitions
              !Then we have to be sure there is the same thickness among the 4 slabs
              !Here we control only SLAB1 at t0 and t because the control with SLAB2 is already done
              !The layering a t and t0 can be done along a different axis
              !I do not know why it could be interesting, but it is possible
              !So vol_meth has to be update for t and t0
              if(options%v0%vol_meth=='X')then
                 param=a
              elseif(options%v0%vol_meth=='Y')then
                 param=b
              elseif(options%v0%vol_meth=='Z')then
                 param=c
              endif
              if(options%vt%max1_r>=options%vt%min1_r)then
                 diff1=options%v0%max1_r-options%v0%min1_r
              else!pbc
                 diff1=param+options%v0%max1_r-options%v0%min1_r
              endif
              diff1=diff1-dble(floor(diff1/param))*param
              if(options%vt%vol_meth=='X')then
                 param=a
              elseif(options%vt%vol_meth=='Y')then
                 param=b
              elseif(options%vt%vol_meth=='Z')then
                 param=c
              endif
              if(options%vt%max1_r>=options%vt%min1_r)then
                 diff2=options%vt%max2_r-options%vt%min2_r
              else!pbc
                 diff2=param+options%vt%max2_r-options%vt%min2_r
              endif
              diff2=diff2-dble(floor(diff2/param))*param
              if(&
                   (abs(diff1-diff2)>epsilon(0.d0)).and.&
                   (abs(diff1)>epsilon(0d0)).and.&
                   (abs(diff2)>epsilon(0d0)))then
                 write(*,*)"Problem for the VVAF ",trim(options%filename)
                 write(*,*)"For incremental layering, the slab thickness has to be the same for all the SLABS"
                 write(*,*)""
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              endif
           end if!If there is an incremental layering


           if((options%v0%c2v_mode>0).or.(options%vt%c2v_mode>0))then !The c2v mode is studied
              !mu_y or alpha_xy or alpha_yz are non zero
              if(  (abs(options%v0%dip_stretch(2))>error) .or. &
                   (abs(options%vt%dip_stretch(2))>error) .or. &
                   (abs(options%v0%pol_stretch(1,2))>error) .or. &
                   (abs(options%vt%pol_stretch(1,2))>error) .or. &
                   (abs(options%v0%pol_stretch(2,3))>error) .or. &
                   (abs(options%vt%pol_stretch(2,3))>error) )then
                 write(*,*)"You do want to study the C_2v modes, but you "
                 write(*,*)"do not respect the symmetry rules."
                 write(*,*)"Some y components have to be set to 0"
                 write(*,*)"    &DIPOLE ?"
                 write(*,*)"      STRETCH_y 0"
                 write(*,*)"    &END DIPOLE"
                 write(*,*)"    &POLARIZABILITY ? ?"
                 write(*,*)"      STRETCH_y 0 ?"
                 write(*,*)"      STRETCH_z ? 0 ?"
                 write(*,*)"    &END POLARIZABILITY"
                 write(*,*)""
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
              end if
           else if(options%v0%c2v_deriv .or. options%vt%c2v_deriv)then
              !If the SIMPLE_C2V / &HESSIAN key word is not used, therefore it is
              !useless to use the bending
              write(*,*)"You do not want to study the normal modes (stretching version)"
              write(*,*)"but you defined some bending contribution."
              write(*,*)"Either use SIMPLE_C2V / &HESSIAN C2V_AX2, or remove the BENDING contribution"
              write(*,*)"&VVAF"
              write(*,*)"  &VELOCITY0 BOND"
              write(*,*)"    SIMPLE_C2V / &HESSIAN C2V_AX2 ... &END HESSIAN"
              write(*,*)"    ..."
              write(*,*)"  &END VELOCITY0"
              write(*,*)"  ..."
              write(*,*)"&VVAF"
              write(*,*)""
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           end if !End of c2v mode


           
           !V_BOND_REF is there that we want to project the velocity
           !on the bond. With C2V symmetry, we want to take into
           !account the whole molecule.
           if(&
                (options%v0%scal_proj .and. (options%v0%c2v_mode>0)).or.&
                (options%vt%scal_proj .and. (options%vt%c2v_mode>0))    &
              )then
              write(*,*)"You do want to study the normal modes of a MOLECULE"
              write(*,*)"but you want to project the velocity on a single BOND."
              write(*,*)"You have to chose between V_BOND_REF, SIMPLE_C2V and &HESSIAN"
              write(*,*)"&VVAF"
              write(*,*)"  &VELOCITY0 BOND"
              write(*,*)"    V_BOND_REF / SIMPLE_C2V / &HESSIAN C2V_AX2"
              write(*,*)"    ..."
              write(*,*)"  &END VELOCITY0"
              write(*,*)"  ..."
              write(*,*)"&VVAF"
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           end if


           !=========================================================
           !When the input file has been read twice (nb_read==2)
           !it means the input file is correctly recorded
           !Therefore, the vvaf which is just read can be calculated.
           !=========================================================
           if(nb_read>=2)then
              end_vvaf_block=.false.
              
              !Initialization of tab_mol which contains the molecule information
              !Allocation even if there is no bond (easier tests)
              if(options%v0%nb_bond>0)then
                 allocate(options%v0%tab_mol(options%v0%nb_bond,natoms))
              else
                 allocate(options%v0%tab_mol(1,natoms))
              endif
              if(options%vt%nb_bond>0)then
                 allocate(options%vt%tab_mol(options%vt%nb_bond,natoms))
              else
                 allocate(options%vt%tab_mol(1,natoms))
              endif
              options%v0%tab_mol(:,:)%index=0
              options%vt%tab_mol(:,:)%index=0

              !=========================================================
              !Several VVAF are done wih only one VVAF block if there is
              !an increment for the layering
              !We calculate do a correlation for the volume V1, then 
              !a second correlation is done for the next volume V2 and
              !the results are summed to have the correlation in the
              !volume V1+V2.
              !etc... with the other volumes
              !=========================================================
              !Initialization of the moving limits
              if(options%v0%inc_param)then  !M1-m2 increment will change
                 if(  (options%v0%min1_r>epsilon(0d0)).or.&
                      (options%v0%max1_r>epsilon(0d0)))then
                    options%v0%min1=options%v0%min1_r-options%inc_val
                    options%v0%max1=options%v0%min1_r
                 end if
                 if(  (options%v0%min2_r>epsilon(0d0)).or.&
                      (options%v0%max2_r>epsilon(0d0)))then
                    options%v0%min2=options%v0%max2_r
                    options%v0%max2=options%v0%max2_r+options%inc_val
                 end if
              else
                 if(  (options%v0%min1_r>epsilon(0d0)).or.&
                      (options%v0%max1_r>epsilon(0d0)))then
                    options%v0%min1=options%v0%max1_r
                    options%v0%max1=options%v0%max1_r+options%inc_val
                 end if
                 if(  (options%v0%min2_r>epsilon(0d0)).or.&
                      (options%v0%max2_r>epsilon(0d0)))then
                    options%v0%min2=options%v0%min2_r-options%inc_val
                    options%v0%max2=options%v0%min2_r
                 end if
              endif
              if(options%vt%inc_param)then  !M1-m2 increment will change
                 if(  (options%vt%min1_r>epsilon(0d0)).or.&
                      (options%vt%max1_r>epsilon(0d0)))then
                    options%vt%min1=options%vt%min1_r-options%inc_val
                    options%vt%max1=options%vt%min1_r
                 end if
                 if(  (options%vt%min2_r>epsilon(0d0)).or.&
                      (options%vt%max2_r>epsilon(0d0)))then
                    options%vt%min2=options%vt%max2_r
                    options%vt%max2=options%vt%max2_r+options%inc_val
                 end if
              else
                 if(  (options%vt%min1_r>epsilon(0d0)).or.&
                      (options%vt%max1_r>epsilon(0d0)))then
                    options%vt%min1=options%vt%max1_r
                    options%vt%max1=options%vt%max1_r+options%inc_val
                 end if
                 if(  (options%vt%min2_r>epsilon(0d0)).or.&
                      (options%vt%max2_r>epsilon(0d0)))then
                    options%vt%min2=options%vt%min2_r-options%inc_val
                    options%vt%max2=options%vt%min2_r
                 end if
              endif
                

              do while(.not.end_vvaf_block)
                 end_vvaf_block=.true.
                 !Preparation of the options for the layering
                 !New filename, new min and max
                 if(options%inc)then
                    call update_layer(options%v0,end_vvaf_block,options%inc_val,a,b,c)
                    call update_layer(options%vt,end_vvaf_block,options%inc_val,a,b,c)
                    call update_filename(options%filename,options%filename_tmp,options%v0,options%vt)
                 else
                    options%filename_tmp=options%filename!The filename does not need to be modified
                    options%v0%min1=options%v0%min1_r
                    options%v0%max1=options%v0%max1_r
                    options%v0%min2=options%v0%min2_r
                    options%v0%max2=options%v0%max2_r
                    options%vt%min1=options%vt%min1_r
                    options%vt%max1=options%vt%max1_r
                    options%vt%min2=options%vt%min2_r
                    options%vt%max2=options%vt%max2_r
                 endif
                 
                 !call cpu_time(start)
                 call vvaf(nmovie,natoms,Di,shift,options,vx(1,1),vy(1,1),vz(1,1),&
                      x(1,1),y(1,1),z(1,1),slab_orient,center,a,b,c,box_size,&
                      surface,corr)
                 !call cpu_time(finish)
                 !write(*,*)"Time=",finish-start
              enddo

              !This arrays are allocated even if nb_bond=0
              deallocate(options%v0%tab_mol,options%vt%tab_mol)

           endif

           !The options derived type are deleted and rewrote at each VVAF block 
           deallocate(options%v0%tab_at_bond,options%vt%tab_at_bond,&
                      options%v0%tab_select ,options%vt%tab_select,&
                      options%v0%tab_center ,options%vt%tab_center)
           exit
        endif
        
     case('')
        
     case default
        call wrong_keyword(chain,block_name)

     end select
  enddo

end subroutine input_vvaf


!================================================
!Subroutine which reads the block VVAF > VELOCITY
!It contains the information about the velocity
!1 or 2 to take into account
!================================================
subroutine input_vvaf_velocity(nb_vel,filepos,natoms,vel,inc,cross,atom_type,fmass,a,b,c)
  use mymod, only : velocity,read_word,wrong_keyword,input_vvaf_velocity_hessian,&
       input_vvaf_velocity_atom,read_tab_pol,read_tab_dipol,error
  implicit none

  !From upper subroutine
  integer natoms,nb_vel
  double precision a,b,c
  double precision, dimension(natoms) :: fmass
  character(len=3), dimension(natoms) :: atom_type
  character(len=100) filepos
  logical :: inc,cross
  type (velocity) :: vel

  !Local
  integer:: i,j,end_file
  double precision:: x
  double precision, dimension(2) :: m
  character(len=100):: chain,block_name
  logical next_new

  
  !Initialization
  i=0 ; j=0 ; end_file=0
  x=0.d0
  m(:)=0.d0
  chain=''
  next_new=.false.
  if(nb_vel==1)then
     block_name='VELOCITY0'
  else
     block_name='VELOCITYT'
  endif

  !Do we study the velocity of the reference atom or
  !The velocity of the reference atom minus those of the bonded atom?
  !There is also the possibility to do specifically an autocorrelation
  call read_word(chain,end_file,next_new)
  if(chain=='BOND')then
     vel%bond_diff=.true.
  elseif(chain=='ATOM') then
     vel%bond_diff=.false.
  elseif(chain=='AUTO') then
     cross=.false.
     if(block_name=='VELOCITY0')then
        write(*,*)"AUTO works only with VELOCITYT"
        call wrong_keyword(chain,block_name)
     endif
  else
     call wrong_keyword(chain,block_name)
  endif

  !Loop which read the block VVAF > VELOCITY
  do while (end_file<=0)
     call read_word(chain,end_file,next_new)

     select case(chain)
     case('&ATOM')
        if(.not.cross.and.(block_name=='VELOCITYT'))then
           write(*,*)"&ATOM block to remove because auto-correlation requiered"
           write(*,*)"&VVAF"
           write(*,*)"  &VELOCITYT AUTO"
           write(*,*)"  &END VELOCITYT"
           write(*,*)"&END VVAF"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        else
           call input_vvaf_velocity_atom(natoms,vel,atom_type,a,b,c,inc)
        endif

     case('NOSYM_STRETCH')!No stretching mode study (independent bonds)
        vel%c2v_mode=0

     case('SIMPLE_C2V')!Study simple/doom/uncorrect bending / antisymmetric / symmetric stretching of water molecules (c2v symmetry)
        vel%c2v_mode=1

     case("&HESSIAN")!Study only the bending / antisymmetric / symmetric stretching of water molecules has defined by the hessian
        vel%c2v_mode=2
        call input_vvaf_velocity_hessian(vel)
        
     case('V_BOND_REF')
        vel%scal_proj=.true.

     case('VEL_PROJ')
        call read_word(chain,end_file,next_new)
        if(chain=='n')then
           vel%vproj=0
        elseif(chain=='X')then
           vel%vproj=1
        elseif(chain=='Y')then
           vel%vproj=2
        elseif(chain=='Z')then
           vel%vproj=3
        else
           write(*,*)"Velocity projection not well defined"
           write(*,*)"&VVAF"
           write(*,*)"  VEL_PROJ 'X'/'Y'/'Z'/'n'"
           write(*,*)"&END VVAF"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        endif

     case('&DIPOLE')
        call read_word(chain,end_file,next_new)
        select case (chain(1:1))
        case ('n')
           vel%R=0
        case ('X')
           vel%R=1
        case ('Y')
           vel%R=2
        case ('Z')
           vel%R=3
        case default
           write(*,*)"Dipole orientation not well defined"
           write(*,*)"&VVAF"
           write(*,*)"  &DIPOLE 'X'/'Y'/'Z'/'n'"
           write(*,*)"    dpx/dpr   dpy/dpr   dpz/dpr"
           write(*,*)"  &END DIPOLE"
           write(*,*)"&END VVAF"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        end select
        call read_tab_dipol(vel,next_new)

     case('&POLARIZABILITY')
        do i=1,2!To define the polarizability, 2 parameters are requiered
           call read_word(chain,end_file,next_new)
           select case (chain(1:1))
           case ('n')
              if(i==1)then
                 vel%P=0
              else 
                 vel%Q=0
              endif
           case ('X')
              if(i==1)then
                 vel%P=1
              else 
                 vel%Q=1
              endif
           case ('Y')
              if(i==1)then
                 vel%P=2
              else 
                 vel%Q=2
              endif
           case ('Z')
              if(i==1)then
                 vel%P=3
              else 
                 vel%Q=3
              endif
           case default
              write(*,*)"Polarizability not well defined"
              write(*,*)"&VVAF"
              write(*,*)"  POLARIZABILITY 'X'/'Y'/'Z'/'n' 'X'/'Y'/'Z'/'n'"
              write(*,*)"&END VVAF"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           end select
        enddo
        !If among the parameters, one is 0 (nothing), it does not work
        if(  ((vel%P/=0).and.(vel%Q==0)).or.&
             ((vel%Q/=0).and.(vel%P==0)))then
           write(*,*)"Two axis have to be defined for the polarizability 'X n', etc... is forbidden."
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        else !If the two parameters are well defined, then we record the the derivate polarizability matrix 
           !The data have to be written like a symmetric matix
           !Some test have to be done
           call read_tab_pol(vel,next_new)
        endif

     case('STRETCH_PROJ')
        call read_word(chain,end_file,next_new)
        if(chain=='y')then
           vel%stretch_mode=.true.
        elseif(chain=='n')then
           vel%stretch_mode=.false.
        else
           write(*,*)"Stretching projection not well defined"
           write(*,*)"&VVAF"
           write(*,*)"  STRETCH_PROJ 'y'/'n'"
           write(*,*)"&END VVAF"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        endif

     case('&END')
        call read_word(chain,end_file,next_new)!&END as to be followed by the name of the block
        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        else
           !Before to leave the VELOCITY block, some verification have to be done
           
           !If some informations about bonds are requiered, the position file as to be specified!!!!
           !It is also necessary to have at least one bond (vel%tab_at_bond=.true.)
           !BUT for this case I do the test at the end of the VVAF block 
           !because I do not want to impose an order between VELOCITY0 and VELOCITYT
           !even if the AUTO flag is added after VELOCITYT
           if((vel%c2v_mode>0).or.(vel%P/=0).or.(vel%R/=0).or.vel%stretch_mode.or.vel%bond_diff)then
              if(filepos=='')then
                 write(*,*)"You need to define some bonds, but you have no position file."
                 write(*,*)"&MAIN"
                 write(*,*)"  POS_FILE 'name' "
                 write(*,*)"&END MAIN"
                 write(*,*)
                 write(*,*)"!!! END OF PROGRAM !!!"
                 call exit()
                 
              endif
           endif

           
           !========================================
           !Have we got normal modes?
           !sqrt(M)*r_i . sqrt(M)*r_j = 0 if i/=j
           !The test are done here, because we need
           !to have some info from the ATOM block
           !Moreover, (sqrt(M)*r_i)^2 may not be 1
           !For the user, it has to be relevant with
           !the variations given in the DIPOLE and
           !the POLARIZABILITY blocks.
           !The normalization which is important is:
           !p_i=M r_i/(sqrt(M)*r_i . sqrt(M)*r_i)
           !Like this we will have:
           !p_i.r=contribution of r to the i-th mode
           !========================================
           !1) We have to check the mass of the different species
           !We expect to have AX2, therefore we have to be sure that
           !all the A (X) have the same mass (necessary if we use
           !index instead of symbol)
           if(vel%c2v_mode>0) then
              do i=1,natoms
                 if(vel%tab_select(i))then
                    m(1)=fmass(i)
                    exit
                 end if
              end do
              do j=i,natoms
                 if(vel%tab_select(j))then
                    if( (m(1)-fmass(j))**2 > error )then
                       write(*,*)"Your central atoms do not have the same mass."
                       write(*,*)"Therefore, you cannot use the the Hessian"
                       write(*,*)"definition of the normal modes"
                       write(*,*)
                       write(*,*)"!!! END OF PROGRAM !!!"
                       call exit()
                    end if
                 end if
              end do

              do i=1,natoms
                 if(vel%tab_at_bond(i))then
                    m(2)=fmass(i)
                    exit
                 end if
              end do
              do j=i,natoms
                 if(vel%tab_at_bond(j))then
                    if( (m(2)-fmass(j))**2 > error )then
                       write(*,*)"Your ligands do not have the same mass."
                       write(*,*)"Therefore, you cannot use the the Hessian"
                       write(*,*)"definition of the normal modes."
                       write(*,*)
                       write(*,*)"!!! END OF PROGRAM !!!"
                       call exit()
                    end if
                 end if
              end do
              !Security to not have too small numbers
              !We will use these coefficients just after and will see
              !if some scalar products are equal to 0
              if(m(1)<m(2))then
                 m(2)=m(2)/m(1)
                 m(1)=1.d0
              else
                 m(1)=m(1)/m(2)
                 m(2)=1.d0
              end if

              !2) If i!=j, the scalar product has to be equal to 0
              do i=1,3
                 do j=i+1,3

                    if((&
                         m(1) * vel%molf2q(1,1,i) * vel%molf2q(1,1,j) + &
                         m(1) * vel%molf2q(2,1,i) * vel%molf2q(2,1,j) + &
                         m(1) * vel%molf2q(3,1,i) * vel%molf2q(3,1,j) + &
                         m(2) * vel%molf2q(1,2,i) * vel%molf2q(1,2,j) + &
                         m(2) * vel%molf2q(2,2,i) * vel%molf2q(2,2,j) + &
                         m(2) * vel%molf2q(3,2,i) * vel%molf2q(3,2,j) + &
                         m(2) * vel%molf2q(1,3,i) * vel%molf2q(1,3,j) + &
                         m(2) * vel%molf2q(2,3,i) * vel%molf2q(2,3,j) + &
                         m(2) * vel%molf2q(3,3,i) * vel%molf2q(3,3,j)   &
                         )**2>error)then
                       write(*,*)"Your modes are not normal"
                       write(*,*)
                       write(*,*)"!!! END OF PROGRAM !!!"
                       call exit()
                    end if
                 end do
              end do

              !3) Normalization: p_i=M r_i/(sqrt(M)*r_i . sqrt(M)*r_i)
              !Each time the water molecule moves of 1 unit
              !(defined by the Hessian), the dipole/polarizability
              !increases of 1 unit (defined by the DIPOLE/POLARIZABILITY
              !block)
              do i=1,3 !Loop of the modes
                 x=   m(1) * vel%molf2q(1,1,i)**2 + &
                      m(1) * vel%molf2q(2,1,i)**2 + &
                      m(1) * vel%molf2q(3,1,i)**2 + &
                      m(2) * vel%molf2q(1,2,i)**2 + &
                      m(2) * vel%molf2q(2,2,i)**2 + &
                      m(2) * vel%molf2q(3,2,i)**2 + &
                      m(2) * vel%molf2q(1,3,i)**2 + &
                      m(2) * vel%molf2q(2,3,i)**2 + &
                      m(2) * vel%molf2q(3,3,i)**2

                 vel%molf2q(:,1,i)=m(1)*vel%molf2q(:,1,i)/x !A
                 vel%molf2q(:,2,i)=m(2)*vel%molf2q(:,2,i)/x !X1
                 vel%molf2q(:,3,i)=m(2)*vel%molf2q(:,3,i)/x !X2
              end do
           end if

           exit !end of the while loop and therefore the subroutine
        endif

     case('')

     case default
        call wrong_keyword(chain,block_name)

     end select

  enddo
  
end subroutine input_vvaf_velocity



!=====================================================
!Subroutine which records the derivative dipole vector
!x y z
!=====================================================
subroutine read_tab_dipol(vel,next_new)
  use mymod, only : velocity,read_word,wrong_keyword,skip_line_tab_dipol,&
       error_deriv_dip_pol
  implicit none

  !From upper subroutine
  logical next_new
  type (velocity) :: vel

  !Local variables
  integer:: i,end_file,mode
  character(len=100) :: chain,block_name

  !Initialization
  i=0 ; end_file=0 ; mode=0
  chain=''
  block_name='DIPOLE'

  !A line is skipped after the keyword &DIPOLE
  call skip_line_tab_dipol(next_new)!Skip line

  !Dipole derivative vectors
  call read_word(chain,end_file,next_new)
  do while(chain/='&END')
     select case (chain)
     case("BEND")
        if(mode==1)call error_deriv_dip_pol()
        mode=2
        vel%c2v_deriv=.true.
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%dip_bend !Along the bissector (z)
     case("AS")
        if(mode==1)call error_deriv_dip_pol()
        mode=2
        vel%c2v_deriv=.true.
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%dip_as !Along perp to the bissector but in the normal plane (x)
     case("SS")
        if(mode==1)call error_deriv_dip_pol()
        mode=2
        vel%c2v_deriv=.true.
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%dip_ss !Along the bissector (z)
     case("STRETCH_x")
        if(mode==2)call error_deriv_dip_pol()
        mode=1
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%dip_stretch(1) !In the molecular plane
     case("STRETCH_y")
        if(mode==2)call error_deriv_dip_pol()
        mode=1
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%dip_stretch(2) !Out of the molecular plane
     case("STRETCH_z")
        if(mode==2)call error_deriv_dip_pol()
        mode=1
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%dip_stretch(3) !Along the bond
     case("")
        call read_word(chain,end_file,next_new)
     case default
        call wrong_keyword(chain,block_name)
     end select

     call skip_line_tab_dipol(next_new)!Skip line
     call read_word(chain,end_file,next_new)
  enddo

  !&END DIPOLE
  call read_word(chain,end_file,next_new)
  if(chain/='DIPOLE')then
     call wrong_keyword(chain,block_name)
  endif
  
end subroutine read_tab_dipol

!============================================================================
!Subroutine which checks if the derivative dipole vector is correctly written
!============================================================================
subroutine skip_line_tab_dipol(next_new)
  use mymod, only : read_word
  implicit none

  !From upper subroutine
  logical next_new

  !Local variables
  integer:: end_file
  character(len=100) :: chain

  !Initialization
  end_file=0
  chain=''
  
  !A line has to be skipped (either the END OF THE LINE is already reached
  !or we have to continue to read the BLANKS until the end of the line)
  if(.not.next_new)then
     call read_word(chain,end_file,next_new)
     if(chain/=' ')then
        !The derivative polarization matrix is not correctly written
        write(*,*)"The derivative dipole vector is not correctly written"
        write(*,*)"&VVAF"
        write(*,*)"  &VELOCITY?"
        write(*,*)"    &DIPOLE ? (X/Y/Z)"
        write(*,*)"      BEND MU_z^MOL"
        write(*,*)"      SS   MU_z^MOL"
        write(*,*)"      AS   MU_x^MOL"
        write(*,*)" or "
        write(*,*)"      STRETCH_x  MU_x^BOND"
        write(*,*)"      STRETCH_y  MU_y^BOND"
        write(*,*)"      STRETCH_z  MU_z^BOND"
        write(*,*)"      (ommited lines are set to 0)"
        write(*,*)"    &END DIPOLE"
        write(*,*)"  &END VELOCITY?"
        write(*,*)"&END VVAF"
        write(*,*)
        write(*,*)"!!! END OF PROGRAM !!!"
        call exit()
     endif
  endif
  call read_word(chain,end_file,next_new)!The line is skipped
end subroutine skip_line_tab_dipol


!===========================================================
!Subroutine which records the derivative polarization matrix
!xx
!yx yy
!zx zy zz
!===========================================================
subroutine read_tab_pol(vel,next_new)
  use mymod, only : velocity,read_word,wrong_keyword,skip_line_tab_pol,&
       error_deriv_dip_pol
  implicit none

  !From upper subroutine
  type (velocity) :: vel
  logical next_new

  !Local variables
  integer:: i,j,end_file,mode
  character(len=100) :: chain,block_name

  !Initialization
  i=0; j=0; end_file=0; mode=0
  chain=''
  block_name='POLARIZABILITY'

  !A line is skipped after the keyword &POLARIZABILITY
  call skip_line_tab_pol(next_new)!Skip line

  !Polarizability derivative tensors
  call read_word(chain,end_file,next_new)
  do while(chain/='&END')
     select case (chain)
     case("BEND")
        if(mode==1)call error_deriv_dip_pol()
        mode=2
        vel%c2v_deriv=.true.
        !xx yy zz
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_bend(1)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_bend(2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_bend(3)
     case("AS")
        if(mode==1)call error_deriv_dip_pol()
        mode=2
        vel%c2v_deriv=.true.
        !xz
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_as
     case("SS")
        if(mode==1)call error_deriv_dip_pol()
        mode=2
        vel%c2v_deriv=.true.
        !xx yy zz
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_ss(1)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_ss(2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_ss(3)
     case("STRETCH_x")
        if(mode==2)call error_deriv_dip_pol()
        mode=1
        !xx
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_stretch(1,1)
     case("STRETCH_y")
        if(mode==2)call error_deriv_dip_pol()
        mode=1
        !yx yy
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_stretch(1,2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_stretch(2,2)
     case("STRETCH_z")
        if(mode==2)call error_deriv_dip_pol()
        mode=1
        !zx zy zz
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_stretch(1,3)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_stretch(2,3)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%pol_stretch(3,3)
     case("")
        call read_word(chain,end_file,next_new)
     case default
        call wrong_keyword(chain,block_name)
     end select

     call skip_line_tab_pol(next_new)!Skip line
     call read_word(chain,end_file,next_new)
  enddo

  !&END POLARIZABILITY
  call read_word(chain,end_file,next_new)
  if(chain/='POLARIZABILITY')then
     call wrong_keyword(chain,block_name)
  endif

  !The upper part of the symmetric matrix is filled
  vel%pol_stretch(2,1)=vel%pol_stretch(1,2)
  vel%pol_stretch(3,1)=vel%pol_stretch(1,3)
  vel%pol_stretch(3,2)=vel%pol_stretch(2,3)

end subroutine read_tab_pol



!==========================================================================
!Subroutine which skips lines when one reads derivative polarization matrix
!==========================================================================
subroutine skip_line_tab_pol(next_new)
  use mymod, only : read_word
  implicit none

  !From upper subroutine
  logical next_new

  !Local variables
  integer:: end_file
  character(len=100) :: chain

  !Initialization
  end_file=0
  chain=''
  
  !A line has to be skipped (either the END OF THE LINE is already reached
  !or we have to continue to read the BLANKS until the end of the line)
  if(.not.next_new)then
     call read_word(chain,end_file,next_new)
     if(chain/=' ')then
        !The derivative polarization matrix is not correctly written
        write(*,*)"The derivative polarizability matrix is not correctly written"
        write(*,*)"&VVAF"
        write(*,*)"  &VELOCITY?"
        write(*,*)"    &POLARIZABILITY ? ?"
        write(*,*)"      BEND xx^MOL yy^MOL zz^MOL"
        write(*,*)"      SS   xx^MOL yy^MOL zz^MOL"
        write(*,*)"      AS   xz^MOL"
        write(*,*)" or "
        write(*,*)"      STRETCH_x ALPHA_xx^bond"
        write(*,*)"      STRETCH_y ALPHA_yx^bond yy"
        write(*,*)"      STRETCH_z ALPHA_zx^bond zy zz"
        write(*,*)"      (ommited lines are set to 0)"
        write(*,*)"    &END POLARIZABILITY"
        write(*,*)"  &END VELOCITY?"
        write(*,*)"&END VVAF"
        write(*,*)
        write(*,*)"!!! END OF PROGRAM !!!"
        call exit()
     endif
  endif
  call read_word(chain,end_file,next_new)!The line is skipped
end subroutine skip_line_tab_pol


!================================================
!Subroutine which send an error message
!if you use at the same time the derivatives
!for the stretching and those of the c2v symmetry
!================================================
subroutine error_deriv_dip_pol()
  write(*,*)""
  write(*,*)"You cannot use at the same time the derivatives"
  write(*,*)"for the stretching (bond frame) and those for "
  write(*,*)"the c2v derivatives"
  write(*,*)
  write(*,*)"STRETCH_x"
  write(*,*)"STRETCH_y"
  write(*,*)"STRETCH_z"
  write(*,*)"or"
  write(*,*)"BEND"
  write(*,*)"SS"
  write(*,*)"AS"
  write(*,*)
  write(*,*)"!!! END OF PROGRAM !!!"
  call exit()

end subroutine error_deriv_dip_pol

!======================================================
!Loop which read the block VVAF > VELOCITY > ATOM
!All the information about the reference atom are there
!======================================================
subroutine input_vvaf_velocity_atom(natoms,vel,atom_type,a,b,c,inc)
  use mymod, only : velocity,read_word,wrong_keyword,input_vvaf_velocity_atom_bondw,&
       input_vvaf_velocity_atom_layer,input_vvaf_velocity_atom_annulus
  implicit none

  !From upper subroutine
  integer natoms
  double precision a,b,c
  character(len=3), dimension(natoms) :: atom_type
  logical inc
  type (velocity) :: vel

  !Local
  integer:: i,end_file
  character(len=100):: chain,block_name
  logical next_new

  !Initialization
  i=0
  end_file=0
  chain=''
  block_name='ATOM'
  next_new=.false.

  do while (end_file<=0)
     call read_word(chain,end_file,next_new)  

     select case(chain)
     case('&BOND_WITH')
        call input_vvaf_velocity_atom_bondw(natoms,vel,atom_type)

     case('&LAYER')!Study of the atoms in a specific layer
        call input_vvaf_velocity_atom_layer(vel,a,b,c,inc)

     case('&ANNULUS')
        call input_vvaf_velocity_atom_annulus(natoms,vel,atom_type)

     case('INDEX')!The atoms we the good index are selected
        do while(.not. next_new)!-1: end of word but not end of line
           call read_word(chain,end_file,next_new)
           if(chain/='')read(chain,*)i
           vel%tab_select(i)=.true.
        enddo

     case('SYMBOL')
        do while(.not. next_new)
           call read_word(chain,end_file,next_new)
           do i=1,natoms
              if(atom_type(i).eq.chain)vel%tab_select(i)=.true.
           enddo
        enddo 

     case('ALL')
        vel%tab_select=.true.

     case('&END')
        call read_word(chain,end_file,next_new)!&END as to be followed by the name of the block
        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        else
           !Before to leave the ATOM block, some verifications have to be done

           !To do a VVAF, some atoms have to be taken into account
           if(.not.any(vel%tab_select))then
              write(*,*)"No atoms are selected"
              write(*,*)"&VVAF"
              write(*,*)"  &VELOCITY?"
              write(*,*)"    &ATOM"
              write(*,*)"      INDEX 1 2 3 / SYMBOL H O / ALL"
              write(*,*)"    &END ATOM"
              write(*,*)"  &END VELOCITY?"
              write(*,*)"&VVAF"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           endif

           exit

        endif

     case('')

     case default
        call wrong_keyword(chain,block_name)
     endselect
  enddo
end subroutine input_vvaf_velocity_atom



!============================================================
!Loop which read the block VVAF > VELOCITY > ATOM > BOND_WITH
!The reference atom can make some bonds with another atom
!============================================================
subroutine input_vvaf_velocity_atom_bondw(natoms,vel,atom_type)
  use mymod, only : velocity,wrong_keyword,read_word
  implicit none

  !From upper subroutine
  integer natoms
  character(len=3), dimension(natoms) :: atom_type
  type (velocity) :: vel

  !Local
  integer:: i,end_file
  character(len=100):: chain,block_name
  logical next_new

  !Initialization
  i=0
  end_file=0
  chain=''
  block_name='BOND_WITH'
  next_new=.false.


  do while (end_file<=0)
     call read_word(chain,end_file,next_new)  
     select case(chain)
     !Atom selection
     case('INDEX')
        do while(.not. next_new)
           call read_word(chain,end_file,next_new)
           if(chain/='')read(chain,*)i
           vel%tab_at_bond(i)=.true.
        enddo

     case('SYMBOL')
        do while(.not. next_new)
           call read_word(chain,end_file,next_new)
           do i=1,natoms
              if(atom_type(i).eq.chain)vel%tab_at_bond(i)=.true.
           enddo
        enddo 

     case('ALL')
        vel%tab_at_bond=.true.

     !Bond characteristics
     case('NB_BOND')
        call read_word(chain,end_file,next_new)
        if(chain==' ')then
           write(*,*)"How many atoms make bonds with the central atom?"
           write(*,*)"&VVAF"
           write(*,*)"  &VELOCITY?"
           write(*,*)"    &ATOM"
           write(*,*)"      &BOND_WITH"
           write(*,*)"        NB_BOND   1/2/3/..."
           write(*,*)"      &END BOND_WITH"
           write(*,*)"    &END ATOM"
           write(*,*)"  &END VELOCITY?"
           write(*,*)"&END VVAF"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        else
           read(chain,*)vel%nb_bond
        endif

     case('DMIN')
        call read_word(chain,end_file,next_new)
        if(chain==' ')then
           write(*,*)"How many atoms make bonds with the central atom?"
           write(*,*)"&VVAF"
           write(*,*)"  &VELOCITY?"
           write(*,*)"    &ATOM"
           write(*,*)"      &BOND_WITH"
           write(*,*)"        DMIN   ???"
           write(*,*)"      &END BOND_WITH"
           write(*,*)"    &END ATOM"
           write(*,*)"  &END VELOCITY?"
           write(*,*)"&END VVAF"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        else
           read(chain,*)vel%dmin2
           vel%dmin2=vel%dmin2**2
        endif

     case('DMAX')
        call read_word(chain,end_file,next_new)
        if(chain==' ')then
           write(*,*)"How many atoms make bonds with the central atom?"
           write(*,*)"&VVAF"
           write(*,*)"  &VELOCITY?"
           write(*,*)"    &ATOM"
           write(*,*)"      &BOND_WITH"
           write(*,*)"        DMAX   ???"
           write(*,*)"      &END BOND_WITH"
           write(*,*)"    &END ATOM"
           write(*,*)"  &END VELOCITY?"
           write(*,*)"&END VVAF"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        else
           read(chain,*)vel%dmax2
           vel%dmax2=vel%dmax2**2
        endif

     !End
     case('&END')
        call read_word(chain,end_file,next_new)!&END as to be followed by the name of the block
        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        else
           exit
        endif

     case('')

     case default
        call wrong_keyword(chain,block_name)
     endselect

  end do

end subroutine input_vvaf_velocity_atom_bondw


!========================================================
!Loop which read the block VVAF > VELOCITY > ATOM > LAYER
!Only the atoms which are in a specific layer (or in two 
!specifi layers are taken into account
!========================================================
subroutine input_vvaf_velocity_atom_layer(vel,a,b,c,inc)
  use mymod, only : velocity,wrong_keyword,recenter,read_word
  implicit none

  !From upper subroutine
  double precision a,b,c
  logical inc
  type (velocity) :: vel

  !Local
  integer:: end_file
  character(len=100):: chain,block_name
  logical next_new

  !Initialization
  end_file=0
  chain=''
  block_name='LAYER'
  next_new=.false.

  !We cannot select at the same time LAYER and ANNULUS
  if(vel%vol_meth=='r')then
     write(*,*)"ANNULUS already selected"
     write(*,*)"LAYER cannot be used"
     write(*,*)
     write(*,*)"!!! END OF PROGRAM !!!"
     call exit()
  endif

  !Normal axis to the layer
  call read_word(chain,end_file,next_new)
  if((chain/='X').and.(chain/='Y').and.(chain/='Z'))then
     write(*,*)"Axis of the layer not defined"
     write(*,*)"&LAYER 'X'/'Y'/'Z'"
     write(*,*)
     write(*,*)"!!! END OF PROGRAM !!!"
     call exit()
  endif
  vel%vol_meth=chain(1:1)

  do while (end_file<=0)
     call read_word(chain,end_file,next_new)  

     select case(chain)
     case('SLAB1')
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%min1_r
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%max1_r

     case('SLAB2')
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%min2_r
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%max2_r

     case('PARAM_INC')
        !An increment will be associated with the layering
        inc=.true.

        !The parameters which will change will be...
        call read_word(chain,end_file,next_new)
        if(chain=='m1-M2')then
           vel%inc_param=.false.
        elseif(chain=='M1-m2')then
           vel%inc_param=.true.
        else
           write(*,*)"Wrong keyword after PARAM_INC"
           write(*,*)"PARAM_INC M1-m2"
           write(*,*)"or"
           write(*,*)"PARAM_INC m1-M2"
           write(*,*)
           write(*,*)"!!! END OF PROGRAM !!!"
           call exit()
        endif

     case('&END')
        call read_word(chain,end_file,next_new)!&END as to be followed by the name of the block
        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        else
           !Before to leave some controls have to be done
           !Centering centering the layers according to the pbc
           if(trim(vel%vol_meth)=='X')then
              call recenter(vel,a)
           elseif(trim(vel%vol_meth)=='Y')then
              call recenter(vel,b)
           elseif(trim(vel%vol_meth)=='Z')then
              call recenter(vel,c)
           endif

           exit
        endif

     case('')

     case default
        call wrong_keyword(chain,block_name)
     end select
  enddo

end subroutine input_vvaf_velocity_atom_layer



!==========================================================
!Loop which read the block VVAF > VELOCITY > ATOM > ANNULUS
!Only the atoms which are in a specific annulus are taken
! into account
!==========================================================
subroutine input_vvaf_velocity_atom_annulus(natoms,vel,atom_type)
  use mymod, only : velocity,wrong_keyword,read_word
  implicit none

  !From upper subroutine
  integer natoms
  character(len=3), dimension(natoms) :: atom_type
  type (velocity) :: vel

  !Local
  integer:: i=0,end_file=0
  character(len=100):: chain='',block_name='ANNULUS'
  logical next_new

  !Initialization
  i=0
  end_file=0
  chain=''
  block_name='ANNULUS'
  next_new=.false.
  allocate(vel%tab_center(natoms))
  vel%tab_center=.false.

  !If the layering is not taken into account we can have a annulus study.
  if(vel%vol_meth/='n')then
     write(*,*)"LAYER already selected"
     write(*,*)"ANNULUS cannot be used"
     write(*,*)
     write(*,*)"!!! END OF PROGRAM !!!"
     call exit()
  else
     vel%vol_meth='r'
  endif
  

  do while (end_file<=0)
     call read_word(chain,end_file,next_new)  

     select case(chain)
     case('INDEX')
        do while(.not. next_new)
           call read_word(chain,end_file,next_new)
           if(chain/='')read(chain,*)i
           vel%tab_center(i)=.true.
        enddo

     case('SYMBOL')
        do while(.not. next_new)
           call read_word(chain,end_file,next_new)
           do i=1,natoms
              if(atom_type(i).eq.chain)vel%tab_center(i)=.true.
           enddo
        enddo 

     case('V0_AT')
           vel%v0_at=.true.

     case('ALL')
        vel%tab_center=.true.

     case('RADIUS')
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%min1_r
        call read_word(chain,end_file,next_new)
        read(chain,*)vel%max1_r

     case('&END')
        call read_word(chain,end_file,next_new)!&END as to be followed by the name of the block
        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        else
           !Before to leave ANNULUS block, some tests have to be done
           
           !Incompatible options
           if(any(vel%tab_center).and.vel%v0_at)then
              write(*,*)"Incompatible options."
              write(*,*)"Either 'V0_AT' or 'INDEX 1 2 3 / SYMBOL H O / ALL'"
              write(*,*)
              write(*,*)"&VVAF"
              write(*,*)"  &VELOCITY?"
              write(*,*)"    &ATOM"
              write(*,*)"      &ANNULUS"
              write(*,*)"        V0_AT / INDEX 1 2 3 / SYMBOL H O / ALL"
              write(*,*)"      &END ANNULUS"
              write(*,*)"    &END ATOM"
              write(*,*)"  &END VELOCITY?"
              write(*,*)"&VVAF"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           elseif((.not.any(vel%tab_center)).and.(.not.vel%v0_at))then
              !At least one center has to be defined
              write(*,*)"You need to define at least one center. E.g."
              write(*,*)"&VVAF"
              write(*,*)"  &VELOCITY?"
              write(*,*)"    &ATOM"
              write(*,*)"      &ANNULUS"
              write(*,*)"        V0_AT / INDEX 1 2 3 / SYMBOL H O / ALL"
              write(*,*)"      &END ANNULUS"
              write(*,*)"    &END ATOM"
              write(*,*)"  &END VELOCITY?"
              write(*,*)"&VVAF"
              write(*,*)
              write(*,*)"!!! END OF PROGRAM !!!"
              call exit()
           endif

           exit
        endif

     case('')

     case default
        call wrong_keyword(chain,block_name)
     end select
  enddo

end subroutine input_vvaf_velocity_atom_annulus


!==================================
!Subroutine which reads the Hessian
!==================================
subroutine input_vvaf_velocity_hessian(vel)
  use mymod, only : velocity,wrong_keyword,read_word,error
  implicit none

  !From upper subroutine
  type (velocity) :: vel

  !Local
  integer:: end_file, i, j
  character(len=100):: chain,block_name
  logical next_new

  !Initialization
  end_file=0 ; i=0 ; j=0
  next_new=.false.
  chain='' ; block_name='HESSIAN'
  
  call read_word(chain,end_file,next_new)
  if(chain/="C2V_AX2")then
     write(*,*)"You have to precise that you want to study the Hessian"
     write(*,*)" of an AX2 molecule with a C2V symmetry. Today only"
     write(*,*)"such symmetry are available (e.g. H2O) but I do not"
     write(*,*)"know what will be the next updates."
     write(*,*)"&VVAF"
     write(*,*)"  &VELOCITY?"
     write(*,*)"    &HESSIAN C2V_AX2"
     write(*,*)"    &END HESSIAN"
     write(*,*)"  &END VELOCITY?"
     write(*,*)"&VVAF"
     write(*,*)
     write(*,*)"!!! END OF PROGRAM !!!"
     call exit()
  end if

  
  do while (end_file<=0)
     call read_word(chain,end_file,next_new)  

     select case(chain)

     case('BEND_A')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,1,1)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,1,1)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,1,1)

     case('BEND_Xl')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,2,1)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,2,1)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,2,1)

     case('BEND_Xr')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,3,1)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,3,1)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,3,1)

     case('AS_A')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,1,2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,1,2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,1,2)

     case('AS_Xl')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,2,2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,2,2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,2,2)

     case('AS_Xr')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,3,2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,3,2)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,3,2)

     case('SS_A')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,1,3)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,1,3)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,1,3)

     case('SS_Xl')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,2,3)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,2,3)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,2,3)

     case('SS_Xr')
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(1,3,3)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(2,3,3)
        call read_word(chain,end_file,next_new) ; read(chain,*)vel%molf2q(3,3,3)

     case('&END')
        call read_word(chain,end_file,next_new)!&END as to be followed by the name of the block
        if(chain/=block_name)then
           call wrong_keyword(chain,block_name)
        else
           !Before to leave this block some tests have to be done
           !One has to check if the given values are in agreement with the
           !c2v symmetry
           do i=1,3
              do j=1,3
                 if(vel%molf2q(2,j,i)**2>error)then
                    write(*,*)"For c2v symmetry and AX2 molecule in the molecular"
                    write(*,*)"framework, there is no mode having contribution"
                    write(*,*)"along y"
                    write(*,*)
                    write(*,*)"!!! END OF PROGRAM !!!"
                    call exit()
                 end if
              end do
           end do

           do i=1,3
              do j=1,3
                 if(vel%molf2q(j,2,i)**2-vel%molf2q(j,3,i)**2>error)then
                    write(*,*)"For c2v symmetry and AX2 molecule, the two X "
                    write(*,*)"should have symmetric properties"
                    write(*,*)
                    write(*,*)"!!! END OF PROGRAM !!!"
                    call exit()
                 end if
              end do
           end do

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !At the end of the VELOCITY block, we test if the
           !modes are normal.
           !It is not done here, because we need the information
           !about the AX2 atoms (mass) which are obtained in
           !another block.
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           exit
        end if

     case('')

     case default
        call wrong_keyword(chain,block_name)
     end select
  enddo


  
end subroutine input_vvaf_velocity_hessian
        



!================================================================
!Subroutine which stops the program when there is a wrong keyword
!================================================================
subroutine wrong_keyword(chain,block_name)
  implicit none
  !From upper subroutine
  character(len=100) block_name, chain

  !ARBEIT WORK TRAVAIL
  !I HAVE TO ADD THE LINE NUMBER

  write(*,*)"Error in the block ",trim(block_name)
  write(*,*)"Problematic keyword:",trim(chain)
  write(*,*)
  write(*,*)"!!! END OF PROGRAM !!!"
  call exit()
end subroutine wrong_keyword



!=========================================
!Read lines word by word (space separator)
!and return some information about iostat
!If a line is skipped, then word=''
!=========================================
subroutine read_word(word,err,next_new)
  implicit none
  !From call subroutine
  character(len=100) word
  integer err
  logical next_new

  !Local
  integer:: i

  !Initialization
  err=0
  word=''
  i=1


  !If previously a line has been skipped, then word=''
  if(next_new)then
     err=-2
     word=''
     next_new=.false.

  else!On the contrary, we read word
     !Passing line
     do while((err==0).and.((word(i:i)==' ').or.(word(i:i)==achar(9))))
        read(10,'(a1)',advance='no',iostat=err)word(i:i)

        if(err==-2)then
           next_new=.true.
        elseif(word(i:i)=='#')then
           read(10,*)
           word(i:i)=' '
           next_new=.true.
           err=-2
        endif
     end do

     !Reading word
     do while(err==0)
        i=i+1
        read(10,'(a1)',advance='no',iostat=err)word(i:i)

        if(err==-2)then
           next_new=.true.
        elseif(word(i:i)=='#')then
           read(10,*)
           word(i:i)=' '
           next_new=.true.
           err=-2
        elseif((word(i:i)==' ').or.(word(i:i)==achar(9)))then
           word(i:i)=' '
           exit
        endif

     end do
  endif

end subroutine read_word





!======================================================
!Fill by .true. all the blocks of an array() associated 
!======================================================
subroutine symb2tab(natoms,tab_at,atom_type)
  implicit none

  !From calling modulus
  integer natoms
  logical, dimension(natoms) :: tab_at
  character(len=3), dimension(natoms) :: atom_type

  !Local
  integer:: i,j,end_line
  character(len=100) :: chain

  !Initialization
  i=1
  j=0
  end_line=0
  chain=''
  

  write(*,*)'Format to use: Symb1,Symb2,...' 
  do while (end_line==0)
     read(*,'(a)',advance='no',iostat=end_line)chain(i:i)!The line is read character by character to interpret ","
     if(chain(i:i)==",")then
        do j=1,natoms
           if(atom_type(j).eq.chain(1:i-1))then
              tab_at(j)=.true.
           endif
        enddo
        i=0
     endif
     i=i+1
  enddo
  do j=1,natoms
     if(atom_type(j).eq.chain(1:i-1))then
        tab_at(j)=.true.
     endif
  enddo

end subroutine symb2tab








!===============================
!Centering data according to pbc
!===============================
subroutine recenter(vel,c)
  use mymod, only : velocity
  implicit none

  !From upper subroutine
  double precision c
  type (velocity) :: vel

  !All the boundaries are recenter
  !This operation is not done so often (max twice per vvaf)
  vel%min1_r=vel%min1_r-dble(floor(vel%min1_r/c))*c
  vel%max1_r=vel%max1_r-dble(floor(vel%max1_r/c))*c
  vel%min2_r=vel%min2_r-dble(floor(vel%min2_r/c))*c
  vel%max2_r=vel%max2_r-dble(floor(vel%max2_r/c))*c
  vel%min1=vel%min1-dble(floor(vel%min1/c))*c
  vel%max1=vel%max1-dble(floor(vel%max1/c))*c
  vel%min2=vel%min2-dble(floor(vel%min2/c))*c
  vel%max2=vel%max2-dble(floor(vel%max2/c))*c

  !If min=max, it means we do not want to take into account any atom
  !So min is in [0;c[ and max also
  !But if min/=max, it is important to have
  !min in [0;c[, and max in ]0;c])
  if((abs(vel%min1_r)>epsilon(0d0)).and.(abs(vel%max1_r)<=epsilon(0d0)))vel%max1_r=c
  if((abs(vel%min2_r)>epsilon(0d0)).and.(abs(vel%max2_r)<=epsilon(0d0)))vel%max2_r=c
  if((abs(vel%min1)>epsilon(0d0)).and.(abs(vel%max1)<=epsilon(0d0)))vel%max1=c
  if((abs(vel%min2)>epsilon(0d0)).and.(abs(vel%max2)<=epsilon(0d0)))vel%max2=c


end subroutine recenter


!==========================================
!Update the layers if there is an increment
!==========================================
subroutine update_layer(opt_vel,end_vvaf_block,inc_val,a,b,c)
  use mymod, only : velocity,recenter
  implicit none

  !From upper subroutine
  double precision inc_val,a,b,c
  logical end_vvaf_block
  type (velocity) :: opt_vel

  !===========================
  !M1-m2 increment will change
  !===========================
  if(opt_vel%inc_param)then

     !Incremented boundaries
     if(&
          (abs(opt_vel%min1_r)<=epsilon(0d0)).and.&
          (abs(opt_vel%max1_r)<=epsilon(0d0)))then!SLAB1 is untouched if it is not defined
        opt_vel%min2=opt_vel%min2-inc_val
        opt_vel%max2=opt_vel%max2-inc_val
     elseif(&
          (abs(opt_vel%min2_r)<=epsilon(0d0)).and.&
          (abs(opt_vel%max2_r)<=epsilon(0d0)))then!SLAB2 is untouched if it is not defined
        opt_vel%min1=opt_vel%min1+inc_val
        opt_vel%max1=opt_vel%max1+inc_val
     else!SLAB1 and SLAB2 are updated because they are both defined
        opt_vel%min1=opt_vel%min1+inc_val
        opt_vel%max1=opt_vel%max1+inc_val
        opt_vel%min2=opt_vel%min2-inc_val
        opt_vel%max2=opt_vel%max2-inc_val
     endif

     !Recenter
     if(trim(opt_vel%vol_meth)=='X')then
        call recenter(opt_vel,a)
     elseif(trim(opt_vel%vol_meth)=='Y')then
        call recenter(opt_vel,b)
     elseif(trim(opt_vel%vol_meth)=='Z')then
        call recenter(opt_vel,c)
     endif

     !Condition to go out of the loop
     !2 times 2 conditions have to be studied
     !1) SLAB1 or SLAB2 is not defined (while the other is)
     !2) min<max or man<min (possible because of pbc)
     if(&
          !SLAB1 is defined
          ((abs(opt_vel%min1_r)>epsilon(0d0)).and.&
          (abs(opt_vel%max1_r)>epsilon(0d0)).and.&
          (opt_vel%max1>=opt_vel%max1_r).and.(&
          !min1<max1
          (opt_vel%min1_r<=opt_vel%max1_r).or.&!Normal case when min<=max
          !max1>min1 
          ((opt_vel%min1_r>opt_vel%max1_r).and.(opt_vel%max1<opt_vel%min1_r))&
          )).or.&
          !SLAB2 is defined
          ((abs(opt_vel%min2_r)>epsilon(0d0)).and.&
          (abs(opt_vel%max2_r)>epsilon(0d0)).and.&
          (opt_vel%min2<=opt_vel%min2_r).and.(&
          !min2<max2
          (opt_vel%min2_r<=opt_vel%max2_r).or.&!Normal case when min<=max
          !max2>min2 
          ((opt_vel%min2_r>opt_vel%max2_r).and.(opt_vel%min2>opt_vel%max2_r))&
          ))&
          )then
        opt_vel%min2=opt_vel%min2_r!Just in case the max-min is not a multiple of the increment 
        opt_vel%max1=opt_vel%max1_r
        end_vvaf_block=.true.
     else
        end_vvaf_block=.false.
     endif
     
  !===========================
  !m1-M2 increment will change
  !===========================
  else
     !Incremented boundaries
     if(&
          (abs(opt_vel%min1_r)<=epsilon(0d0)).and.&
          (abs(opt_vel%max1_r)<=epsilon(0d0)))then!SLAB1 is untouched if it is not defined
        opt_vel%min2=opt_vel%min2+inc_val
        opt_vel%max2=opt_vel%max2+inc_val
     elseif(&
          (abs(opt_vel%min2_r)<=epsilon(0d0)).and.&
          (abs(opt_vel%max2_r)<=epsilon(0d0)))then!SLAB2 is untouched if it is not defined
        opt_vel%min1=opt_vel%min1-inc_val
        opt_vel%max1=opt_vel%max1-inc_val
     else!SLAB1 and SLAB2 are updated because they are both defined
        opt_vel%min1=opt_vel%min1-inc_val
        opt_vel%max1=opt_vel%max1-inc_val
        opt_vel%min2=opt_vel%min2+inc_val
        opt_vel%max2=opt_vel%max2+inc_val
     endif

     !Recenter
     if(trim(opt_vel%vol_meth)=='X')then
        call recenter(opt_vel,a)
     elseif(trim(opt_vel%vol_meth)=='Y')then
        call recenter(opt_vel,b)
     elseif(trim(opt_vel%vol_meth)=='Z')then
        call recenter(opt_vel,c)
     endif

     !Condition to go out of the loop
     !2 times 2 conditions have to be studied
     !1) SLAB1 or SLAB2 is not defined (while the other is)
     !2) min<max or man<min (possible because of pbc)
     if(&
          !SLAB1 is defined
          ((abs(opt_vel%min1_r)>epsilon(0d0)).and.&
          (abs(opt_vel%max1_r)>epsilon(0d0)).and.&
          (opt_vel%min1<=opt_vel%min1_r).and.(&
          !min1<max1
          (opt_vel%min1_r<=opt_vel%max1_r).or.&!Normal case when min<=max
          !max1<min1 
          ((opt_vel%min1_r>opt_vel%max1_r).and.(opt_vel%min1>opt_vel%max1_r))&
          )).or.&
          !SLAB2 is defined
          ((abs(opt_vel%min2_r)>epsilon(0d0)).and.&
          (abs(opt_vel%max2_r)>epsilon(0d0)).and.&
          (opt_vel%max2>=opt_vel%max2_r).and.(&
          !min2<max2
          (opt_vel%min2_r<=opt_vel%max2_r).or.&!Normal case when min<=max
          !max2<min2 
          ((opt_vel%min2_r>opt_vel%max2_r).and.(opt_vel%max2<opt_vel%min2_r))&
          ))&
          )then
        opt_vel%min1=opt_vel%min1_r
        opt_vel%max2=opt_vel%max2_r
        end_vvaf_block=.true.
     else
        end_vvaf_block=.false.
     endif
  endif

end subroutine update_layer



!===========================================
!Subroutine to modify the name of the output
!if one has used the incremental layering
!===========================================
subroutine update_filename(name,name_tmp,v0,vt)
  use mymod, only : velocity
  implicit none

  !From upper subroutine
  character(len=100) name,name_tmp
  type (velocity) :: v0,vt

  !Local
  integer i,j,vel
  character(len=100) min1_0,min2_0,max1_0,max2_0
  character(len=100) min1_t,min2_t,max1_t,max2_t

  !Initialization
  i=0 ; vel=0
  min1_0='' ; min2_0='' ; max1_0='' ; max2_0=''
  min1_t='' ; min2_t='' ; max1_t='' ; max2_t=''

  !================================
  !Record of the min and max for v0
  !================================
  !1)Conversion of the min/max from rounded double (real) to character 
  if(v0%inc_param)then
     write(min1_0,*)real(v0%min1_r) ; min1_0=adjustl(min1_0)
     write(max1_0,*)real(v0%max1)   ; max1_0=adjustl(max1_0)
     write(min2_0,*)real(v0%min2)   ; min2_0=adjustl(min2_0)
     write(max2_0,*)real(v0%max2_r) ; max2_0=adjustl(max2_0)
  else
     write(min1_0,*)real(v0%min1)   ; min1_0=adjustl(min1_0)
     write(max1_0,*)real(v0%max1_r) ; max1_0=adjustl(max1_0)
     write(min2_0,*)real(v0%min2_r) ; min2_0=adjustl(min2_0)
     write(max2_0,*)real(v0%max2)   ; max2_0=adjustl(max2_0)
  endif

  !2)Remove the last 0
  i=len_trim(min1_0)
  do while(min1_0(i:i)=='0')
     min1_0(i:i)=''
     i=i-1
  end do
  i=len_trim(max1_0)
  do while(max1_0(i:i)=='0')
     max1_0(i:i)=''
     i=i-1
  end do
  i=len_trim(min2_0)
  do while(min2_0(i:i)=='0')
     min2_0(i:i)=''
     i=i-1
  end do
  i=len_trim(max2_0)
  do while(max2_0(i:i)=='0')
     max2_0(i:i)=''
     i=i-1
  end do

  !===========
  !Same for vt
  !===========
  !1bis)Conversion of the min/max from rounded double (real) to character 
  if(vt%inc_param)then
     write(min1_t,*)real(vt%min1_r) ; min1_t=adjustl(min1_t)
     write(max1_t,*)real(vt%max1)   ; max1_t=adjustl(max1_t)
     write(min2_t,*)real(vt%min2)   ; min2_t=adjustl(min2_t)
     write(max2_t,*)real(vt%max2_r) ; max2_t=adjustl(max2_t)
  else
     write(min1_t,*)real(vt%min1)   ; min1_t=adjustl(min1_t)
     write(max1_t,*)real(vt%max1_r) ; max1_t=adjustl(max1_t)
     write(min2_t,*)real(vt%min2_r) ; min2_t=adjustl(min2_t)
     write(max2_t,*)real(vt%max2)   ; max2_t=adjustl(max2_t)
  endif

  !2bis)Remove the last 0
  i=len_trim(min1_t)
  do while(min1_t(i:i)=='0')
     min1_t(i:i)=''
     i=i-1
  end do
  i=len_trim(max1_t)
  do while(max1_t(i:i)=='0')
     max1_t(i:i)=''
     i=i-1
  end do
  i=len_trim(min2_t)
  do while(min2_t(i:i)=='0')
     min2_t(i:i)=''
     i=i-1
  end do
  i=len_trim(max2_t)
  do while(max2_t(i:i)=='0')
     max2_t(i:i)=''
     i=i-1
  end do


  !==========================================
  !3)Fill step by step name_tmp
  !The '%v0m1', '%vtm1',...  in the title are
  !replaced by their true values
  !==========================================
  i=1!increment for name
  j=1!j = increment for name_tmp
  name_tmp=''
  do while (i<=len_trim(name)-4)!The size of the variables is 1+4 (%v0m1,...)

     if((name(i:i)/='%'))then!If it is a normal character
        name_tmp(j:j)=name(i:i)
        i=i+1
        j=j+1
     else
        select case(name(i:i+4))
        case('%v0m1')
           name_tmp(j:)=trim(min1_0)
           i=i+5
           j=j+len_trim(min1_0)

        case('%v0M1')
           name_tmp(j:)=trim(max1_0)
           i=i+5
           j=j+len_trim(max1_0)

        case('%v0m2')
           name_tmp(j:)=trim(min2_0)
           i=i+5
           j=j+len_trim(min2_0)

        case('%v0M2')
           name_tmp(j:)=trim(max2_0)
           i=i+5
           j=j+len_trim(max2_0)

        case('%vtm1')
           name_tmp(j:)=trim(min1_t)
           i=i+5
           j=j+len_trim(min1_t)

        case('%vtM1')
           name_tmp(j:)=trim(max1_t)
           i=i+5
           j=j+len_trim(max1_t)

        case('%vtm2')
           name_tmp(j:)=trim(min2_t)
           i=i+5
           j=j+len_trim(min2_t)

        case('%vtM2')
           name_tmp(j:)=trim(max2_t)
           i=i+5
           j=j+len_trim(max2_t)

        case default!If '%' is not followed by the good letters it is considered as a normal character
           name_tmp(j:j)=name(i:i)
           i=i+1
           j=j+1
        end select
     endif
     
  enddo
  
  !Last numbers
  do while (i<=len_trim(name))
     name_tmp(j:j)=name(i:i)
     i=i+1
     j=j+1
  end do

end subroutine update_filename

