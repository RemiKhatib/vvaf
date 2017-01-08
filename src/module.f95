!=============================================================================================
!Module
!=============================================================================================
module mymod
  implicit none

  double precision, parameter :: error = 1d-6
  
  !Definition of a structure to store the informations about common elements (already used atoms) and new elements (to define when the program runs)
  type :: info_atoms
     character(len=3) :: name!Symbol
     double precision :: mw!Weight
     integer :: amount!Number of atoms "name"
  end type info_atoms

  !Information about a specific bond
  !Used to store efficiently the main information about an atom which makes a bond
  type :: info_bond
     integer :: index !Index of the bonded atom. 
     integer :: third !Index of the third atom which will alow to determine the molecule framework
     double precision, dimension(3) :: cos !Cosine of the angle between the bond and the X,Y,Z axis
     double precision :: dist !Distance square between the studied atom and the bonded atom 
  end type info_bond
  
  !Information about the velocities treatment for v0 or vt
  type :: velocity
     integer :: vproj!For projection of the velocity on an xyz axis (scalar not vector)
     logical :: scal_proj!For projection of the velocity on the bond axis (scalar and not vector)
     logical :: bond_diff!false if atomic velocities are studied, true if bond velocities
     logical, dimension(:), pointer :: tab_select!if tab_select(i) true, the velocity of this atom will be studied
     integer :: P,Q,R!Polarization of the 3 different beams (0=nothing, 1=X, 2=Y, 3=Z)
     double precision, dimension(3) :: dip_stretch!dipole derivative vector (stretching)
     double precision :: dip_bend, dip_ss, dip_as!dipole derivative vector (c2v symmetry: A1 RI only the component along z, B1 RI only along x)
     double precision, dimension(3,3) :: pol_stretch!polarizability derivative matrix (stretching, the matrix is symmetric)
     double precision, dimension(3) :: pol_bend, pol_ss!polarizability derivative matrix (bending and SS, only the diagonal)
     double precision :: pol_as!polarizability derivative matrix (AS, only XZ)
     double precision, dimension(3,3,3) :: molf2q!Normal modes in the moleculare framework
     logical, dimension(:), pointer :: tab_at_bond!Record the atoms liable to be bonded to the atoms whom the velocity is studied
     integer :: nb_bond!Number of bond that the central atom can do
     type (info_bond), dimension(:,:), pointer :: tab_mol! tab_mol(:,i)=j means that atom j is bonded to the atom i
     double precision :: dmin2,dmax2 !Square of the minimal and the maximal distances of the bonds
     character :: vol_meth!method to select the interesting volume (depends on interface distance or atom distance) 
     logical, dimension(:), pointer :: tab_center!It is the atoms which serve as center if vol_meth="r". In this case only an annulus with min1<=r<=max1 is selected.
     logical :: v0_at!When true the center of the ANNULUS at t is the atom studied at t0  (works only with vol_meth="r")
     double precision :: min1_r,max1_r,min2_r,max2_r!Reference.Take into account only the atoms in the good layer(s). The orientation depends vol_meth
     double precision :: min1,max1,min2,max2!Idem than previously, but they are the values really used if one of this parameter has to be update with an increment.
     logical :: inc_param!Increment for the layering, inc_val is associated with the value (min or max) which will be incremented
     logical :: stretch_mode!Stretching prefactor (vector not scalar. Different from scal_proj)
     integer :: c2v_mode!For molecule with water symmetry (c2v): 0 nothing (bond), 1 simple c2v mode, 2 hessian c2v mode
     logical :: c2v_deriv!check if we have defined c2v derivatives (if false, it is stretching derivatives or nothing)
  end type velocity

  !Record all the options associated with the vvaf number X (name,selected atoms,)
  type :: vvaf_options
     character(len=100) :: filename,filename_tmp
     double precision :: inc_val!Increment value for the layering
     logical :: inc,cross!Increment used or not for the layering - General case: cross-correlation
     type (velocity) :: v0,vt
  end type vvaf_options
  
  interface
     !====================================
     !Information about the weight of atom
     !(for mass center)
     !====================================
     subroutine weight_info(ptr_info,kindmax)
       import :: info_atoms
       integer, intent(inout) :: kindmax
       type (info_atoms), dimension(:), pointer, intent(inout) :: ptr_info 
     end subroutine weight_info

     !===================================================
     !Reading of the parameters from standard input
     !Question about the calculation which has to be done
     !===================================================
     subroutine read_input(nb_read,filetr,filepos,nmovie,natoms,istep,kindmax,Di,shift,&
          ptr_info,atom_type,fmass,a,b,c,box_size,center,surface,slab_orient,step_min,step_max,&
          vx,vy,vz,x,y,z)
       import :: info_atoms,vvaf_options
       integer, intent(inout) :: nb_read
       integer, intent(inout) :: nmovie,natoms,istep,kindmax,Di,shift,step_min,step_max
       double precision, intent(inout) ::  a,b,c,box_size,center,surface
       double precision, dimension(:), allocatable, intent(inout) :: fmass
       double precision,dimension(:,:), intent(in):: vx,vy,vz,x,y,z
       character, intent(inout) ::  slab_orient
       character(len=100), intent(inout) :: filetr,filepos
       character(len=3), dimension(:), allocatable, intent(inout) :: atom_type
       type (info_atoms), dimension(:), pointer, intent(inout) :: ptr_info 
     end subroutine read_input
     subroutine input_main(nb_read,nmovie,natoms,istep,kindmax,Di,shift,step_min,step_max,a,b,c,box_size,center,surface,&
          fmass,filetr,filepos,slab_orient,atom_type,ptr_info)
       import :: vvaf_options,info_atoms
       integer, intent(in):: nb_read
       integer, intent(inout) :: nmovie,natoms,istep,kindmax,Di,shift,step_min,step_max
       double precision, intent(inout) :: a,b,c,box_size,center,surface
       double precision, dimension(:), allocatable, intent(inout) :: fmass
       character, intent(inout) ::  slab_orient
       character(len=100), intent(inout):: filetr,filepos
       character(len=3), dimension(:), allocatable, intent(inout) :: atom_type
       type (info_atoms), dimension(:), pointer, intent(inout) :: ptr_info 
     end subroutine input_main
     subroutine input_main_atweight(kindmax,ptr_info)
       import::info_atoms
       integer, intent(inout)::kindmax
       type (info_atoms),dimension(:), pointer,intent(inout) :: ptr_info
     end subroutine input_main_atweight
     subroutine input_vvaf(nb_read,filepos,nmovie,natoms,options,atom_type,fmass,a,b,c,box_size,center,shift,Di,slab_orient,&
          surface,vx,vy,vz,x,y,z)
       import :: vvaf_options
       integer, intent(in) :: nb_read,nmovie,natoms,Di,shift
       double precision,intent(in) :: a,b,c,box_size,center,surface
       double precision,dimension(:,:),intent(in):: vx,vy,vz,x,y,z
       double precision, dimension(natoms), intent(in) :: fmass
       character, intent(in) :: slab_orient
       character(len=3), dimension(natoms), intent(in) :: atom_type
       character(len=100), intent(in) :: filepos
       type (vvaf_options), intent(inout) :: options
     end subroutine input_vvaf
     subroutine input_vvaf_velocity(nb_vel,filepos,natoms,vel,inc,cross,atom_type,fmass,a,b,c)
       import :: velocity
       integer, intent(in) :: natoms,nb_vel
       double precision, intent(in) :: a,b,c
       double precision, dimension(natoms), intent(in) :: fmass
       character(len=3), dimension(natoms), intent(in) :: atom_type
       character(len=100), intent(in) :: filepos
       logical, intent(inout) :: inc,cross
       type (velocity), intent(inout) :: vel
     end subroutine input_vvaf_velocity
     subroutine read_tab_dipol(vel,next_new)
       import :: velocity
       logical, intent(inout) :: next_new
       type (velocity), intent(inout) :: vel
     end subroutine read_tab_dipol
     subroutine skip_line_tab_dipol(next_new)
       logical, intent(inout) :: next_new
     end subroutine skip_line_tab_dipol
     subroutine read_tab_pol(vel,next_new)
       import :: velocity
       logical, intent(inout) :: next_new
       type (velocity), intent(inout) :: vel
     end subroutine read_tab_pol
     subroutine skip_line_tab_pol(next_new)
       logical, intent(inout) :: next_new
     end subroutine skip_line_tab_pol
     subroutine error_deriv_dip_pol()
     end subroutine error_deriv_dip_pol
     subroutine input_vvaf_velocity_atom(natoms,vel,atom_type,a,b,c,inc)
       import :: velocity
       integer, intent(in):: natoms
       double precision, intent(in) :: a,b,c
       character(len=3), dimension(natoms), intent(in) :: atom_type
       logical, intent(inout) :: inc
       type (velocity), intent(inout) :: vel
     end subroutine input_vvaf_velocity_atom
     subroutine input_vvaf_velocity_atom_bondw(natoms,vel,atom_type)
       import :: velocity
       integer, intent(in):: natoms
       character(len=3), dimension(natoms), intent(in) :: atom_type
       type (velocity), intent(inout) :: vel
     end subroutine input_vvaf_velocity_atom_bondw
     subroutine input_vvaf_velocity_atom_layer(vel,a,b,c,inc)
       import :: velocity
       double precision, intent(in) :: a,b,c
       logical, intent(inout) :: inc
       type (velocity), intent(inout) :: vel
     end subroutine input_vvaf_velocity_atom_layer
     subroutine input_vvaf_velocity_atom_annulus(natoms,vel,atom_type)
       import :: velocity
       integer, intent(in) :: natoms
       character(len=3), dimension(natoms), intent(in) :: atom_type
       type (velocity), intent(inout) :: vel
     end subroutine input_vvaf_velocity_atom_annulus
     subroutine input_vvaf_velocity_hessian(vel)
       import :: velocity
       type (velocity), intent(inout) :: vel
     end subroutine input_vvaf_velocity_hessian
     subroutine wrong_keyword(chain,block_name)
       character(len=100) block_name, chain
     end subroutine wrong_keyword
     subroutine read_word(word,err,next_new)
       character(len=100), intent(inout) :: word
       integer, intent(inout) :: err
       logical, intent(inout) :: next_new
     end subroutine read_word
     
     
     !======================
     !Getting information:
     !-natoms
     !-starting step (istep)
     !-filling fmass
     !-filling atom_type
     !======================
     subroutine read_step1(natoms,istep,kindmax,fmass,atom_type,ptr_info)
       import :: info_atoms
       integer, intent(inout) :: natoms,istep,kindmax
       double precision, dimension(:), allocatable, intent(inout) :: fmass  
       character(len=3), dimension(:), allocatable, intent(inout) :: atom_type
       type (info_atoms), dimension(:), pointer, intent(inout) :: ptr_info 
     end subroutine read_step1



     !===========================================
     !Tables allocation and velocity file reading
     !===========================================
     subroutine reading(filepos,natoms,nmovie,istep,fmass,a,b,c,step_min,step_max,&
          vx,vy,vz,x,y,z)
       integer, intent(in) :: natoms,nmovie,istep,step_min,step_max
       double precision, intent(in) :: a,b,c
       double precision, dimension(:), allocatable, intent(inout) :: fmass  
       double precision, dimension(natoms,nmovie) :: vx,vy,vz,x,y,z 
       character(len=100), intent(in) :: filepos

     end subroutine reading

     !===============================
     !Centering data according to pbc
     !===============================
     subroutine recenter(vel,c)
       import :: velocity
       double precision, intent(in) :: c
       type (velocity), intent(inout) :: vel
     end subroutine recenter

     !==========================================
     !Update the layers if there is an increment
     !==========================================
     subroutine update_layer(opt_vel,end_vvaf_block,inc_val,a,b,c)
       import :: velocity
       double precision, intent(in):: inc_val,a,b,c
       logical,intent(inout) :: end_vvaf_block
       type (velocity),intent(inout) :: opt_vel
     end subroutine update_layer

     !===========================================
     !Subroutine to modify the name of the output
     !if one has used the incremental layering
     !===========================================
     subroutine update_filename(name,name_tmp,v0,vt)
       import :: velocity
       character(len=100),intent(inout):: name,name_tmp
       type (velocity),intent(in) :: v0,vt
     end subroutine update_filename


     !================
     !VVAF calculation
     !================
     subroutine vvaf(nmovie,natoms,Di,shift,options,vx,vy,vz,x,y,z,&
          slab_orient,center,a,b,c,box_size,surface,corr)
       import :: vvaf_options
       integer, intent(in) :: nmovie,natoms,Di,shift
       double precision,intent(in) :: vx,vy,vz,x,y,z
       double precision,intent(in) :: a,b,c,box_size,center,surface
       double precision,dimension(-Di:Di)::corr
       character, intent(in) :: slab_orient
       type (vvaf_options), intent(inout) :: options
     end subroutine vvaf

     !=========================================================================
     !Modification of the velocity of the atom according to different criteria
     !All the prefactors are calculated for the studied velocity (ready to use)
     !=========================================================================
     subroutine bond_calc(natoms,i,j,box_size,center,slab_orient,opt_vel,vout,&
          vx,vy,vz,x,y,z)
       import :: velocity
       integer, intent(in) :: natoms,i,j
       double precision, intent(in) :: box_size,center
       double precision, dimension(3), intent(inout) ::vout
       double precision, dimension(natoms),intent(in)::vx,vy,vz,x,y,z
       character, intent(in):: slab_orient
       type (velocity),intent(in) :: opt_vel
     end subroutine bond_calc

     !======================================================================
     !Calculation of the time derivatives for molecule having a c2v symmetry
     !but with a simple description of the different modes
     !======================================================================
     subroutine simple_c2v_calc(natoms,i,box_size,center,slab_orient,opt_vel,vout,vx,vy,vz,x,y,z)
       import :: velocity
       integer, intent(in) :: natoms,i
       double precision, intent(in) :: box_size,center
       double precision, dimension(3), intent(inout) ::vout
       double precision, dimension(natoms),intent(in)::vx,vy,vz,x,y,z
       character, intent(in):: slab_orient
       type (velocity),intent(in) :: opt_vel
     end subroutine simple_c2v_calc

     !======================================================================
     !Calculation of the time derivatives for molecule having a c2v symmetry
     !but with a hessian description of the different modes
     !======================================================================
     subroutine hessian_c2v_calc(natoms,i,box_size,center,slab_orient,opt_vel,vout,vx,vy,vz,x,y,z)
       import :: velocity
       integer, intent(in) :: natoms,i
       double precision, intent(in) :: box_size,center
       double precision, dimension(3), intent(inout) ::vout
       double precision, dimension(natoms),intent(in)::vx,vy,vz,x,y,z
       character, intent(in):: slab_orient
       type (velocity),intent(in) :: opt_vel
     end subroutine hessian_c2v_calc

     !=============================================================
     !Function which checks if an atom is inside the desired volume
     !=============================================================
     logical function f_in_vol(natoms,i0,a,b,c,x,y,z,tab_x,tab_y,tab_z,opt_vel)
       import :: velocity
       integer, intent(in) :: natoms,i0
       double precision, intent(in) :: a,b,c,x,y,z
       double precision, dimension(natoms), intent(in):: tab_x,tab_y,tab_z
       type (velocity), intent(in) :: opt_vel
     end function f_in_vol


     !===================================
     !Calculation of the norm of a vector
     !===================================
     subroutine vect2proj(v,bond)
       import :: info_bond
       double precision, dimension(3), intent(inout):: v
       type (info_bond), intent(in) :: bond
     end subroutine vect2proj


     !=============================================
     !Projection of the velocity on a specific axis
     !=============================================
     subroutine vel_proj(axis,vout)
       integer,intent(in):: axis
       double precision, dimension(3), intent(inout):: vout
     end subroutine vel_proj

     !==================================================================================================
     !Subroutine which determines the atoms associated to the same molecule
     !The bonded atoms are selected among those of options%tab_at_bond (except the iatom itself)
     !==================================================================================================
     subroutine make_molecule(natoms,x,y,z,a,b,c,vel,tab_new_mol,record)
       import :: velocity
       integer, intent(in) :: natoms
       double precision, intent(in) :: a,b,c
       double precision, dimension(natoms), intent(in) :: x,y,z
       logical, intent(in) :: record
       logical, dimension(natoms), intent(inout) :: tab_new_mol
       type (velocity), intent(inout) :: vel
     end subroutine make_molecule


     !===============================================================
     !Subroutine which determines which atom k will be used to define
     !xyz (the molecule framework)
     !===============================================================
     subroutine third_at(natoms,opt_vel,i,j,x,y,z,a,b,c)
       import :: velocity
       integer, intent(in) :: natoms,i,j
       double precision, intent(in):: a,b,c
       double precision, dimension(natoms), intent(in) :: x,y,z
       type (velocity), intent(inout) :: opt_vel
     end subroutine third_at


     !============================================================================================
     !Calculation of the prefactor associated with the dipole moment along R (if PQR poalrization)
     !============================================================================================
     double precision function f_dipole(c,center,at_center,ligand,third,opt_vel,z)
       import :: velocity
       integer, intent(in) :: at_center,ligand,third
       double precision, intent(in) :: c,center,z
       type (velocity), intent(in) :: opt_vel
     end function f_dipole


     !========================================================================================
     !Calculation of the prefactor associated with the PQ-polarizability (if PQR poalrization)
     !========================================================================================
     double precision function f_polariz(at_center,ligand,third,opt_vel)
       import :: velocity
       integer, intent(in) :: at_center,ligand,third
       type (velocity), intent(in) :: opt_vel
     end function f_polariz


     !======================================================
     !Subroutine which dertermines the molecule framework
     !which corresponds with the direction cosine matrix
     !The inversion center of the cell is taken into account
     !xyz is the bond framework projected on XYZ
     !the lab framework
     !======================================================
     subroutine bond_frame(cosx,cosy,cosz,opt_vel,at_center,ligand,third)
       import :: velocity
       integer, intent(in):: at_center,ligand,third
       double precision, dimension(3), intent(inout) :: cosx,cosy,cosz
       type (velocity), intent(in) :: opt_vel
     end subroutine bond_frame

     !======================================================
     !Subroutine which dertermines the molecular framework
     !which corresponds with the direction cosine matrix
     !The inversion center of the cell is taken into account
     !xyz is the molecular framework projected on XYZ
     !the lab framework
     !======================================================
     subroutine mol_frame(cosb1,cosb2,n1,n2,cosm)
       double precision, dimension(3), intent(inout) :: cosb1, cosb2, n1, n2
       double precision, dimension(3,3), intent(inout) :: cosm 
     end subroutine mol_frame
     
     !======================================================================
     !Subroutine which calculates the projection of the velocity on the bond
     !======================================================================
     subroutine stretch(vout,bond)
       import :: info_bond
       double precision, dimension(3),intent(inout) :: vout
       type (info_bond), intent(in) :: bond
     end subroutine stretch

     subroutine bond_update(dx,dy,dz,bond,a,b,c)
       import :: info_bond
       double precision, intent(in) ::  a,b,c,dx,dy,dz
       type (info_bond), intent(inout)  :: bond
     end subroutine bond_update

     !================================
     !Calculation of the bond velocity
     !================================
     subroutine vel_diff(opt_vel,vatx,vaty,vatz,vrefx,vrefy,vrefz,vout)
       import :: velocity
       double precision, intent(in) :: vatx,vaty,vatz,vrefx,vrefy,vrefz
       double precision, dimension(3), intent(inout) :: vout
       type (velocity), intent(in) :: opt_vel
     end subroutine vel_diff



     !=================================================================================
     !Subroutine  which does some velocity treatments associated with bonds
     !-Update of the bonds at t (the molecules are defined at each t0)
     !-Vecocity difference between two atoms
     !-Calculation of the dipole moment prefactor (projectiopn of the bond on the axis)
     !-Calculation of polarizability (with our level of theory: only 0 or 1)
     !=================================================================================
     subroutine bond_only(natoms,opt_vel,i,ligand,vout,vx,vy,vz,x,y,z,box_size,center,slab_orient)
       import :: velocity
       integer, intent(in):: natoms,i,ligand
       double precision, intent(in):: box_size,center
       double precision, dimension(3), intent(inout):: vout
       double precision, dimension(natoms), intent(in)::vx,vy,vz,x,y,z
       type (velocity), intent(in) :: opt_vel
       character, intent(in):: slab_orient
     end subroutine bond_only

     !=============
     !Cross product
     !=============
     subroutine cross_prod(v1,v2,v3)
       double precision, dimension(3), intent(in)::v1,v2
       double precision, dimension(3), intent(inout)::v3
     end subroutine cross_prod


     !=============================================
     !Multiplication of dq/dt by the coefficients
     !associated with the transition dipole moments
     !and polarizabilities
     !=============================================
     subroutine velq2dippol(opt_vel,vout,box_size,center,x,y,z,slab_orient,cosm)
       import :: velocity
       type (velocity), intent(in) :: opt_vel
       double precision, intent(in):: x,y,z
       double precision, dimension(3), intent(inout)::vout 
       double precision, dimension(3,3), intent(in) :: cosm 
       double precision, intent(in) :: box_size,center
       character, intent(in):: slab_orient
     end subroutine velq2dippol



     
     !==============================================================
     !This subroutine is especially done in order to multiply
     !a rotation matrix by a vector having only one non-nul element
     !and extracting a single element of the resulting vector
     !
     !E.g. (Rotation matrix) x (transition dipole moment of bending along z)
     !projected on R (from PQR polarization)
     !==============================================================
     double precision function tdip_proj(cosm,dip,R)
       integer ::R !Polarization of the IR beam (X=1, Y=2, Z=3, nothing=0)
       double precision dip !transition dipole moment
       double precision, dimension(3) :: cosm !Column of the direction cosine matrix
     end function tdip_proj

     !==============================================================
     !This subroutine is especially done in order to multiply
     !a rotation matrix by a sparse tensor (from a normal mode)
     !and extracting a single element of the resulting vector
     !
     !E.g. Element PQ (from PQR polarization)
     !(Rotation matrix) x (transition polarizability of bending along z)
     !==============================================================
     double precision function tpol_proj_A1(cosm,pol,P,Q)
       integer ::P !First element of the polarizability tensor (X=1, Y=2, Z=3, nothing=0)
       integer ::Q !Second element of the polarizability tensor (X=1, Y=2, Z=3, nothing=0)
       double precision, dimension(3) :: pol ! transition polarizability (only three values are requiered)
       double precision, dimension(3,3) :: cosm ! direction cosine matrix
     end function tpol_proj_A1


     !==============================================================
     !The same for B1 symmetry (alpha_XZ=alpha_ZX, the rest is 0)
     !alpha_PQ = alpha_xz (Pz xQ + Px zQ)
     !==============================================================
     double precision function tpol_proj_B1(cosm,pol,P,Q)
       integer ::P !First element of the polarizability tensor (X=1, Y=2, Z=3, nothing=0)
       integer ::Q !Second element of the polarizability tensor (X=1, Y=2, Z=3, nothing=0)
       double precision :: pol ! transition polarizability (only three values are requiered)
       double precision, dimension(3,3) :: cosm ! direction cosine matrix
     end function tpol_proj_B1

  end interface

end module mymod
