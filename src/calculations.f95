!This file contains the subroutine doing the calculations of X(t).Y(0)
!-bond_calc for the bond version
!-simple_c2v_calc for the version with "pure/brute/simple" SS, AS and bending
!-hessian_c2v_calc for the version using the modes defined by the Hessian

!=======================================================================
!Subroutine which takes the velocity and the position of a specific atom
!to calculate its contribution (depends on opt_vel)
!Only vout is modified while vx, vy and vz are safe for next VVAF!!!
!=======================================================================
subroutine bond_calc(natoms,i,j,box_size,center,slab_orient,opt_vel,vout,vx,vy,vz,x,y,z)
  use mymod, only : velocity,vel_proj,vect2proj,f_in_vol,bond_only
  implicit none

  !From upper subroutine
  integer natoms,i,j
  double precision :: box_size,center
  double precision, dimension(3)::vout
  double precision, dimension(natoms)::vx,vy,vz,x,y,z
  character slab_orient
  type (velocity) :: opt_vel

  !=======================================================!
  !=========================REMARK========================!
  !At the end of the subroutine, the vout MUST be updtated!
  !If no update is requiered, vout(:)=0.d0                !
  !To do less calculation if "shift < 2*Di+1", the        !
  !previous values of vout are not deleted at the end of  !
  !the step. It can be problematic to propagate these     !
  !value if it is not requiered.                          !
  !=======================================================!


  !=======================================================================================
  !Treatment of the velocity which are pointless if no bond are involved (opt_vel%nb_bond)
  !=======================================================================================
  if(opt_vel%nb_bond>0)then!Bond defined
     !No treatment if there is not the good number of bonds (tab_mol(j,i)%index=0)
     if(opt_vel%tab_mol(j,i)%index>0)then
        call bond_only(natoms,opt_vel,i,j,vout,vx,vy,vz,x,y,z,&
             box_size,center,slab_orient)

        !==========================================================
        !Conversion from a vector to a scalar
        !The a projection of the velocity on the bond (as a scalar)
        ! cannot be done if there is also projection on an axis
        !The scalar is stored in the first row of the array
        !==========================================================
        if(opt_vel%scal_proj)then
           !The velocity is projected on the bond axis
           call vect2proj(vout(1),opt_vel%tab_mol(j,i))
        elseif(opt_vel%vproj/=0)then
           !The velocity is projected on a specific axis
           call vel_proj(opt_vel%vproj,vout(1))
        endif
        
     else!The bond disappears
        vout(:)=0.d0
     endif
  else!Atomic treatment
     !If no bond is requiered, only the atomic velocity is taken into account
     !The data are sent to vout which is the vector which will be send back to the upper subroutine
     vout(1)=vx(i)
     vout(2)=vy(i)
     vout(3)=vz(i)
          
     !The velocity is projected on a specific axis
     if(opt_vel%vproj/=0)call vel_proj(opt_vel%vproj,vout(1))
  endif

end subroutine bond_calc



!=======================================================================
!Calculation of the time derivatives for molecule having a c2v symmetry
!but with a simple description of the different modes
!vout(1)= bending along z (depends only on theta)
!vout(2)= AS along x (depends on r1-r2)
!vout(3)= SS along z (depends on r1+r2-2rO)
!=======================================================================
subroutine simple_c2v_calc(natoms,i,box_size,center,slab_orient,opt_vel,vout,vx,vy,vz,x,y,z)
  use mymod, only : velocity, cross_prod, vel_diff, vect2proj, mol_frame, velq2dippol
  implicit none

  !From upper subroutine
  integer natoms,i
  double precision :: box_size,center
  double precision, dimension(3)::vout 
  double precision, dimension(natoms)::vx,vy,vz,x,y,z
  character slab_orient
  type (velocity) :: opt_vel

  !Local variables
  double precision, dimension(3):: v1, v2
  !Direction cosine matrix of the molecule
  !z along the bissector
  !x in the molecular plane
  !y out of the molecular plane
  double precision, dimension(3,3) :: cosm 
  !Vectors normal to he bonds, in the molecular plane
  !and pointing in the same direction than z (z.n1>0 and z.n2>0)
  double precision, dimension(3):: n1, n2
  
  !Initialization
  v1(:)=0.d0 ; v2(:)=0.d0
  n1(:)=0.d0 ; n2(:)=0.d0 
  cosm(:,:)=0.d0

  !No treatment if there is not the good number of bonds (tab_mol(j,i)%index=0)
  !We only check the index of the first ligand, because if the first is good,
  !the second is necessarly good
  if(opt_vel%tab_mol(1,i)%index>0)then
     !=============================================================
     !At this point we are sure that there are one central atom (i)
     !and twice the same ligand
     !=============================================================
     !Molecular framework (direction cosine matrix)
     call mol_frame(opt_vel%tab_mol(1,i)%cos(1),opt_vel%tab_mol(2,i)%cos(1),&
          n1(1),n2(1),cosm(1,1))
     
     !Calculation of the velocity difference (bond velocity instead of atomic velocity)
     call vel_diff(opt_vel,&
          vx(opt_vel%tab_mol(1,i)%index),vy(opt_vel%tab_mol(1,i)%index),vz(opt_vel%tab_mol(1,i)%index),&
          vx(i),vy(i),vz(i),&
          v1(1))
     call vel_diff(opt_vel,&
          vx(opt_vel%tab_mol(2,i)%index),vy(opt_vel%tab_mol(2,i)%index),vz(opt_vel%tab_mol(2,i)%index),&
          vx(i),vy(i),vz(i),&
          v2(1))

     !===============================================================
     !Calculation of the angle speed
     !d theta / dt = (v1.n1)/l1 + (v2.n2)/l2
     !vi.ni is the projection of the velocity vector perpendicular
     !to the bond and in the molecular plane
     !Reminder: be careful to the sign, in both cases, n1 and n2 point
     !in the opposite direction of the first bissector.
     !Moreover, if we have - instead of + we have the librational mode
     !IMPORTANT REMARK: the angle unit is the RADIAN
     !===============================================================
     vout(1)=( &
          v1(1)*n1(1) + &
          v1(2)*n1(2) + &
          v1(3)*n1(3)   &
          )/opt_vel%tab_mol(1,i)%dist + &
          ( &
          v2(1)*n2(1) + &
          v2(2)*n2(2) + &
          v2(3)*n2(3)   &
          )/opt_vel%tab_mol(2,i)%dist

     !===========================================
     !Projection of the velocity on the bond axis
     !===========================================
     call vect2proj(v1(1),opt_vel%tab_mol(1,i))
     call vect2proj(v2(1),opt_vel%tab_mol(2,i))
     !Deduction of VAS and VSS (vout(2) and vout(3) respectively)
     vout(2)=v1(1)-v2(1)
     vout(3)=v1(1)+v2(1)

     !=============================================
     !Multiplication of dq/dt by the coefficients
     !associated with the transition dipole moments
     !and polarizabilities
     !Here, q is theta, r1-r2, r1+r2
     !=============================================
     call velq2dippol(opt_vel,vout,box_size,center,x(i),y(i),z(i),slab_orient,cosm(1,1))

  else!Central atoms which do not have the good number of ligands
     vout(:)=0.d0
  endif

end subroutine simple_c2v_calc




!=======================================================================
!Calculation of the time derivatives for molecule having a c2v symmetry
!but with a hessian description of the different modes
!vout(1)= bending along z 
!vout(2)= AS along x
!vout(3)= SS along z
!=======================================================================
subroutine hessian_c2v_calc(natoms,i,box_size,center,slab_orient,opt_vel,vout,vx,vy,vz,x,y,z)
  use mymod, only : velocity, mol_frame, velq2dippol
  implicit none

  !From upper subroutine
  integer natoms,i
  double precision :: box_size,center
  double precision, dimension(3)::vout 
  double precision, dimension(natoms)::vx,vy,vz,x,y,z
  character slab_orient
  type (velocity) :: opt_vel

  !Local variables
  integer j,k,l
  !Direction cosine matrix of the molecule
  !z along the bissector
  !x in the molecular plane
  !y out of the molecular plane
  double precision, dimension(3,3) :: cosm, labf2q
  !Vectors normal to he bonds, in the molecular plane
  !and pointing in the same direction than z (z.n1>0 and z.n2>0)
  double precision, dimension(3):: n1, n2
  
  !Initialization
  j=0 ; k=0 ; l=0
  cosm(:,:)=0.d0
  n1(:)=0.d0 ; n2(:)=0.d0

  !No treatment if there is not the good number of bonds (tab_mol(j,i)%index=0)
  !We only check the index of the first ligand, because if the first
  ! the second is necessarly good
  if(opt_vel%tab_mol(1,i)%index>0)then
     !=============================================================
     !At this point we are sure that there are one central atom (i)
     !and twice the same ligand
     !=============================================================
     !Molecular framework (direction cosine matrix)
     call mol_frame(opt_vel%tab_mol(1,i)%cos(1),opt_vel%tab_mol(2,i)%cos(1),&
          n1(1),n2(1),cosm(1,1))

     !============================================
     !Contribution of the modes in the molecule
     !At the end we obtain the mode speed (scalar)
     !dq/dt
     !============================================
     do j=1,3 !Loop for the modes
        !Rotation of molf2q in order to have the normal modes
        !expressed in the lab framework
        do k=1,3 !Loop for the atoms
           do l=1,3 !Loop for coordinates (lab frame)
              labf2q(l,k)=&
                   cosm(l,1)*opt_vel%molf2q(1,k,j) +&
                   cosm(l,2)*opt_vel%molf2q(2,k,j) +&
                   cosm(l,3)*opt_vel%molf2q(3,k,j)
           end do
        end do
     
        !Scalar product in order to determine the contribution of the mode
        !(mol velocities in lab frame).(labf2q)
        vout(j)=&
             vx(i)*labf2q(1,1)+&
             vy(i)*labf2q(2,1)+&
             vz(i)*labf2q(3,1)+&
             vx(opt_vel%tab_mol(1,i)%index)*labf2q(1,2)+&
             vy(opt_vel%tab_mol(1,i)%index)*labf2q(2,2)+&
             vz(opt_vel%tab_mol(1,i)%index)*labf2q(3,2)+&
             vx(opt_vel%tab_mol(2,i)%index)*labf2q(1,3)+&
             vy(opt_vel%tab_mol(2,i)%index)*labf2q(2,3)+&
             vz(opt_vel%tab_mol(2,i)%index)*labf2q(3,3)
        
     end do

     !=============================================
     !Multiplication of dq/dt by the coefficients
     !associated with the transition dipole moments
     !and polarizabilities
     !Here, q is the hessian B, AS and SS
     !=============================================
     call velq2dippol(opt_vel,vout,box_size,center,x(i),y(i),z(i),slab_orient,cosm(1,1))
     
  else!Central atoms which do not have the good number of ligands
     vout(:)=0.d0
  endif

end subroutine hessian_c2v_calc


  


!=============================================
!Multiplication of dq/dt by the coefficients
!associated with the transition dipole moments
!and polarizabilities
!=============================================
subroutine velq2dippol(opt_vel,vout,box_size,center,x,y,z,slab_orient,cosm)
  use mymod, only : velocity, tdip_proj, tpol_proj_A1, tpol_proj_B1
  implicit none

  !From upper subroutines
  type (velocity) :: opt_vel
  double precision, dimension(3):: vout 
  double precision, dimension(3,3) :: cosm 
  double precision :: box_size,center,x,y,z
  character :: slab_orient

  !Local variables
  double precision:: ptr_z

  !Initialization
  ptr_z=0.d0

  !Projection of the transition dipole moment in the mol frame
  !on the R axis (if PQR polarization)
  !Since we have normal modes:
  !-for bending and SS (A1) are along z (the prime bissector)
  !-for AS (B1) is along x (in the molecular plane)
  if(opt_vel%R/=0)then
     vout(1)=vout(1)*tdip_proj(cosm(1,3),opt_vel%dip_bend,opt_vel%R)
     vout(2)=vout(2)*tdip_proj(cosm(1,1),opt_vel%dip_as,opt_vel%R)
     vout(3)=vout(3)*tdip_proj(cosm(1,3),opt_vel%dip_ss,opt_vel%R)

     !===================================================================
     !If the molecule is on a certain side, we suppose that there is
     !an inversion center 
     !Remark: an inversion center has no impact on the the polarizability
     !So the equivalent does not exist for f_polariz
     !===================================================================
     if(slab_orient=='X')then
        ptr_z = x
     elseif(slab_orient=='Y')then
        ptr_z = y
     else
        ptr_z = z
     end if

     if(&
          ((center>box_size/2).and.(ptr_z<center).and.(ptr_z<=center-box_size/2)).or.&
          ((center<=box_size/2) .and.((ptr_z<center).or. (ptr_z>=center+box_size/2)))&
          )then
        vout(:)=-vout(:)
     endif

  endif

  !Projections of xyz (molecule framework) in PQ (axis of XYZ -- lab framework)
  !Since we have normal modes:
  !-for bending and SS (A1), we have a diagonal matrix
  !-for AS (B1), we have only XZ=ZX/=0
  if(opt_vel%P/=0)then
     vout(1)=vout(1)*tpol_proj_A1(cosm(1,1),opt_vel%pol_bend(1),opt_vel%P,opt_vel%Q)
     vout(2)=vout(2)*tpol_proj_B1(cosm(1,1),opt_vel%pol_as,opt_vel%P,opt_vel%Q)
     vout(3)=vout(3)*tpol_proj_A1(cosm(1,1),opt_vel%pol_ss(1),opt_vel%P,opt_vel%Q)
  endif

end subroutine velq2dippol
