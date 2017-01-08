!This is a file with small subroutines,
!most of time associated with math

!=============================================================
!Function which checks if an atom is inside the desired volume
!=============================================================
logical function f_in_vol(natoms,i0,a,b,c,x,y,z,tab_x,tab_y,tab_z,opt_vel)
  use mymod, only : velocity
  implicit none
  
  !From subroutine vvaf
  integer natoms,i0
  double precision a,b,c
  double precision x,y,z!Position of the studied atom
  double precision,dimension(natoms):: tab_x,tab_y,tab_z !The value x,y and z are already included into tab_x,y,z but I think it is better for the reading 
  type (velocity) :: opt_vel

  !Local
  integer i
  double precision d,dx,dy,dz
  
  !Initialization
  i=0
  d=0.d0 ; dx=0.d0 ; dy=0.d0 ; dz=0.d0
  f_in_vol=.false.

  if(opt_vel%vol_meth=='n')then!If there is no methodology to determine if an atom has to be selected, we select all of them
     f_in_vol=.true.

  elseif(opt_vel%vol_meth=='r')then!Selection according to a distance
     if(opt_vel%v0_at)then!To select the atom at t, the center of the annulus is the atom studied at t0 
        !Calculation of the ditance with periodic boundary conditions
        dx=x-tab_x(i0)
        dy=y-tab_y(i0)
        dz=z-tab_z(i0)
        dx=dx-anint(dx/a)*a
        dy=dy-anint(dy/b)*b
        dz=dz-anint(dz/c)*c
        d=sqrt(dx**2+dy**2+dz**2)
        if((d>=opt_vel%min1).and.(d<=opt_vel%max1))then!The atom in (x,y,z) will be selected only if it is in the annulus of the atom i0
           f_in_vol=.true.
        endif
        
     else!The center of the annulus is a specific atom
        do i=1,natoms
           !We check if the atom in (x,y,z) is at the good distance
           !of an atom which is selected as a center 
           if(opt_vel%tab_center(i))then
              !Calculation of the ditance with periodic boundary conditions
              dx=x-tab_x(i)
              dy=y-tab_y(i)
              dz=z-tab_z(i)
              dx=dx-anint(dx/a)*a
              dy=dy-anint(dy/b)*b
              dz=dz-anint(dz/c)*c
              d=sqrt(dx**2+dy**2+dz**2)
              if((d>=opt_vel%min1).and.(d<=opt_vel%max1))then!The atom in (X,Y,Z) will be selected only if it is in tthe annulus of the atom i
                 f_in_vol=.true.
                 exit
              endif
           endif
        enddo
     endif

  elseif(opt_vel%vol_meth=='X')then!Selection according to the the X-axis
     if(  ((opt_vel%max1>opt_vel%min1).and.(x>=opt_vel%min1) .and.(x<=opt_vel%max1)) .or.&
          ((opt_vel%max2>opt_vel%min2).and.(x>=opt_vel%min2) .and.(x<=opt_vel%max2)) .or.&
          ((opt_vel%max1<opt_vel%min1).and.((x>=opt_vel%min1).or. (x<=opt_vel%max1))).or.&
          ((opt_vel%max2<opt_vel%min2).and.((x>=opt_vel%min2).or. (x<=opt_vel%max2)))&
          )then
        f_in_vol=.true.
     endif
  elseif(opt_vel%vol_meth=='Y')then!Selection according to the the Y-axis
     if(  ((opt_vel%max1>opt_vel%min1).and.(y>=opt_vel%min1) .and.(y<=opt_vel%max1)) .or.&
          ((opt_vel%max2>opt_vel%min2).and.(y>=opt_vel%min2) .and.(y<=opt_vel%max2)) .or.&
          ((opt_vel%max1<opt_vel%min1).and.((y>=opt_vel%min1).or. (y<=opt_vel%max1))).or.&
          ((opt_vel%max2<opt_vel%min2).and.((y>=opt_vel%min2).or. (y<=opt_vel%max2)))&
          )then
        f_in_vol=.true.
     endif
  elseif(opt_vel%vol_meth=='Z')then!Selection according to the the Z-axis
     if(  ((opt_vel%max1>opt_vel%min1).and.(z>=opt_vel%min1) .and.(z<=opt_vel%max1)) .or.&
          ((opt_vel%max2>opt_vel%min2).and.(z>=opt_vel%min2) .and.(z<=opt_vel%max2)) .or.&
          ((opt_vel%max1<opt_vel%min1).and.((z>=opt_vel%min1).or. (z<=opt_vel%max1))).or.&
          ((opt_vel%max2<opt_vel%min2).and.((z>=opt_vel%min2).or. (z<=opt_vel%max2)))&
          )then
        f_in_vol=.true.
     endif

  endif

  return
end function f_in_vol



!==================================
!Projection of the norm on the bond
!==================================
subroutine vect2proj(v,bond)
  use mymod, only : info_bond
  implicit none  

  !From upper subroutine
  double precision, dimension(3):: v
  type (info_bond) :: bond

  v(1)=v(1)*bond%cos(1)+v(2)*bond%cos(2)+v(3)*bond%cos(3)
  v(2)=0.d0
  v(3)=0.d0

end subroutine vect2proj



!=============================================
!Projection of the velocity on a specific axis
!=============================================
subroutine vel_proj(axis,vout)
  implicit none  
  !From upper subroutine
  integer axis
  double precision, dimension(3):: vout
  
  !The velocity is projected along the axis "axis"
  !But the data is stored in vout(1) (for the x axis)
  vout(1)=vout(axis)
  vout(2)=0.d0
  vout(3)=0.d0

end subroutine vel_proj





!==============================================================
!Calculation of the prefactor associated with the dipole moment
!along R (if PQR poalrization)
!Depends on the molecule orientation
!==============================================================
double precision function f_dipole(c,center,at_center,ligand,third,opt_vel,z)
  use mymod, only : velocity,bond_frame
  implicit none  

  !From the main
  integer at_center,ligand,third
  double precision c,center,z
  type (velocity):: opt_vel
  
  !Local
  double precision, dimension(3) :: cosx,cosy,cosz

  !Initialization
  f_dipole=0.d0
  cosx(:)=0.d0 ; cosy(:)=0.d0 ; cosz(:)=0.d0

  !=======================================
  !Calculation of the molecule framework
  !=======================================
  call bond_frame(cosx,cosy,cosz,opt_vel,at_center,ligand,third)

  !==================
  !Dipole calculation
  !==================
  f_dipole= &
       cosx(opt_vel%R)*opt_vel%dip_stretch(1)+&
       cosy(opt_vel%R)*opt_vel%dip_stretch(2)+&
       cosz(opt_vel%R)*opt_vel%dip_stretch(3)


  !===================================================================
  !If the molecule is on a certain side, we suppose that there is
  !an inversion center 
  !Remark: an inversion center has no impact on the the polarizability
  !So the equivalent does not exist for f_polariz
  !===================================================================
  if(&
       ((center<=c/2).and.((z>=center).and.(z<center+c/2))).or.&
       ((center>c/2) .and.((z>=center).or. (z<center-c/2)))&
       )then
     f_dipole=-f_dipole
  endif

  
end function f_dipole



!==========================================================================
!Calculation of the prefactor associated with bond projection on a lab axis
!If there is a bond A-B (where A is the central atom)
!Therefore, the molecule plan will depend on the orientation of the atom C
!which make bond with A and which is the closest to B
!If there is no C atom (linear molecule) we take a random orientation
!==========================================================================
double precision function f_polariz(at_center,ligand,third,opt_vel)
  use mymod, only : velocity,bond_frame
  implicit none  
  
  !From the main
  integer at_center,ligand,third
  type (velocity):: opt_vel

  !Local
  integer i,j
  double precision tmp
  double precision, dimension(3) :: cosPi,cosQj,cosx,cosy,cosz

  !Initialization
  i=0 ; j=0
  tmp=0.d0 ; f_polariz=0.d0
  cosPi(:)=0.d0 ; cosQj(:)=0.d0 ; cosx(:)=0.d0 ; cosy(:)=0.d0 ; cosz(:)=0.d0

  !=======================================
  !Calculation of the molecule framework
  !=======================================
  call bond_frame(cosx,cosy,cosz,opt_vel,at_center,ligand,third)

  !=====================================================
  !Attribution of the the good values to cosPi and cosQj
  !=====================================================
  cosPi(1)=cosx(opt_vel%P) ; cosPi(2)=cosy(opt_vel%P) ; cosPi(3)=cosz(opt_vel%P)
  cosQj(1)=cosx(opt_vel%Q) ; cosQj(2)=cosy(opt_vel%Q) ; cosQj(3)=cosz(opt_vel%Q)

  !========================
  !Polarization calculation
  !========================
  do i=1,3
     tmp=0.d0
     do j=1,3
        tmp=tmp+cosQj(j)*opt_vel%pol_stretch(j,i)
     enddo
     f_polariz=f_polariz+cosPi(i)*tmp
  enddo


end function f_polariz



!======================================================
!Subroutine which dertermines the bond framework
!which corresponds with the direction cosine matrix
!The inversion center of the cell is taken into account
!xyz is the bond framework projected on XYZ
!the lab framework
!======================================================
subroutine bond_frame(cosx,cosy,cosz,opt_vel,at_center,ligand,third)
  use mymod, only : velocity,cross_prod
  implicit none

  !From upper subroutine
  integer at_center,ligand,third
  double precision, dimension(3) :: cosx,cosy,cosz
  type (velocity):: opt_vel

  !====
  !cosz
  !====
  cosz(:)=opt_vel%tab_mol(ligand,at_center)%cos

  !===============================================
  !cosy
  !A guess is used depending on the number of atom
  !===============================================
  if(third/=0)then!If the molecule contains 3 atoms at least
     cosx(:)=opt_vel%tab_mol(third,at_center)%cos
  else!The molecule has only 2 atoms
     !The vector (1,0,0) is used as a guess to calculate cosy()
     cosx(1)=1.d0
  endif
  call cross_prod(cosz,cosx,cosy)

  !====
  !cosx
  !====
  call cross_prod(cosy,cosz,cosx)
end subroutine bond_frame



!======================================================
!Subroutine which dertermines the molecular framework
!which corresponds with the direction cosine matrix
!The inversion center of the cell is taken into account
!xyz is the molecular framework projected on XYZ
!the lab framework
!======================================================
subroutine mol_frame(cosb1,cosb2,n1,n2,cosm)
  use mymod, only : error, cross_prod
  implicit none

  !Declaration from upper subroutine
  double precision, dimension(3):: cosb1, cosb2, n1, n2
  double precision, dimension(3,3) :: cosm

  !Local variables
  double precision :: norm

  !Initialization
  norm=0.d0
  
  !============================================
  !Determination of the direction cosine matrix
  !============================================
  !Case of a linear molecule
  if( abs( &
       cosb1(1)*cosb2(1) + &
       cosb1(2)*cosb2(2) + &
       cosb1(3)*cosb2(3) &
       ) > 1-error )then

     if( cosb1(3)**2 > 1-error) then !The molecule is fully along Z (lab framework)
        !Basic direction cosine matrix (with z along X, x along Y and y along z)
        cosm(1,3)=1.d0 !Xz
        cosm(2,3)=0.d0 !Yz
        cosm(3,3)=0.d0 !Zz

        !To be coherent with the general case, x is along the first bond
        if(cosb1(3)>0) then
           cosm(1,1)=0.d0 ; cosm(1,2)= 0.d0 !Xx Xy
           cosm(2,1)=0.d0 ; cosm(2,2)=-1.d0 !Yx Yy
           cosm(3,1)=1.d0 ; cosm(3,2)= 0.d0 !Zx Zy
        else
           cosm(1,1)= 0.d0 ; cosm(1,2)=0.d0 !Xx Xy
           cosm(2,1)= 0.d0 ; cosm(2,2)=1.d0 !Yx Yy
           cosm(3,1)=-1.d0 ; cosm(3,2)=0.d0 !Zx Zy
        end if

        !n1 = n2 = -z = -X
        n1(1)=-1.d0 ; n2(1)=-1.d0
        n1(2)= 0.d0 ; n2(2)= 0.d0
        n1(3)= 0.d0 ; n2(3)= 0.d0

     else !If not
        !1) x is along the first bond
        cosm(:,1)=cosb1(:)

        !2) z=x^Z
        norm=sqrt(cosm(1,1)**2+cosm(2,1)**2)
        cosm(1,3)= cosm(2,1)/norm!Xz
        cosm(2,3)=-cosm(1,1)/norm!Yz
        cosm(3,3)=0.d0           !Zz

        !3) y=z^x
        call cross_prod(cosm(1,3),cosm(1,1),cosm(1,2))

        !n1 = n2 = -z
        n1(1)=-cosm(1,3) ; n2(1)=-cosm(1,3)
        n1(2)=-cosm(2,3) ; n2(2)=-cosm(2,3)
        n1(3)= 0.d0      ; n2(3)= 0.d0

     end if

  else !Non linear water molecule
     !z (molecular framework) is along the first bissector
     cosm(:,3)=cosb1(:) + cosb2(:)
     norm=sqrt( cosm(1,3)**2 + cosm(2,3)**2 + cosm(3,3)**2 )
     cosm(:,3)=cosm(:,3)/norm

     !x is in plane
     !IMPORTANT REMARK: x.OH2 > 0
     !One may have some sign problems if the orientation of x is the opposite
     cosm(:,1)= -cosb1(:) + cosb2(:)
     norm=sqrt( cosm(1,1)**2 + cosm(2,1)**2 + cosm(3,1)**2 )
     cosm(:,1)=cosm(:,1)/norm

     !y is out of the plane
     call cross_prod(cosm(1,3),cosm(1,1),cosm(1,2))

     !n1 =-y_mol ^ z_bond1
     call cross_prod(cosm(1,2),cosb1(1),n1(1))
     n1(:)=-n1(:)
     !n2 =y_mol ^ z_bond2
     call cross_prod(cosm(1,2),cosb2(1),n2(1))

  end if
  
end subroutine mol_frame







!======================================================================================================
!Subroutine which determines which atoms are inside a molecule
!The atoms have to respect the rules defined in the block BOND_WITH
!When the atom i and j (i<j and nb_bond=2) makes bonds with the atom k: opt_vel%bond(k)%index=j
!Atomic diatances will be calculated for all the central atoms (opt_vel%tab_select(i))
!If all the criteria are respected:
!*good symbol: tab_at_bond(j)==.true.
!*dmin <= d <= dmax
!*i_ligand==nb_bond (at the end)
!therefore tab_mol(i,:) is filled
!======================================================================================================
subroutine make_molecule(natoms,x,y,z,a,b,c,opt_vel,tab_new_mol,record)
  use mymod, only : velocity,third_at,bond_update
  implicit none

  !From vvaf subroutine
  integer natoms
  double precision a,b,c
  double precision, dimension(natoms) :: x,y,z
  logical :: record !true if t0, false if t
  logical, dimension(natoms) :: tab_new_mol
  type (velocity) :: opt_vel

  !Local
  integer i,j,i_ligand
  double precision dx,dy,dz,d2

  !Initialization
  i=0 ; j=0 ; i_ligand=0
  dx=0.d0 ; dy=0.d0 ; dz=0.d0 ; d2=0.d0

  !Loop of the central atoms
  do i=1,natoms
     if(opt_vel%tab_select(i))then
        do j=1,natoms
           if(opt_vel%tab_at_bond(j).and.(j/=i))then!j has to be allowed to make bonds with i (and be different from i)
              dx=x(j)-x(i)
              dy=y(j)-y(i)
              dz=z(j)-z(i)
              dx=dx-anint(dx/a)*a
              dy=dy-anint(dy/b)*b
              dz=dz-anint(dz/c)*c
              d2=dx**2+dy**2+dz**2

              !Is d2 in the bonds?
              if((opt_vel%dmin2<=d2).and.(opt_vel%dmax2>=d2))then
                 if(i_ligand<opt_vel%nb_bond)then!The upper limit of ligand is not yet reached
                    i_ligand=i_ligand+1
                    !The ligands are different at t0-shift and t0, this is a new molecule
                    if(opt_vel%tab_mol(i_ligand,i)%index/=j)then
                       opt_vel%tab_mol(i_ligand,i)%index=j
                       tab_new_mol(i)=.true.
                    endif
                 else!If there are too much ligands, the atom i does not do any molecule
                    i_ligand=0
                    exit
                 endif!Number of ligand
              endif!Bond with good length

           endif!Atom j can make bonds
        enddo!Ligand loop (j)

        
        if(i_ligand==opt_vel%nb_bond)then
           !The molecule is good, therefore we have to:
           !1)define the 3rd atom which will allow to build the molecule framework
           !2)calculate the bond orientation (for v0 only)
           do j=1,opt_vel%nb_bond
              call third_at(natoms,opt_vel,i,j,x,y,z,a,b,c)
              if(record) call bond_update(&
                   x(opt_vel%tab_mol(j,i)%index)-x(i),&
                   y(opt_vel%tab_mol(j,i)%index)-y(i),&
                   z(opt_vel%tab_mol(j,i)%index)-z(i),&
                   opt_vel%tab_mol(j,i),a,b,c)
           enddo
           i_ligand=0
        elseif(i_ligand<opt_vel%nb_bond)then!There are not enougth bonds (or too much because i_ligand has been set to 0)
           opt_vel%tab_mol(:,i)%index=0
           i_ligand=0
        endif

     endif
  enddo!Central atom loop (i)

end subroutine make_molecule



!=======================================================================================
!Subroutine which determines which atom k will be used to define
!xyz (the molecule framework)
!There are already 2 atoms (the central atom i and the bounded atom j).
!The 3rd atom is those which is the closest from the atom j (and different from i and j)
!The bond will be along z
!The atoms i j and k will be in Oxz
!=======================================================================================
subroutine third_at(natoms,opt_vel,i,j,x,y,z,a,b,c)
  use mymod, only : velocity
  implicit none

  !From the vvaf subroutine
  integer :: natoms,i,j
  double precision:: a,b,c
  double precision, dimension(natoms) :: x,y,z
  type (velocity) :: opt_vel

  !Local
  integer k
  double precision dx,dy,dz,d2,d2min

  !Initialization
  opt_vel%tab_mol(j,i)%third=0
  k=0
  dx=0.d0 ; dy=0.d0 ; dz=0.d0 ; d2=0.d0 ; d2min=1000000000000.d0

  do k=1,opt_vel%nb_bond
     if(k/=j)then
        !Distance calculation
        dx=x(opt_vel%tab_mol(k,i)%index)-x(opt_vel%tab_mol(j,i)%index)
        dy=y(opt_vel%tab_mol(k,i)%index)-y(opt_vel%tab_mol(j,i)%index)
        dz=z(opt_vel%tab_mol(k,i)%index)-z(opt_vel%tab_mol(j,i)%index)
        
        !pbc
        dx=dx-anint(dx/a)*a
        dy=dy-anint(dy/b)*b
        dz=dz-anint(dz/c)*c

        !Minimum
        if(d2min>d2)then
           opt_vel%tab_mol(j,i)%third=k
           d2min=d2
        endif
     endif
  enddo
  
end subroutine third_at





!======================================================================
!Subroutine which calculates the projection of the velocity on the bond
!======================================================================
subroutine stretch(vout,bond)
  use mymod, only : info_bond
  implicit none

  !From the vvaf subroutine
  double precision, dimension(3) :: vout
  type (info_bond) :: bond

  !Local
  double precision dotp

  !Initialization
  dotp=0.d0

  !Dot product normalized according to the bond length
  dotp=vout(1)*bond%cos(1)+vout(2)*bond%cos(2)+vout(3)*bond%cos(3)

  !Projection on X,Y,Z (lab framework)
  vout(1)=dotp*bond%cos(1)
  vout(2)=dotp*bond%cos(2)
  vout(3)=dotp*bond%cos(3)
  
end subroutine stretch




!================================================================
!Subroutine which updates options%v?%bond
!The atoms involved are the same but
!dist, dx, dy and dz are updated 
!IMPORTANT REMARK
!The projected cosine are not calculated now,
!only when a polarizability calculation is requiered
!However, to calculate such matrix, it is important to know first
!the direction of all the molecule bonds
!================================================================
subroutine bond_update(dx,dy,dz,bond,a,b,c)
  use mymod, only : info_bond
  implicit none

  !From vvaf subroutine
  double precision a,b,c,dx,dy,dz
  type (info_bond) :: bond

  !Projection of the bond on the normal axis
  bond%cos(1)=dx-anint(dx/a)*a
  bond%cos(2)=dy-anint(dy/b)*b
  bond%cos(3)=dz-anint(dz/c)*c

  !Distance
  bond%dist=sqrt(bond%cos(1)**2+bond%cos(2)**2+bond%cos(3)**2)

  !Normalization (cosine)
  bond%cos=bond%cos/bond%dist

end subroutine bond_update



!==========================================================
!Subroutine which calculates a difference of velocities
!if it is requiered by options%velocity
!The velocities from the velocity file are not modified
!No velocity difference is done if options%velocity=.false.
!==========================================================
subroutine vel_diff(opt_vel,vatx,vaty,vatz,vrefx,vrefy,vrefz,vout)
  use mymod, only : velocity
  implicit none

  !From vvaf suboutine
  double precision vatx,vaty,vatz,vrefx,vrefy,vrefz
  double precision,dimension(3):: vout
  type (velocity) :: opt_vel

  if(opt_vel%bond_diff)then
     vout(1)=vatx-vrefx
     vout(2)=vaty-vrefy
     vout(3)=vatz-vrefz
  else
     vout(1)=vatx
     vout(2)=vaty
     vout(3)=vatz
  endif

end subroutine vel_diff






!=================================================================================
!Subroutine  which does some velocity treatments associated with bonds
!-Update of the bonds at t (the molecules are defined at each t0)
!-Vecocity difference between two atoms
!-Calculation of the dipole moment prefactor (projectiopn of the bond on the axis)
!-Calculation of polarizability (with our level of theory: only 0 or 1)
!=================================================================================
!i is the central atom (those which is studied by the f_in_vol)
!j is the number of the bond
subroutine bond_only(natoms,opt_vel,i,ligand,vout,vx,vy,vz,x,y,z,box_size,center,slab_orient)
  use mymod, only : velocity,vel_diff,stretch,f_dipole,f_polariz
  implicit none
  
  !From upper subroutine
  integer natoms,i,ligand
  double precision box_size,center
  double precision, dimension(3):: vout
  double precision, dimension(natoms)::vx,vy,vz,x,y,z
  character slab_orient
  type (velocity) :: opt_vel

  !Local
  integer j
  double precision ptr_z

  !Initialization
  ptr_z=0.d0
  j=opt_vel%tab_mol(ligand,i)%index!To simplify the reading


  !Calculation of the velocity difference (bond velocity instead of atomic velocity)
  call vel_diff(opt_vel,vx(j),vy(j),vz(j),&
       vx(i),vy(i),vz(i),&
       vout(1))
  
  !Projection of the velocity on the bond (stretching only)
  !If only the stretching mode is studied, then we do a projection of the velocity over the bonds
  !Mathematically: projection of vout on bond%cos
  if(opt_vel%stretch_mode)then
     call stretch(vout,opt_vel%tab_mol(ligand,i))
  endif

  !Calculation of the projection of the molecule on the normal axis
  !Necessary for the sign (due to the cell symmetry) of the prefactor
  if(slab_orient=='X')then
     ptr_z = x(i)
  elseif(slab_orient=='Y')then
     ptr_z = y(i)
  else
     ptr_z = z(i)
  end if

  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !Calculation of the prefactors
  !Remark: the two prefactors (dipole moment and polarizability)
  !can be used at same time
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !Calculation of the prefactor for dipole moment
  !Projection of the bond on the R axis (if PQR polarization)
  if(opt_vel%R/=0)then
     vout=vout*f_dipole(box_size,center,i,ligand,opt_vel%tab_mol(ligand,i)%third,opt_vel,ptr_z)
  endif
  !Calculation of the prefactor for the polarizability
  !Projections of xyz (molecule framework) in PQ (axis of XYZ -- lab framework)
  if(opt_vel%P/=0)then
     vout=vout*f_polariz(i,ligand,opt_vel%tab_mol(ligand,i)%third,opt_vel)
  endif
  
end subroutine bond_only


!=============
!Cross product
!=============
subroutine cross_prod(v1,v2,v3)
  implicit none

  !From upper subroutine
  double precision, dimension(3)::v1,v2,v3

  !local
  double precision norm
  
  !Initialization
  norm=0.d0

  !Vectorial product
  v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
  v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
  v3(3)=v1(1)*v2(2)-v1(2)*v2(1)

  !Normalization
  norm=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
  if(norm>epsilon(0d0))then
     v3=v3/norm
  else
     v3=0.d0
  endif

end subroutine cross_prod



!==============================================================
!This subroutine is especially done in order to multiply
!a rotation matrix by a vector having only one non-nul element
!and extracting a single element of the resulting vector
!
!E.g. (Rotation matrix) x (transition dipole moment of bending along z)
!projected on R (from PQR polarization)
!==============================================================
double precision function tdip_proj(cosm,dip,R)
  implicit none  

  !From upper subroutine
  integer R !Polarization of the IR beam (X=1, Y=2, Z=3, nothing=0)
  double precision dip !transition dipole moment
  double precision, dimension(3) :: cosm !Column of the direction cosine matrix

  tdip_proj=dip*cosm(R)

  
end function


!==============================================================
!This subroutine is especially done in order to rotate
!a diagonal tensor (from a normal mode)
!and extracting a single element of the resulting vector
!alpha_PQ = Px alpha_xx xQ + Py alpha_yy yQ + Pz alpha_zz zQ 
!
!E.g. Element PQ (from PQR polarization)
!(Rotation matrix) x (transition polarizability of bending along z) x (Rotation matrix)t
!==============================================================
double precision function tpol_proj_A1(cosm,pol,P,Q)
  implicit none  

  !From upper subroutine
  integer ::P !First element of the polarizability tensor (X=1, Y=2, Z=3, nothing=0)
  integer ::Q !Second element of the polarizability tensor (X=1, Y=2, Z=3, nothing=0)
  double precision, dimension(3) :: pol ! transition polarizability (only three values are requiered)
  double precision, dimension(3,3) :: cosm ! direction cosine matrix

  tpol_proj_A1= &
       pol(1)*cosm(P,1)*cosm(Q,1)+ &
       pol(2)*cosm(P,2)*cosm(Q,2)+ &
       pol(3)*cosm(P,3)*cosm(Q,3)
  
end function tpol_proj_A1
  
!==============================================================
!The same for B1 symmetry (alpha_XZ=alpha_ZX, the rest is 0)
!alpha_PQ = alpha_xz (Pz xQ + Px zQ)
!==============================================================
double precision function tpol_proj_B1(cosm,pol,P,Q)
  implicit none  

  !From upper subroutine
  integer ::P !First element of the polarizability tensor (X=1, Y=2, Z=3, nothing=0)
  integer ::Q !Second element of the polarizability tensor (X=1, Y=2, Z=3, nothing=0)
  double precision :: pol ! transition polarizability (only three values are requiered)
  double precision, dimension(3,3) :: cosm ! direction cosine matrix

  tpol_proj_B1= pol *( &
       cosm(P,1)*cosm(Q,3) + &
       cosm(P,3)*cosm(Q,1)   &
       )
  
end function tpol_proj_B1



