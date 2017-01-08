!This file contaiin the subroutine doing the loops for
!in order to do the average <X(t).Y(0)>

!================
!VVAF calculation
!================
subroutine vvaf(nmovie,natoms,Di,shift,options,vx,vy,vz,x,y,z,slab_orient,center,a,b,c,&
     box_size,surface,corr)
  use mymod, only : vvaf_options,bond_calc,make_molecule,f_in_vol,bond_update,simple_c2v_calc,&
       hessian_c2v_calc
  implicit none

  !From the main
  integer nmovie,natoms,Di,shift
  double precision,dimension(natoms,nmovie):: vx,vy,vz,x,y,z
  double precision a,b,c,box_size,center,surface
  double precision,dimension(-Di:Di,9)::corr
  character slab_orient
  type (vvaf_options) :: options

  !Local
  integer :: i,j,k,l,m,t,t_vt,t0,nb_shift,nb_v0_mol,nb_vt_mol
  integer, dimension(natoms):: t_last
  double precision, dimension(3,3) :: v0
  double precision, dimension(:,:), allocatable :: vt_all,vt_mol
  double precision, dimension(:,:,:,:), allocatable :: vt_bond
  logical, dimension(natoms) :: tab_new_mol!True if the molecule at t0-shift and t0 is different
  
  !Arrays
  !corr(:,1) is the autocorrelation (the first index is for t)
  !corr(:,2) is the intramolecular correlation 
  !corr(:,3) is the intermolecular correlation 
  !The same for v0(:,:) (the first index is for x/y/z)
  !For vt it is more complicated. There are 3 vt array for the:
  !  1)bond auto-correlation
  !  2)intramolecular cross corelation
  !  3)inter-molecular cross correlation
  !The order of the index are almost the same:
  !  1)x/y/z
  !  2)t
  !  3)the bonds (disapears for intra/inter-molecular correlation)
  !  4)the central atom (disapears for inter-molecular correlation)
  if(options%vt%nb_bond==0)then
     allocate(&
          vt_bond(3,-Di:Di,1,natoms),&
          vt_mol(3,-Di:Di),& !Useless if no bond but requiered to not add one more "if blocks" inside the code
          vt_all(3,-Di:Di)& !idem
          )
  elseif(options%v0%c2v_mode>0) then !Special case for c2v molecules
     !Only 3 time derivative are requiered: bending (along z), AS (along x), SS (along z)
     allocate(&
          vt_bond(3,-Di:Di,options%v0%nb_bond,natoms),&!
          vt_mol(0,0),&
          vt_all(0,0)&
          )
  else
     allocate(&
          vt_bond(3,-Di:Di,options%v0%nb_bond,natoms),&
          vt_mol(3,-Di:Di),&
          vt_all(3,-Di:Di)&
          )
  end if

  !Initialization
  i=0 ; j=0 ; k=0 ; l=0 ; m=0 ; t=0; t_vt=0 ; t0=0 ; nb_shift=0 ; nb_v0_mol=0 ; nb_vt_mol=0
  v0=0.d0  ; vt_bond(:,:,:,:)=0.d0
  if(options%vt%nb_bond/=0)then
     vt_mol(:,:)=0.d0 ; vt_all(:,:)=0.d0
  endif
  tab_new_mol(:)=.false.
  open(11,file=options%filename_tmp)

  !============
  !Optimization
  !============
  !Final step recorded for vt with the previous t0
  !So if the molecule is the same between [t0-shift:t0], it is useless to recalculate vt
  !between [t0-Di:t_last]
  !BE CAREFULL, "t_last/=t0-shift+Di" because of the very first t0 step
  t_last(:)=0

  !==========================
  !Loop for t0 and the atom i
  !==========================
  do t0=Di+1,nmovie-Di,shift
     nb_shift=nb_shift+1!For VVAF normalization

     !Molecule definition
     !It is done once for bothe velocities at each t0 even if we studying the stretching at t...
     !... because  it is meaningless to study the study the stretching of 2 dfferent bonds
     !e.g. O1-H + O2 (at t0) -> O1 + H-O2 (at t)
     !The stetching of O1-H compared with O2-H is meaningless, so we keep the same molecule between t0-Dt and t0+Dt
     !Molecules for v0
     if(options%v0%nb_bond>0)then
        !If we need to study dipole moment or stretching or velocities difference or the polarizability, we need to calculate bonds
        call make_molecule(natoms,x(1,t0),y(1,t0),z(1,t0),a,b,c,options%v0,tab_new_mol(1),.true.)!What is the closest atom among options%v0%tab_at_bond (except the i itself)?
        nb_v0_mol=nb_v0_mol+count(options%v0%tab_mol(1,:)%index>0)!Number of molecule taken into account (important to check results)
     endif
     !Molecules for vt
     if(options%vt%nb_bond>0)then
        if(.not.options%cross)then!Auto-correlation
           options%vt%tab_mol=options%v0%tab_mol
        else!Cross correlation (molecules for v(0) and v(t) can be different)
           call make_molecule(natoms,x(1,t0),y(1,t0),z(1,t0),a,b,c,options%vt,tab_new_mol(1),.false.)
        endif
        nb_vt_mol=nb_vt_mol+count(options%vt%tab_mol(1,:)%index>0)!Number of molecule taken into account (important to check results)
     endif

     i=1!Associated with the central atoms
     j=1!Associated with the bonds
     do while (i<=natoms)
        !We select only the atoms requiered in options%v0%tab_select
        if(options%v0%tab_select(i))then
           !Selection of the atoms in the good volume
           if(f_in_vol(natoms,i,a,b,c,x(i,t0),y(i,t0),z(i,t0),x(1,t0),y(1,t0),z(1,t0),options%v0))then

              !Calculation of the velocity at t0 for the i atom with all the prefactors
              v0(:,1)=0.d0
              if(options%v0%c2v_mode==1)then
                 call simple_c2v_calc(natoms,i,box_size,center,slab_orient,options%v0,v0(1,1),vx(1,t0),vy(1,t0),vz(1,t0),&
                      x(1,t0),y(1,t0),z(1,t0))
              elseif(options%v0%c2v_mode==2) then
                 call hessian_c2v_calc(natoms,i,box_size,center,slab_orient,options%v0,v0(1,1),vx(1,t0),vy(1,t0),vz(1,t0),&
                      x(1,t0),y(1,t0),z(1,t0))
              else
                 call bond_calc(natoms,i,j,box_size,center,slab_orient,options%v0,v0(1,1),vx(1,t0),vy(1,t0),vz(1,t0),&
                   x(1,t0),y(1,t0),z(1,t0))
                 v0(:,2)=v0(:,2)+v0(:,1)!Calculation of the molecular contribution (sum of the bonds contribution)
              endif

              !=========================
              !Loop for t and the atom k
              !=========================
              do t = t0-Di,t0+Di
                 t_vt=modulo(t,2*Di+1)-Di!increment for the array vt_bond(:,t_vt,:,:). Increase the speed if "shift<2*Di+1" 

                 !WORK ARBEIT TRAVAIL
                 !THINK ABOUT A 2D TAB LIKE PREVIOUSLY
                 !BUT IT HAS TO BE ADJSUTED IF THE CROSS PARAMETER IS NOT SELECTED

                 !+++++++++++++++++++++++++++++++++++++++++++++++++
                 !General case: Cross correlation (time consumming)
                 !+++++++++++++++++++++++++++++++++++++++++++++++++
                 if(options%cross)then
                    k=1!Associated with the central atom
                    l=1!Associated with the bonds
                    do while (k<=natoms)

                       !We select only the atoms requiered in options%vt%tab_select
                       if(options%vt%tab_select(k))then
                          !Selection of the atoms in the good volume
                          if(f_in_vol(natoms,i,a,b,c,x(k,t),y(k,t),z(k,t),x(1,t),y(1,t),z(1,t),options%vt))then
                             if((t>t_last(k)).or.tab_new_mol(k))then!It is useless to recalculate vt if it is already calculated for the previous t0
                                if(options%vt%tab_mol(1,k)%index>0)then
                                   !Update of the bonds
                                   !The bonds are recorded in "make_molecule" subroutine only for v0
                                   do m=1,options%vt%nb_bond
                                      call bond_update(&
                                           x(options%vt%tab_mol(m,k)%index,t)-x(k,t),&
                                           y(options%vt%tab_mol(m,k)%index,t)-y(k,t),&
                                           z(options%vt%tab_mol(m,k)%index,t)-z(k,t),&
                                           options%vt%tab_mol(m,k),a,b,c)
                                   enddo
                                   !t_last(k)=t!As soon as the bond-update is done, we update the time of the last update 
                                endif
                                !Calculation of the velocity at t for the j atom with all the prefactors
                                if(options%v0%c2v_mode==1) then
                                   call simple_c2v_calc(natoms,k,box_size,center,slab_orient,options%vt,vt_bond(1,t_vt,l,k),&
                                        vx(1,t),vy(1,t),vz(1,t),x(1,t),y(1,t),z(1,t))
                                elseif(options%v0%c2v_mode==2) then
                                   call hessian_c2v_calc(natoms,i,box_size,center,slab_orient,options%vt,vt_bond(1,t_vt,l,k),&
                                        vx(1,t),vy(1,t),vz(1,t),x(1,t),y(1,t),z(1,t))
                                else
                                   call bond_calc(natoms,k,l,box_size,center,slab_orient,options%vt,&
                                        vt_bond(1,t_vt,l,k),vx(1,t),vy(1,t),vz(1,t),x(1,t),y(1,t),z(1,t))
                                end if

                             endif
                             
                             if(options%v0%c2v_mode>0) then
                                do m =1,3 !1=X_B, 2=X_AS, 3=X_SS (X may be mu or alpha)
                                   corr(t-t0,3*m-2) = corr(t-t0,3*m-2) + v0(m,1)*vt_bond(1,t_vt,l,k) !Y_B
                                   corr(t-t0,3*m-1) = corr(t-t0,3*m-1) + v0(m,1)*vt_bond(2,t_vt,l,k) !Y_AS
                                   corr(t-t0,3*m)   = corr(t-t0,3*m)   + v0(m,1)*vt_bond(3,t_vt,l,k) !Y_SS
                                end do
                                
                                k=k+1!Next central atom
                                l=1!Restart to read the bonds

                             else 
                                corr(t-t0,1) = corr(t-t0,1) &
                                  + v0(1,1)*vt_bond(1,t_vt,l,k) &
                                  + v0(2,1)*vt_bond(2,t_vt,l,k) &
                                  + v0(3,1)*vt_bond(3,t_vt,l,k)!auto-correlation

                                !Once all the bonds are read for the atom k
                                if(l>=options%vt%nb_bond)then
                                   k=k+1!Next central atom
                                   l=1!Restart to read the bonds
                                else
                                   l=l+1
                                endif
                             endif

                          else!k is not in the volume
                             k=k+1
                             l=1
                          endif
                       else!k is not a central atom
                          k=k+1
                          l=1
                       endif
                    enddo!do while k

                    !+++++++++++++++++++++++++++++++++++++
                    !C2V symmetry (bending, AS, SS modes)
                    !The modes are simply defined on theta
                    !r1+r2-2rO and r1-r2
                    !+++++++++++++++++++++++++++++++++++++
                 elseif(options%v0%c2v_mode==1) then
                    if((t>t_last(i)).or.tab_new_mol(i))then!It is useless to recalculate vt if it is already calculated for the previous t0
                       !Update of the bonds if the atom i is a central atom
                       !The bonds are recorded in "make_molecule" subroutine only for v0
                       if(options%vt%tab_mol(1,i)%index>0)then
                          do m=1,options%vt%nb_bond
                             call bond_update(&
                                  x(options%vt%tab_mol(m,i)%index,t)-x(i,t),&
                                  y(options%vt%tab_mol(m,i)%index,t)-y(i,t),&
                                  z(options%vt%tab_mol(m,i)%index,t)-z(i,t),&
                                  options%vt%tab_mol(m,i),a,b,c)
                          enddo
                       end if
                       !Calculation of the velocity at t for the i atom with all the prefactors
                       call simple_c2v_calc(natoms,i,box_size,center,slab_orient,options%vt,vt_bond(1,t_vt,j,i),&
                         vx(1,t),vy(1,t),vz(1,t),x(1,t),y(1,t),z(1,t))
                    end if

                    do k =1,3 !1=X_B, 2=X_AS, 3=X_SS (X may be mu or alpha)
                       corr(t-t0,3*k-2) = corr(t-t0,3*k-2) + v0(k,1)*vt_bond(1,t_vt,j,i) !Y_B
                       corr(t-t0,3*k-1) = corr(t-t0,3*k-1) + v0(k,1)*vt_bond(2,t_vt,j,i) !Y_AS
                       corr(t-t0,3*k)   = corr(t-t0,3*k)   + v0(k,1)*vt_bond(3,t_vt,j,i) !Y_SS
                    end do

                    
                    !++++++++++++++++++++++++++++++++++++
                    !C2V symmetry (bending, AS, SS modes)
                    !The modes are based on the Hessian
                    !++++++++++++++++++++++++++++++++++++
                 elseif(options%v0%c2v_mode==2) then
                    if((t>t_last(i)).or.tab_new_mol(i))then!It is useless to recalculate vt if it is already calculated for the previous t0
                       !Update of the bonds if the atom i is a central atom
                       !The bonds are recorded in "make_molecule" subroutine only for v0
                       if(options%vt%tab_mol(1,i)%index>0)then
                          do m=1,options%vt%nb_bond
                             call bond_update(&
                                  x(options%vt%tab_mol(m,i)%index,t)-x(i,t),&
                                  y(options%vt%tab_mol(m,i)%index,t)-y(i,t),&
                                  z(options%vt%tab_mol(m,i)%index,t)-z(i,t),&
                                  options%vt%tab_mol(m,i),a,b,c)
                          enddo
                       end if
                       !Calculation of the velocity at t for the i atom with all the prefactors
                       call hessian_c2v_calc(natoms,i,box_size,center,slab_orient,options%vt,vt_bond(1,t_vt,j,i),&
                         vx(1,t),vy(1,t),vz(1,t),x(1,t),y(1,t),z(1,t))
                    end if

                    do k =1,3 !1=X_B, 2=X_AS, 3=X_SS (X may be mu or alpha)
                       corr(t-t0,3*k-2) = corr(t-t0,3*k-2) + v0(k,1)*vt_bond(1,t_vt,j,i) !Y_B
                       corr(t-t0,3*k-1) = corr(t-t0,3*k-1) + v0(k,1)*vt_bond(2,t_vt,j,i) !Y_AS
                       corr(t-t0,3*k)   = corr(t-t0,3*k)   + v0(k,1)*vt_bond(3,t_vt,j,i) !Y_SS
                    end do

                    
                    !++++++++++++++++++++++++++++++++++++++
                    !Autocorrelation (less time consumming)
                    !++++++++++++++++++++++++++++++++++++++
                 else
                    if((t>t_last(i)).or.tab_new_mol(i))then!It is useless to recalculate vt if it is already calculated for the previous t0
                       !Update of the bonds if the atom i is a central atom
                       !The bonds are recorded in "make_molecule" subroutine only for v0
                       if(options%vt%tab_mol(1,i)%index>0)then
                          do m=1,options%vt%nb_bond
                             call bond_update(&
                                  x(options%vt%tab_mol(m,i)%index,t)-x(i,t),&
                                  y(options%vt%tab_mol(m,i)%index,t)-y(i,t),&
                                  z(options%vt%tab_mol(m,i)%index,t)-z(i,t),&
                                  options%vt%tab_mol(m,i),a,b,c)
                          enddo
                       end if
                       !Calculation of the velocity at t for the i atom with all the prefactors
                       call bond_calc(natoms,i,j,box_size,center,slab_orient,options%vt,&
                            vt_bond(1,t_vt,j,i),vx(1,t),vy(1,t),vz(1,t),x(1,t),y(1,t),z(1,t))
                    end if

                    vt_mol(:,t_vt)=vt_mol(:,t_vt)+vt_bond(:,t_vt,j,i)!Calculation of the molecular contribution (sum of the bonds contribution)

                    corr(t-t0,1) = corr(t-t0,1) &
                         + v0(1,1)*vt_bond(1,t_vt,j,i) &
                         + v0(2,1)*vt_bond(2,t_vt,j,i) &
                         + v0(3,1)*vt_bond(3,t_vt,j,i)!Auto-correlation
                    
                    
                 endif!Cross-/Auto-correlation
                    
              enddo!t
              
              !The bond j of the atom i has been treated for every t
              if(options%v0%c2v_mode>0)then !Nothing has to be done if we study the normal modes of water (both bonds have been treated)
                 i=i+1
              else if(j>=options%v0%nb_bond)then !All the bonds have been treated
                 do t = -Di,Di
                    t_vt=modulo(t0+t,2*Di+1)-Di!Calculation of the increment for vt
                    corr(t,2) = corr(t,2) + v0(1,2)*vt_mol(1,t_vt) + v0(2,2)*vt_mol(2,t_vt) + v0(3,2)*vt_mol(3,t_vt)!Intra-molecular correlation

                    !Calculation of the inter-molecular contribution (sum of the molecular contributions)
                    !v0 is out of the do-loop
                    vt_all(:,t_vt)=vt_all(:,t_vt)+vt_mol(:,t_vt)

                 enddo
                 v0(:,3)=v0(:,3)+v0(:,2)        

                 !Re-initialization
                 v0(:,2)=0.d0
                 vt_mol(:,:)=0.d0

                 !Once the molecule is treated, it is not necessary anymore
                 !To consider it as a NEW molecule
                 tab_new_mol(i)=.false.

                 i=i+1!Next central atom
                 j=1!Restart to read the bonds

              else !We pass to the next bond of the atom i
                 j=j+1
              endif

           else!The central atom i is not in the volume
              !If the molecule is not in the volume
              !the molecule will be consider as new for next time
              tab_new_mol(i)=.true.
 
              i=i+1!Next central atom
              j=1!Restart to read the bonds
           endif!f_in_vol
        else
           i=i+1!Next central atom
           j=1!Restart to read the bonds
        endif!options%v0%tab_select(i)
     end do!do while i

     !Full-molecular contribution
     if((options%v0%c2v_mode==0) .and. (.not. options%cross)) then
        do t = -Di,Di
           t_vt=modulo(t0+t,2*Di+1)-Di!Calculation of the increment for vt
           corr(t,3) = corr(t,3) + v0(1,3)*vt_all(1,t_vt) + v0(2,3)*vt_all(2,t_vt) + v0(3,3)*vt_all(3,t_vt)
        end do
     end if

     !Re-initialization
     v0(:,3)=0.d0
     vt_all(:,:)=0.d0
     !Update of t_last
     t_last(:)=t0+Di

  enddo!t0

  !============
  !VVAF writing
  !============
  !BE CAREFULL, DO NOT store the normalized data (if there is a layering, the normalization will be done several times and the result will be wrong)
  if(options%v0%c2v_mode>0) then !Bending / Antisymmetric / Symmetric modes
     write(11,'(a)',advance='no') "#Step "
     write(11,'(a)',advance='no') "Y_B(t)X_B(0) Y_AS(t)X_B(0) Y_SS(t)X_B(0) "
     write(11,'(a)',advance='no') "Y_B(t)X_AS(0) Y_AS(t)X_AS(0) Y_SS(t)X_AS(0) "
     write(11,'(a)',advance='yes')"Y_B(t)X_SS(0) Y_AS(t)X_SS(0) Y_SS(t)X_SS(0)"
     do t = -Di,Di
        write(11,'(i8)',advance='no') t
        do i = 1,8
           write(11,'(es14.3)',advance='no') corr(t,i)/surface/dble(nb_shift)
        end do
        write(11,'(es14.3)',advance='yes')corr(t,9)/surface/dble(nb_shift)
     enddo
  elseif(options%inc) then !There is a problem with the layering if we want the cross-correlation
     write(11,*)"#Step    Bond-correlation    Intra-molecular-correlation     Nothing_due_to_layering"
     do t = -Di,Di
        write(11,*)t,corr(t,1)/surface/dble(nb_shift),corr(t,2)/surface/dble(nb_shift),"0"
     enddo
  elseif(options%cross) then !If cross correlation is taken into account
     write(11,*)"#Step    Cross-correlation Nothing Nothing"
     do t = -Di,Di
        write(11,*)t,corr(t,1)/surface/dble(nb_shift),"0 0"
     enddo
  else !General case
     write(11,*)"#Step    Bond-correlation    Intra-molecular-correlation     Full-correlation"
     do t = -Di,Di
        write(11,*)t,corr(t,1)/surface/dble(nb_shift),corr(t,2)/surface/dble(nb_shift),corr(t,3)/surface/dble(nb_shift)
     enddo
  endif
   
  
  !Finishing
  deallocate(vt_bond,vt_mol,vt_all)
  close(11)
  write(*,*)'Correlation results written in file ',trim(options%filename_tmp)

  !Information
  if(options%v0%nb_bond>0)write(*,'(a,f7.3,a)')"Central atoms making molecules for t0: ",&
       100.d0*dble(nb_v0_mol)/dble(nb_shift)/dble(count(options%v0%tab_select)) ,"%"
  if(options%vt%nb_bond>0)write(*,'(a,f7.3,a)')"Central atoms making molecules for t:  ",&
       100.d0*dble(nb_vt_mol)/dble(nb_shift)/dble(count(options%vt%tab_select)) ,"%"
  return
end subroutine vvaf

