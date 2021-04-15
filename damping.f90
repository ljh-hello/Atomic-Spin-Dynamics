Module global
    implicit none
    INTEGER, PARAMETER:: DP = selected_real_kind(15,307) 
    integer,parameter :: k15 = selected_int_kind(15)
    Real(DP) :: J1=0.001
    Integer,Parameter :: dim=3
    real(DP):: Bi(dim)
    real(DP) :: damping(3,3,2,2),A11(3,3),A12(3,3),A21(3,3),A22(3,3)
    real(DP) :: mi=1, Gamma=41.39,dt=0.1,damping_onsite=0.1
    !Data Bi /0,0,-4.24726E-5/
    Data Bi /0.d0,0.d0,-4.24726e-05/
    INTEGER,PARAMETER :: row=2,col=3
    INTEGER,PARAMETER :: step=800000
    DATA A12 /9*0.d0/
    DATA A21 /9*0.d0/
    DATA A11 /1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
    DATA A22 /1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/

end Module

module spin_dynamics
    use global
    PRIVATE
    PUBLIC :: spin_update
    contains
    ! function hamitonian(spin) result(energy)
    !     use global
    !     implicit none
    !     Real(DP),dimension(row,col),INTENT(IN) :: spin
    !     real(dp) :: energy
    !     energy = -J1*DOT_PRODUCT(spin(1,:),spin(2,:)) - J1*DOT_PRODUCT(spin(2,:),spin(1,:))
    ! end function

    function find_Beff(s_new,s_old) result(B_eff)
        implicit none
        real(DP),dimension(row,col) :: B_eff,B_eff1,B_eff2,B_eff3
        Real(DP),dimension(row,col),INTENT(IN) :: s_new,s_old
        INTEGER :: ii,jj,kk
        !specify the damping tensor
         damping(:,:,1,1)=A11*damping_onsite  
         damping(:,:,1,2)=A12*damping_onsite 
         damping(:,:,2,1)=A21*damping_onsite 
         damping(:,:,2,2)=A22*damping_onsite  
        !
         DATA B_eff1 /6*0.d0/
         DATA B_eff2 /6*0.d0/
         DATA B_eff3 /6*0.d0/
        do kk=1,row
            if (kk==1) then
                B_eff1(kk,:)= J1*s_new(2,:) + J1*s_new(2,:)
            else
                B_eff1(kk,:)= J1*s_new(1,:) + J1*s_new(1,:)
            end if
        end do
        !effective field from hamitonian and external field (-gamma*B_eff)
        B_eff2=-Gamma*(B_eff1+spread(Bi,1,row)) !hamitonian+Bi
        !!!damping tensor
        do ii=1,row
            do jj=1,row
                B_eff3(ii,:)=B_eff3(ii,:)+(matmul(damping(:,:,ii,jj),((s_new(jj,:)-s_old(jj,:))/dt))/norm2(s_new(jj,:)))
             
            end do
        end do
        !!!damping is a scalar 
        
        ! B_eff3=damping_onsite*((s_new-s_old)/dt)
        B_eff=B_eff2+B_eff3 
        
    end function


    function spin_update(s_new,s_old) RESULT(spin_update1)
        real(DP),DIMENSION(2*row,col) :: spin_update1
        real(DP),dimension(row,col) :: s_new,s_old
        real(DP),dimension(1,row):: old_spin_length,new_spin_length,old_spin_length2,new_spin_length2
        real(DP), DIMENSION(row,col) :: result
        real(DP),dimension(row,col) :: B_eff_new
        real(DP),dimension(2,1) :: normalization,normalization2
        INTEGER :: ii
        DATA result  /6*0.d0/
        old_spin_length2=reshape(norm2(s_old,dim=2),[1,row])
        old_spin_length=reshape(norm2(s_new,dim=2),[1,row])
       
        !effective field
        B_eff_new=find_Beff(s_new,s_old)
        !update m^(t-1)
        s_old=s_new
        new_spin_length2=reshape(norm2(s_old,dim=2),[1,row])
        !update m^(t+1)
        do ii=1,row
            result(ii,:)=cross(s_new(ii,:),B_eff_new(ii,:))
            s_new(ii,:)=s_new(ii,:)+dt*result(ii,:)
        end do
        new_spin_length=reshape(norm2(s_new,dim=2),[1,row])
 
        !Normalization
        normalization=reshape(old_spin_length/new_spin_length,[2,1])
        spin_update1(1,:)=s_new(1,:)*normalization(1,1)
        spin_update1(2,:)=s_new(2,:)*normalization(2,1)
        normalization2=reshape(old_spin_length2/new_spin_length2,[2,1])
        spin_update1(3,:)=s_old(1,:)*normalization2(1,1)
        spin_update1(4,:)=s_old(2,:)*normalization2(2,1)
        
    end function

    !cross product
    FUNCTION cross(a, b) RESULT(c) 
        real(DP), DIMENSION(3) :: c
        real(DP), DIMENSION(3), INTENT(IN) :: a, b
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
      END FUNCTION 
end module

program main
    use global
    use spin_dynamics
    implicit none
    INTEGER t,i
    real(DP),dimension(2,dim) :: spinlattice_old,spinlattice_new,delta
    real(DP),dimension(4,dim) :: spinlattice
    real(DP), allocatable :: mx(:,:),my(:,:),mz(:,:),time(:,:),Mag_z(:,:),Mag_y(:,:),Mag_x(:,:)
    allocate(mx(row,step),my(row,step),mz(row,step),time(step,1),Mag_z(step,1),Mag_y(step,1),Mag_x(step,1))
    !initial guess of m^(-1) and m^(0)
    spinlattice_old=reshape([real(DP):: 0,0,0.01,-0.01,1,1],[row,col])
    delta=reshape([real(DP):: 0,0,0.01,0.01,-0.1,0.1],[row,col])
    !m^(0)  = m(−1) + δm
    spinlattice_new=spinlattice_old+delta
    !main loop
    do t=1,step
            spinlattice=spin_update(spinlattice_new,spinlattice_old)
            spinlattice_new=spinlattice(1:2,:)
            spinlattice_old=spinlattice(3:4,:)
            mx(:,t)=spinlattice_new(:,1)
            my(:,t)=spinlattice_new(:,2)
            mz(:,t)=spinlattice_new(:,3)
            time(t,1)=t*dt/1000
        end do


   Mag_z=reshape((mz(1,:)+mz(2,:))/2,[step,1])
   Mag_x=reshape((mx(1,:)+mx(2,:))/2,[step,1])
   Mag_y=reshape((my(1,:)+my(2,:))/2,[step,1])

   !see the output mx,my,mz in data1.dat
    open(1, file = 'data1.dat')  
       do i=1,step  
          write(1,*) time(i,1), Mag_x(i,1),Mag_y(i,1),Mag_z(i,1) 
       end do  

    close(1) 
end program main


