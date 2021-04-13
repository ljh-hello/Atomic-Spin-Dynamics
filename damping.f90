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

    function find_Beff(spin) result(B_eff )
        implicit none
        real(DP),dimension(row,col) :: B_eff,B_eff1,B_eff2,B_eff3
        Real(DP),dimension(row,col),INTENT(IN) :: spin
        INTEGER :: ii,jj,kk
         damping(:,:,1,1)=A11*damping_onsite  
         damping(:,:,1,2)=A12*damping_onsite 
         damping(:,:,2,1)=A21*damping_onsite 
         damping(:,:,2,2)=A22*damping_onsite  

        do kk=1,row
            if (kk==1) then
                B_eff1(kk,:)= J1*spin(2,:) + J1*spin(2,:)
            else
                B_eff1(kk,:)= J1*spin(1,:) + J1*spin(1,:)
            end if
        end do
        B_eff2=(B_eff1+spread(Bi,1,row)) !hamitonian+Bi
        ! do ii=1,row
        !     do jj=1,row
        !     B_eff3(ii,:)=B_eff3(ii,:)+(matmul(damping(:,:,ii,jj),cross(spin(jj,:),B_eff2(jj,:)))/norm2(spin(jj,:)))
        !     end do
        ! end do
        do ii=1,row
            B_eff3(ii,:)=damping_onsite*cross(spin(ii,:),B_eff2(ii,:))
        end do
        B_eff=B_eff2+B_eff3 
    end function

    function derivate(spin_single,B_single) RESULT(deriv)
        real(DP), DIMENSION(3) :: deriv
        real(DP),dimension(3),intent(in) :: spin_single
        real(DP),dimension(3),intent(in) :: B_single
        deriv= -1*Gamma*cross(spin_single,B_single)
        end function

    function spin_update(spin,n) RESULT(spin_update1)
        real(DP),DIMENSION(row,col) :: spin_update1
        real(DP),dimension(row,col) :: spin
        real(DP),dimension(1,row):: old_spin_length,new_spin_length
        real(DP), DIMENSION(3) :: result
        real(DP),dimension(row,col) :: B_eff_new
        old_spin_length=reshape(norm2(spin,dim=2),[1,row])
        B_eff_new=find_Beff(spin)
        result=derivate(spin(n,:),B_eff_new(n,:))
        spin_update1(n,1)=spin(n,1)+dt*result(1) 
        spin_update1(n,2)=spin(n,2)+dt*result(2) 
        spin_update1(n,3)=spin(n,3)+dt*result(3) 
        new_spin_length=reshape(norm2(spin_update1,dim=2),[1,row])
    end function
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
    INTEGER t,s,i
    real(DP),dimension(2,dim) :: spinlattice
    real(DP), allocatable :: mx(:,:),my(:,:),mz(:,:),time(:,:),Mag_z(:,:),Mag_y(:,:),Mag_x(:,:)
    allocate(mx(row,step),my(row,step),mz(row,step),time(step,1),Mag_z(step,1),Mag_y(step,1),Mag_x(step,1))
    spinlattice=reshape([real(DP):: 0,0,0.01,-0.01,1,1],[row,col])
    
    do t=1,step
        do s=1 ,row
            spinlattice=spin_update(spinlattice,s)
            mx(s,t)=spinlattice(s,1)
            my(s,t)=spinlattice(s,2)
            mz(s,t)=spinlattice(s,3)
            time(t,1)=t*dt/1000
        end do
    end do

   Mag_z=reshape(mz(1,:)+mz(2,:),[step,1])
   Mag_x=reshape(mx(1,:)+mx(2,:),[step,1])
   Mag_y=reshape(my(1,:)+my(2,:),[step,1])
    open(1, file = 'data1.dat')  
       do i=1,step  
          write(1,*) time(i,1), Mag_x(i,1),Mag_y(i,1),Mag_z(i,1) 
       end do  

    close(1) 
end program main


! command+shift+/是注释
