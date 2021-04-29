Module global
    implicit none

    INTEGER, PARAMETER:: DP = selected_real_kind(15,307) 
    integer,parameter :: k15 = selected_int_kind(15)

    INTEGER ,               PARAMETER  :: N_J = 1
    Real(DP), DIMENSION(3,3,N_J)       :: J

    INTEGER ,               PARAMETER  :: N_damping = 2
    Real(DP), DIMENSION(3,3,N_damping) :: damping

    INTEGER , DIMENSION(:,:),ALLOCATABLE:: pair2index_damping

    real(DP), DIMENSION(:,:),ALLOCATABLE:: Bext

    ! Set parameters by hand
    !size of supercell
    INTEGER,PARAMETER          ::     Nx=1,Ny=1,Nz=1
    !Read the unit cell file
    CHARACTER(LEN =20) :: posfile='Dimer'
    !Choose a type of damping
    integer(DP)                       ::  d_type=1
    !Set the length of onsite local damping
    real(DP)      ::     damping_onsite=0.1,non_scale=0.1
    !Set steps
    INTEGER,PARAMETER :: step=800000
    INTEGER,PARAMETER :: tmod=500
    integer,parameter :: Niter_midpoint = 100
    real(DP),parameter:: error_midpoint = 1E-8
    !
    real(DP) :: mi=1.d0, Gamma=41.39d0,dt=0.1d0,J_size=0.001d0
    
    
end Module

!      build neighbor list
module find_neighbor_list
    use global
    public :: find_neighbor,read_poscar
    private
    contains
    subroutine read_poscar(latt_const,atom_pos,Natom)
        implicit none
        
        Real(DP), DIMENSION(3,3)       :: latt_const
        integer(DP),PARAMETER ::          fileid=99,latt_row=3
        INTEGER(DP)             ::        error,row
        real(DP), DIMENSION(:,:),ALLOCATABLE:: atom_pos
        INTEGER(DP)             ::        Natom
        open(fileid,file=posfile,status="old",iostat=error)
        if(error/=0) then
            WRITE(*,*) "Fail opening file"
            stop
        end if
        read (fileid,*)
        read (fileid,*)
        DO row = 1,latt_row
            READ(fileid,*) latt_const(row,:)
        END DO
        read (fileid,*)
        read (fileid,*) Natom
        read (fileid,*)
        allocate(atom_pos(Natom,3))
        DO row = 1,Natom
            READ(fileid,*) atom_pos(row,:)
        END DO
        
        atom_pos=matmul(atom_pos,latt_const)
        close(fileid)
    
    end subroutine
  
    subroutine expand_cell(sc_pos)
        implicit NONE
        !expand_cell
        INTEGER(DP)             ::        n,ii,jj,kk,mm
        Real(DP), DIMENSION(3)       :: latt_s
        INTEGER(DP), DIMENSION(3)       :: N_cell,Ncell
        real(DP), DIMENSION(:,:),ALLOCATABLE:: atom_pos,sc_pos
        Real(DP), DIMENSION(3,3)       :: latt_const,latt_const2
        INTEGER(DP)                     ::Natom
        DATA  N_cell /Nx,Ny,Nz/
        
    !expand_cell
        call read_poscar(latt_const,atom_pos,Natom)
        n=1
        do ii=1,3
            latt_const2(ii,:)=latt_const(ii,:)*N_cell(ii)
        end do
        latt_const2=latt_const2*latt_const2
        latt_s=sqrt(sum(latt_const,dim=2))
        allocate(sc_pos(Natom*Nx*Ny*Nz,3))
        do ii=1,Nx
            do jj=1,Ny
                do kk=1,Nz
                    do mm=1,Natom
                        Ncell(1)=ii-1
                        Ncell(2)=jj-1
                        Ncell(3)=kk-1
                        sc_pos(n,:)=matmul(Ncell,latt_const)+atom_pos(mm,:)
                        n=n+1
                    end do
                end do
            end do
        end do
  
    end subroutine
  
!find_neighbor
    subroutine find_neighbor(NN,NL1)
        implicit none

        INTEGER(DP), DIMENSION(:),ALLOCATABLE :: NN
        INTEGER(DP), DIMENSION(:,:),ALLOCATABLE :: NL1
        real(DP)                        ::       dis1,temp,rc,r12
        real(DP), DIMENSION(:),ALLOCATABLE :: N_rc
        INTEGER(DP)             ::        ii,jj
        real(DP), DIMENSION(:,:),ALLOCATABLE:: atom_pos,sc_pos
        Real(DP), DIMENSION(3,3)       :: latt_const
        INTEGER(DP)                          ::  Natom

        call expand_cell(sc_pos)
        call read_poscar(latt_const,atom_pos,Natom)
  
        Allocate(N_rc(Natom*Nx*Ny*Nz-1))
        do ii=2,Natom*Nx*Ny*Nz
            dis1=sqrt(dot_product((sc_pos(1,:)-sc_pos(ii,:)),(sc_pos(1,:)-sc_pos(ii,:))))
            ! print*, dis1
            N_rc(ii-1)=dis1
        end do
        do ii=1,Natom*Nx*Ny*Nz-1
            do jj=ii+1,Natom*Nx*Ny*Nz
                if (N_rc(ii) .lt. N_rc(jj)) then
                    temp = N_rc(ii)
                    N_rc(ii) = N_rc(jj)
                    N_rc(jj) = temp
                end if
            end do
        end do
        Allocate(NN(Natom*Nx*Ny*Nz))
        Allocate(NL1(Natom*Nx*Ny*Nz,Natom*Nx*Ny*Nz))
        NN=0
        NL1=0
        rc=N_rc(Natom*Nx*Ny*Nz-1)+0.1
        do ii=1,Natom*Nx*Ny*Nz
            do jj=1,Natom*Nx*Ny*Nz
                r12=sqrt(dot_product((sc_pos(ii,:)-sc_pos(jj,:)),(sc_pos(ii,:)-sc_pos(jj,:))))
                if ((r12 < rc) .AND. (r12 > 0)) then
                    NN(ii)=NN(ii)+1
                    NL1(ii,NN(ii))=jj
                end if
            end do
        end do
    end subroutine

  end module

  module spin_dynamics
    use global

    PUBLIC :: spin_update
    PUBLIC :: spin_dyn
    PUBLIC :: set_spin
    PUBLIC :: set_calc
    PUBLIC :: hamiltonian

    PRIVATE

    contains
     function hamiltonian(spin,N_atome_tot,NN,NL1) result(energy)
        implicit none
  
        !< input
        Real(DP),dimension(3,N_atome_tot),INTENT(IN)               :: spin
        INTEGER(DP), DIMENSION(N_atome_tot),INTENT(IN)             :: NN
        INTEGER(DP), DIMENSION(N_atome_tot,N_atome_tot),INTENT(IN) :: NL1
        INTEGER(DP),INTENT(IN)                                     :: N_atome_tot
        !< output
        real(dp) :: energy
        !< local
        INTEGER(DP)                                                :: iatom, jatom
        
        energy = 0.d0
        DO iatom = 1, N_atome_tot
          !< effective field from spin-Hamiltonian
          DO jatom = 1, NN(iatom)
              energy = energy - DOT_PRODUCT(spin(:,iatom),MATMUL(J(:,:,1),spin(:,NL1(iatom,jatom))))
          END DO
          energy = energy - DOT_PRODUCT(Bext(:,iatom),spin(:,iatom))
        END DO
        energy = 1.d0/DBLE(N_atome_tot)*energy
      end function
    
      subroutine set_calc(damping_type,N_atome_tot)
        IMPLICIT NONE
        INTEGER(DP)     ::     ii,jj,kk,damping_type
        INTEGER(DP),intent(in)          ::        N_atome_tot
        !< units in Ryd
        !< initialize damping and pair-index,spin-Hamiltonian and pair-index
        
        ALLOCATE(pair2index_damping(N_atome_tot,N_atome_tot))
        ALLOCATE(Bext(3,N_atome_tot))
        do ii=1,N_atome_tot
            Bext(:,ii) =4.24726e-05*(/0.d0,0.d0,-1.d0/)
            do jj=1,N_atome_tot
                if(ii==jj) then
                    pair2index_damping(ii,jj)=1
                else 
                    pair2index_damping(ii,jj)=2
  
                end if
            end do
        end do
        J = 0.d0
        J(1,1,1) = J_size
        J(2,2,1) = J_size
        J(3,3,1) = J_size
        !diagonalonsite damping
        if(damping_type==1) then
          damping = 0.d0
          damping(1,1,1) = damping_onsite
          damping(2,2,1) = damping_onsite
          damping(3,3,1) = damping_onsite
  
      !different onsite damping
        elseif(damping_type==2) then
          damping = 0.d0
          damping(1,1,1) = damping_onsite
          damping(2,2,1) = damping_onsite*non_scale
          damping(3,3,1) = damping_onsite*non_scale*0.9
  
      !full onsite damping
        elseif(damping_type==3) then
          damping(:,:,2) = 0.d0
          damping(:,:,1)=damping_onsite*non_scale
          do kk=1,3
              damping(kk,kk,1)=damping_onsite
          end do
      
      ! diagonal nonlocal damping
        elseif(damping_type==4) then
          damping = 0.d0
          damping(1,1,1) = damping_onsite
          damping(2,2,1) = damping_onsite
          damping(3,3,1) = damping_onsite
          damping(1,1,2) = damping_onsite*non_scale
          damping(2,2,2) = damping_onsite*non_scale
          damping(3,3,2) = damping_onsite*non_scale
        ! full nonlocal damping
        elseif(damping_type==5) then
          damping = 0.d0
          damping(1,1,1) = damping_onsite
          damping(2,2,1) = damping_onsite
          damping(3,3,1) = damping_onsite
          damping(:,:,2) = damping_onsite*non_scale*non_scale
          damping(1,1,2) = damping_onsite*non_scale
          damping(2,2,2) = damping_onsite*non_scale
          damping(3,3,2) = damping_onsite*non_scale
        end if
      end subroutine

      subroutine find_Beff(t, s_new, delta, B_eff, B_damp,N_atome_tot,NN,NL1)
        implicit none

        !< input
        integer(DP)                       :: t
        real(DP),dimension(3,N_atome_tot) :: delta
        real(DP), DIMENSION(3,N_atome_tot):: s_new
        INTEGER(DP),intent(in)            ::        N_atome_tot
        INTEGER(DP), DIMENSION(N_atome_tot):: NN
        INTEGER(DP), DIMENSION(N_atome_tot,N_atome_tot) :: NL1
        !< output
        real(DP),dimension(3,N_atome_tot) :: B_eff
        real(DP),dimension(3,N_atome_tot) :: B_damp
        !< local
        REAL(DP),DIMENSION(3)       :: B_eff1,B_eff2,B_eff3, dmdt, mxB
        INTEGER(DP) :: iatom, jatom, jp,nu
        
        
        
        DO iatom = 1, N_atome_tot
          !< effective field from spin-Hamiltonian
          B_eff1 = 0.d0
          DO jatom = 1, NN(iatom)
              DO nu=1,3,1
                B_eff1(nu) = B_eff1(nu) + DOT_PRODUCT(J(nu,1:3,1),s_new(1:3,NL1(iatom,jatom)))+ &
                DOT_PRODUCT(s_new(1:3,NL1(iatom,jatom)),J(1:3,nu,1))
              END DO
         
          END DO

          !< effective field from external field
          B_eff2 = Bext(:,iatom)

          B_eff(:,iatom)=-gamma * (B_eff1 + B_eff2)
          IF(MOD(t,tmod) == 0) THEN
            WRITE(12,"(ES20.8,I8,6ES20.8)",ADVANCE="YES") t*dt, iatom, B_eff1, B_eff2
          END IF
        END DO

        DO iatom = 1, N_atome_tot
          !< effective field from dissipation
          B_eff3 = 0.d0
          DO jatom = 1, N_atome_tot
            jp = pair2index_damping(iatom,jatom)
            IF( jp /= 0) THEN
              mxB=cross(s_new(:,jatom),B_eff(:,jatom))
              dmdt = (delta(:,jatom))/dt

              DO nu=1,3,1
!                B_eff3(nu) = B_eff3(nu) + DOT_PRODUCT(damping(nu,1:3,jp),mxB(1:3))
                B_eff3(nu) = B_eff3(nu) + DOT_PRODUCT(damping(nu,1:3,jp),dmdt(1:3))
!                B_eff3(nu) = B_eff3(nu) + DOT_PRODUCT(damping(nu,1:3,jp),dmdt(1:3))+DOT_PRODUCT(dmdt(1:3),damping(1:3,nu,jp))
!                 WRITE(99,*) iatom, jatom, DOT_PRODUCT(damping(nu,1:3,jp),mxB(1:3))
              END DO

            END IF
          END DO

          B_damp(:,iatom)=B_eff3
          IF(MOD(t,tmod) == 0) THEN
            WRITE(13,"(ES20.8,I8,3ES20.8)",ADVANCE="YES") t*dt/1000, iatom, 1.d0/gamma*B_eff3
          END IF
        END DO

!STOP
       end subroutine


       subroutine spin_dyn(t,s_new,s_update,N_atome_tot,NN,NL1)
        implicit none

        !-- input
        integer(DP)                         :: t
        real(DP), DIMENSION(3,N_atome_tot) :: s_new
        INTEGER(DP),intent(in)            ::        N_atome_tot
        INTEGER(DP), DIMENSION(N_atome_tot),intent(in) :: NN
        INTEGER(DP), DIMENSION(N_atome_tot,N_atome_tot),intent(in) :: NL1
        !-- output
        real(DP),DIMENSION(3,N_atome_tot) :: s_update
        !-- local
        real(DP),dimension(3,N_atome_tot) :: B_eff, B_damp, delta, tmp_delta,s_mid
        real(DP),dimension(3) :: mxB, tmp_s
        real(DP)              :: norm_s
        INTEGER(DP)   :: iatom, iiter
        logical :: converged
        
        delta=0.d0
        tmp_delta=0.d0

        DO iiter=1,Niter_midpoint,1
          converged=.true.
          DO iatom=1,N_atome_tot
            tmp_s  = (delta(:,iatom)+2.d0*s_new(:,iatom))/2.d0             !s_new fixed
            norm_s = norm2(tmp_s)
            s_mid(:,iatom) = mi/norm_s*tmp_s
          END DO

          !effective field
          CALL find_Beff(t, s_mid, delta, B_eff, B_damp,N_atome_tot,NN,NL1)

          !update m^(t+1)
          do iatom=1,N_atome_tot
            mxB(:)= cross(s_mid(:,iatom),B_eff(:,iatom)+B_damp(:,iatom))
            delta(:,iatom) = dt*mxB(:)
          end do

          do iatom=1,N_atome_tot
            if(.NOT.norm2(tmp_delta(:,iatom)-delta(:,iatom))<error_midpoint) then
              converged=.false.
              tmp_delta = delta
              exit
            end if
          end do
          if(converged) exit
        end do
        if(iiter == Niter_midpoint) WRITE(6,*) "warrning!!!!"

        DO iatom=1,N_atome_tot
          tmp_s  = delta(:,iatom)+s_new(:,iatom)
          norm_s = norm2(tmp_s)
          s_update(:,iatom) = mi/norm_s*tmp_s
        END DO

       end subroutine

       subroutine spin_update(s_new,s_update,N_atome_tot)
        implicit none
        !-- input
        INTEGER(DP)                         :: N_atome_tot
        real(DP), DIMENSION(3,N_atome_tot)  :: s_new
        !-- output
        real(DP), DIMENSION(3,N_atome_tot)  :: s_update
        
        s_new = s_update

       end subroutine

       FUNCTION cross(a, b) RESULT(c) 
        real(DP), DIMENSION(3) :: c
        real(DP), DIMENSION(3), INTENT(IN) :: a, b
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
       END FUNCTION 
       
       SUBROUTINE set_spin(s_new,N_atome_tot)
        IMPLICIT NONE
        !-- input
        real(DP), DIMENSION(3,N_atome_tot) :: s_new
        !-- local
        real(DP) :: norm_s
        integer(DP)                                :: ii
        INTEGER(DP),intent(in)             ::        N_atome_tot
        !initial guess of m^(-1) and m^(0)
        do ii=1,N_atome_tot
          s_new(:,ii)=((/0.d0, (0.001d0*ii),1.d0/))
          norm_s = norm2(s_new(:,ii))
          s_new(:,ii) = mi/norm_s * s_new(:,ii)
        end do
  
  
       END SUBROUTINE
  end module



program main
    use global
    use find_neighbor_list, ONLY: read_poscar,find_neighbor
    use spin_dynamics, ONLY: set_calc,spin_dyn, spin_update, set_spin, hamiltonian
    
    implicit none
    INTEGER(DP)                       :: t,iatom
    real(DP),dimension(3)       :: M
    real(DP)                    :: time, E

    Real(DP), DIMENSION(3,3)       :: latt_const
    real(DP), DIMENSION(:,:),ALLOCATABLE:: atom_pos,s_new, s_update
    INTEGER(DP), DIMENSION(:),ALLOCATABLE :: NN
    INTEGER(DP), DIMENSION(:,:),ALLOCATABLE :: NL1
    INTEGER(DP)                           ::  Natom,N_atome_tot



    call read_poscar(latt_const,atom_pos,Natom)
    call find_neighbor(NN,NL1)
    !< set relevant interaction parameters

    !total number of atome in supercell
    N_atome_tot=Natom*Nx*Ny*Nz
    ALLOCATE(s_new(3,N_atome_tot))
    ALLOCATE(s_update(3,N_atome_tot))
    CALL set_spin(s_new,N_atome_tot)
    CALL set_calc(d_type,N_atome_tot)

    !< open output file
    open(10, file = 'moment_tot.dat')  
    write(10,"(6A20)") "#_time (fs)", "E", "|M|", "M%x", "M%y", "M%z"
    open(11, file = 'moment.dat')  
    write(11,"(A20,A8,4A20)") "#_time (fs)", "site", "|m|", "m%x", "m%y", "m%z"
    open(12, file = 'beff.dat')  
    write(12,"(A20,A8,6A20)") "#_time (fs)", "site", "Bexch%x", "Bexch%y", "Bexch%z", "Bext%x", "Bext%y", "Bext%z"
    open(13, file = 'bdamp.dat')  
    write(13,"(A20,A8,3A20)") "#_time (fs)", "site", "Bdamp%x", "Bdamp%y", "Bdamp%z"

    ! main loop
    time=0.d0
    do t=1,step
        !< do one spin dynamics step
        CALL spin_dyn(t,s_new,s_update,N_atome_tot,NN,NL1)
  
        !< determine the average magnetic moment
        M = 1.d0/DBLE(N_atome_tot)*SUM(s_update(:,:),dim=2)
  
        !< total energy density
        E = hamiltonian(s_update,N_atome_tot,NN,NL1)
  
        !< print the average magnetic moment
        IF(MOD(t,tmod) == 0) THEN
        write(10,"(7ES20.8)") time, E, norm2(M),M, DOT_PRODUCT(s_update(:,1),s_update(:,2))
        do iatom=1,N_atome_tot
          write(11,"(ES20.8,I8,4ES20.8)") time, iatom, norm2(s_update(:,iatom)),s_update(:,iatom)
        end do
        END IF
  
        !< update spin arrays
        CALL spin_update(s_new,s_update,N_atome_tot)
  
        !< increase tims step
        time=t*dt/1000
      end do

      close(10)
      close(11)
      close(12)
      close(13)

end program main
