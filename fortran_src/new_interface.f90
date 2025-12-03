program FACIT_interface
    implicit none
    integer, parameter :: rkind = kind(1.0d0)

    !===========================================================
    ! PARAMETERS
    !===========================================================
    integer :: nx, nth, nions

    !===========================================================
    ! INPUT ARRAYS
    !===========================================================
    real(rkind), allocatable, dimension(:) :: xn, theta
    real(rkind), allocatable, dimension(:) :: Za
    real(rkind), allocatable, dimension(:) :: Ai
    real(rkind), allocatable, dimension(:) :: Te, Ti, Ne, Na
    real(rkind), allocatable, dimension(:) :: Machi, gradTi, gradNa, deltaM
    real(rkind), allocatable, dimension(:) :: qmag, dpsidx, FV, BV, RV
    real(rkind), allocatable, dimension(:) :: regulopt
    real(rkind), allocatable, dimension(:,:) :: Zi, Ni, gradNi
    real(rkind), allocatable, dimension(:,:) :: jacob, AsymPhi
    real(rkind), allocatable, dimension(:,:,:) :: AsymN, NV

    !===========================================================
    ! OUTPUT ARRAYS FROM FACIT
    !===========================================================
    real(rkind), allocatable, dimension(:) :: Flux_imp, Da, Vconv,dmin,dmaj
    real(rkind), allocatable, dimension(:,:) :: Da_M, Vconv_M
    real(rkind), allocatable, dimension(:,:) :: Da_BP_M, Da_PS_M, Da_CL_M
    real(rkind), allocatable, dimension(:,:) :: Ka_BP_M, Ka_PS_M, Ka_CL_M
    real(rkind), allocatable, dimension(:,:) :: Ha_BP_M, Ha_PS_M, Ha_CL_M
    real(rkind), allocatable, dimension(:,:) :: Va_BP_M, Va_PS_M, Va_CL_M, nn

    real(rkind), allocatable, dimension(:) :: Da_PS, Da_BP, Da_CL
    real(rkind), allocatable, dimension(:) :: Ka_PS, Ka_BP, Ka_CL
    real(rkind), allocatable, dimension(:) :: Ha_PS, Ha_BP, Ha_CL
    real(rkind), allocatable, dimension(:) :: Va_PS, Va_BP, Va_CL

    !===========================================================


    logical :: pol_asym, full_geom, rotation
    integer :: solution
    real(rkind) :: invaspct, B0, R0, Aa
    real(rkind) :: qa, mp,q_e

    real(rkind),allocatable, dimension(:) :: term_ICRH, term_rot, delta_phi, Deltam_phi, anis, epsi
    real(rkind), allocatable, dimension(:) :: Cs, S1, S2, meff, Omega
    integer :: i, j, k

    nx    = 1000
    nth   = 1000
    nions = 5

    ! ALLOCATE ARRAYS

    allocate(xn(nx), theta(nth))
    allocate(Za(nx), Ai(nions))
    allocate(Te(nx), Ti(nx), Ne(nx), Na(nx), Machi(nx), deltaM(nx))
    allocate(gradTi(nx), gradNa(nx))
    allocate(qmag(nx), dpsidx(nx), FV(nx), BV(nx), RV(nx))
    allocate(regulopt(4))
    allocate(Zi(nx,nions), Ni(nx,nions), gradNi(nx,nions))
    allocate(jacob(nx,nth), AsymPhi(nx,nth))
    allocate(AsymN(nx,nth,2),NV(nx,nth,nions))

    !Outputs
    allocate(Flux_imp(nx), Da(nx), Vconv(nx),dmin(nx),dmaj(nx))
    allocate(Da_M(nx,nions), Vconv_M(nx,nions))
    allocate(Da_BP_M(nx,nions), Da_PS_M(nx,nions), Da_CL_M(nx,nions))
    allocate(Ka_BP_M(nx,nions), Ka_PS_M(nx,nions), Ka_CL_M(nx,nions))
    allocate(Ha_BP_M(nx,nions), Ha_PS_M(nx,nions), Ha_CL_M(nx,nions))
    allocate(Va_BP_M(nx,nions), Va_PS_M(nx,nions), Va_CL_M(nx,nions), nn(nx,nth))
    allocate(Da_PS(nx), Da_BP(nx), Da_CL(nx))
    allocate(Ka_PS(nx), Ka_BP(nx), Ka_CL(nx))
    allocate(Ha_PS(nx), Ha_BP(nx), Ha_CL(nx))
    allocate(Va_PS(nx), Va_BP(nx), Va_CL(nx))

    !extras

    allocate(term_ICRH(nx), term_rot(nx), delta_phi(nx), Deltam_phi(nx), anis(nx), epsi(nx))
    allocate(Cs(nx), S1(nx), S2(nx), meff(nx), Omega(nx))


    !===========================================================
    ! GENERATE PROFILES (You fill this section)
    !===========================================================
    invaspct = 0.2_rkind
    B0 = 3.7_rkind
    R0 = 2.5_rkind
    qa = 4.0_rkind
    mp = 1.67262192369e-27_rkind 
    q_e = 1.602176634e-19_rkind 

    do i = 1, nx
        xn(i) = 0.001*(i-1)
    end do

    do i = 1, nth
        theta(i) = (i-1) * 2*(3.14159) / (nth-1)
    end do

    S1= 0.0_rkind
    S2= 0.0_rkind

    do j = 1, nions
        S1 = S1 + Ni(:,j)*Ai(j) *mp 
        S2 = S2 + Ni(:,j)
    end do

    meff = S1/S2

    Cs = sqrt((Te + Ti*S2)/meff)

    Omega = Za * q_e * B0/(Aa *mp )

    Machi = (R0*Omega) / Cs


    Ai = (/2._rkind, 3._rkind, 4._rkind, 14._rkind, 16._rkind/)
    do i=1,nx
        Ni(i,1) = 4e19_rkind*(1-xn(i)**2)**2
        Ni(i,2) = 1e19_rkind*(1-xn(i)**2)**2
        Ni(i,3) = 1e14_rkind*(1-xn(i)**2)**2
        Ni(i,4) = 5e9_rkind*(0.3 + exp(6*(xn(i)-1)))
        Ni(i,5) = 2e10_rkind*(0.3 + exp(6*(xn(i)-1)))
        Na(i) = 6e12_rkind*(0.3 + exp(6*(xn(i)-1)))
        Ne(i) = 6e19_rkind*(1-xn(i)**2)**2

        Zi(i,1)=1; Zi(i,2)=1
        if (xn(i)>=0.6_rkind) then
            Zi(i,3)=1
        else
            Zi(i,3)=2
        end if
        Zi(i,4)= 6 - 5*exp(-((xn(i)-1)/0.25)**2)
        Zi(i,5)= 7 - 6*exp(-((xn(i)-1)/0.25)**2)
        Za(i)= 60 - 55*exp(-((xn(i)-1)/0.25)**2)

        Ti(i) = 1.5_rkind*(1-xn(i)**2)**2
        Te(i) = 3.0_rkind*(1-xn(i)**2)**2

        qmag(i) =1 + (qa-1)*xn(i)**2
        FV(i) = R0*B0
        RV(i) = xn(i) +invaspct*cos(theta(i))
        BV(i) = FV(i)**2 /(RV(i)**2) * (1+ invaspct**2 /(RV(i)**2) * xn(i)**2 /(qmag(i)**2) )

        do j = 1, nth
            jacob(i,j) = xn(i) /(B0*(1.0_rkind-invaspct*cos(theta(j))))
        enddo


        epsi(i) = xn(i) / R0

        
        anis(i) =  (Te(i)/Ti(i) - 1.0_rkind)

        term_ICRH(i) = anis(i)  / ( (Te(i)/Ti(i)))
        term_rot(i)  = Machi(i) * 2.0_rkind            

        delta_phi(i) = epsi(i)/(1.0_rkind + Za(i)*Te(i)/Ti(i)) * (term_ICRH(i) + term_rot(i))

        Deltam_phi(i) = -0.165*epsi(i)         
           
        AsymPhi(i,1) = delta_phi(i)
        AsymPhi(i,2) = Deltam_phi(i)

        do j = 1, nions
            AsymN(i,j,1) = - Za(i) * (Te(i)/Ti(i)) * delta_phi(i)
            AsymN(i,j,2) = - Za(i) * (Te(i)/Ti(i)) * Deltam_phi(i)
            do k = 1, nth
                NV(i,j,k) = 1.0_rkind + AsymN(i,j,1)*cos(theta(k)) + AsymN(i,j,2)*sin(theta(k))
            enddo
        enddo

        deltaM(i) = 2*epsi(i)*meff(i)/(Aa*mp) *Machi(i)**2

    enddo


    !grad Ni

    do j= 1, nions
        gradNi(1,j) = (Ni(2,j) - Ni(1,j)) / (xn(2) - xn(1))
        do i = 2, nx-1
            gradNi(i,j) = (Ni(i+1,j) - Ni(i-1,j)) / (xn(i+1) - xn(i-1))
        end do
        gradNi(nx,j) = (Ni(nx,j) - Ni(nx-1,j)) / (xn(nx) - xn(nx-1))
    enddo

    !Grad Na

    gradNa(1) = (Na(2) - Na(1)) / (xn(2) - xn(1))

    do i = 2, nx-1
        gradNa(i) = (Na(i+1) - Na(i-1)) / (xn(i+1) - xn(i-1))
    end do

    gradNa(nx) = (Na(nx) - Na(nx-1)) / (xn(nx) - xn(nx-1))

    !Grad Ti

    gradTi(1) = (Ti(2) - Ti(1)) / (xn(2) - xn(1))

    do i = 2, nx-1
        gradTi(i) = (Ti(i+1) - Ti(i-1)) / (xn(i+1) - xn(i-1))
    end do

    gradTi(nx) = (Ti(nx) - Ti(nx-1)) / (xn(nx) - xn(nx-1))



    ! Options
    pol_asym  = .true.
    rotation  = .true.
    full_geom = .true.
    solution  = 3
    regulopt  = 0.0_rkind

    !===========================================================
    ! CALL FACIT SUBROUTINE
    !===========================================================
    call FACIT(nx, nth, xn,nions, theta, Za, Aa, Zi, Ai,Te, Ti, Ne, Ni, Na, Machi,gradTi, gradNi, gradNa,invaspct, B0, R0, qmag, dpsidx, FV, BV, RV, jacob, &
               AsymPhi, AsymN, pol_asym, full_geom,solution, rotation, regulopt, Flux_imp, Da, Da_M, Vconv, nn, dmin, dmaj, &
               Da_PS_M, Da_BP_M, Da_CL_M, Ka_PS_M, Ka_BP_M, Ka_CL_M,Ha_PS_M, Ha_BP_M, Ha_CL_M,Va_PS_M,Va_BP_M,Va_CL_M,&
               Da_PS, Da_BP, Da_CL, Ka_PS, Ka_BP, Ka_CL, Ha_PS, Ha_BP, Ha_CL,Va_PS,Va_BP,Va_CL)

    ! supplementary outputs: PS, BP, CL & centrifugal components


    !===========================================================
    ! WRITE INPUTS AND OUTPUTS TO DAT
    !===========================================================
    call write_inputs(nx,xn,Te,Ti,Ne,Ni)
    call write_outputs(nx,nions,Flux_imp,Da,Vconv,Da_M,Vconv_M)

contains

    subroutine write_inputs(nx,xn,Te,Ti,Ne,Ni)
        integer, intent(in) :: nx
        real(rkind), intent(in) :: xn(:), Te(:), Ti(:), Ne(:)
        real(rkind), intent(in) :: Ni(:,:)
        integer :: i
        open(10,file="inputsfa.dat")
        do i=1,nx
            write(10,'(E15.6,1x,E15.6,1x,E15.6,1x,E15.6,1x,E15.6)') xn(i),Te(i),Ti(i),Ne(i),Ni(i,1)
        end do
        close(10)
    end subroutine

    subroutine write_outputs(nx,nions,Flux_imp,Da,Vconv,Da_M,Vconv_M)
        integer, intent(in) :: nx, nions
        real(rkind), intent(in) :: Flux_imp(:), Da(:), Vconv(:)
        real(rkind), intent(in) :: Da_M(:,:), Vconv_M(:,:)
        integer :: i, j

        !=== Vectores 1D ===
        open(11,file="outputs_vectors.dat")
        write(11,*) "# Flux_imp   Da   Vconv"
        do i=1,nx
            write(11,'(E15.6,1x,E15.6,1x,E15.6)') Flux_imp(i), Da(i), Vconv(i)
        end do
        close(11)

        !=== Matriz Da_M ===
        open(12,file="outputs_Da_M.dat")
        write(12,*) "# Da_M(i,j)"
        do j=1,nions
            do i=1,nx
                write(12,'(I5,1x,I3,1x,E15.6)') i,j,Da_M(i,j)
            end do
        end do
        close(12)

        !=== Matriz Vconv_M ===
        open(13,file="outputs_Vconv_M.dat")
        write(13,*) "# Vconv_M(i,j)"
        do j=1,nions
            do i=1,nx
                write(13,'(I5,1x,I3,1x,E15.6)') i,j,Vconv_M(i,j)
            end do
        end do
        close(13)
    end subroutine

end program FACIT_interface
