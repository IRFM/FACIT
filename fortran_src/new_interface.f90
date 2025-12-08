program facit_try
    use facit_mod
    implicit none

    !===========================================================
    ! PARAMETERS
    !===========================================================
    integer :: nx, nth, nions, solution

    !===========================================================
    ! INPUT ARRAYS
    !===========================================================
    real(rkind), allocatable, dimension(:) :: theta
    real(rkind), allocatable, dimension(:) :: Te, Ti, Na, gradTi, Ne, gradNa
    real(rkind), allocatable, dimension(:,:) :: Ni, gradNi
    real(rkind), allocatable, dimension(:) :: qmag, xn, dpsidx, FV, Za
    real(rkind), allocatable, dimension(:) :: Machi2
    real(rkind), allocatable, dimension(:) :: Ai
    real(rkind), allocatable, dimension(:,:) :: Zi
    real(rkind), allocatable, dimension(:,:) :: jacob, BV,RV
    real(rkind), allocatable, dimension(:,:) :: AsymPhi
    real(rkind), allocatable, dimension(:,:,:) :: AsymN
    real(rkind), allocatable, dimension(:) :: regulopt
    logical  :: pol_asym, rotation, full_geom
    real(rkind)   :: B0, R0, invaspct, Aa

    !===========================================================
    ! OUTPUT ARRAYS
    !===========================================================

    real(rkind), allocatable, dimension(:) :: Da, Vconv, dmin, dmaj, Flux_imp
    real(rkind), allocatable, dimension(:,:) :: Da_M, Vconv_M
    real(rkind), allocatable, dimension(:,:) :: nn
    real(rkind), allocatable, dimension(:,:) ::  Da_BP_M, Da_PS_M, Da_CL_M, Ka_BP_M, Ka_PS_M, Ka_CL_M, Ha_BP_M, Ha_PS_M, Ha_CL_M, Va_BP_M, Va_PS_M, Va_CL_M
    real(rkind),  allocatable, dimension(:) :: Da_BP, Da_PS, Da_CL, Ka_BP, Ka_PS, Ka_CL, Ha_BP, Ha_PS, Ha_CL, Va_BP, Va_PS, Va_CL

    !===========================================================
    ! EXTRA ARRAYS
    !===========================================================

    real(rkind), allocatable, dimension(:) :: term_ICRH, term_rot, delta_phi, Deltam_phi, anis, epsi
    real(rkind), allocatable, dimension(:) :: S1, S2, Omega
    real(rkind)   :: ma
    integer :: i, j

    nx      = 100
    nth     = 100
    nions   = 5

    !===========================================================
    ! Allocate
    !===========================================================

    allocate(theta(nth),Te(nx), Ti(nx), Na(nx), gradTi(nx), Ne(nx), gradNa(nx))
    allocate(Ni(nx,nions), gradNi(nx,nions))
    allocate(qmag(nx), xn(nx), dpsidx(nx), FV(nx), Za(nx),Machi2(nx))
    allocate(Ai(nions))
    allocate(Zi(nx,nions))
    allocate(jacob(nx,nth), BV(nx,nth),RV(nx,nth))
    allocate(AsymPhi(nx,nth))
    allocate(AsymN(nx,nions,2))
    allocate(regulopt(4))
    allocate(Da(nx), Vconv(nx), dmin(nx), dmaj(nx), Flux_imp(nx))
    allocate(Da_M(nx,nions), Vconv_M(nx,nions))
    allocate(nn(nx,nth))
    allocate(Da_BP_M(nx,nions), Da_PS_M(nx,nions), Da_CL_M(nx,nions), Ka_BP_M(nx,nions), Ka_PS_M(nx,nions), Ka_CL_M(nx,nions), Ha_BP_M(nx,nions), Ha_PS_M(nx,nions), Ha_CL_M(nx,nions), Va_BP_M(nx,nions), Va_PS_M(nx,nions), Va_CL_M(nx,nions))
    allocate(Da_BP(nx), Da_PS(nx), Da_CL(nx), Ka_BP(nx), Ka_PS(nx), Ka_CL(nx), Ha_BP(nx), Ha_PS(nx), Ha_CL(nx), Va_BP(nx), Va_PS(nx), Va_CL(nx))
    allocate(term_ICRH(nx), term_rot(nx), delta_phi(nx), Deltam_phi(nx), anis(nx), epsi(nx))
    allocate(S1(nx), S2(nx), Omega(nx))

    !===========================================================
    ! Profiles
    !===========================================================

    invaspct = 0.2_rkind  ! Inverse aspect ratio 
    B0       = 3.7_rkind  ! Toroidal magnetic field on axis 
    R0       = 2.5_rkind  ! Major radius
    Aa       = 184.0_rkind ! Atomic mass number of the main impurity
    ma       = Aa*mp  

    do i = 1, nx
        xn(i) = 0.01_rkind*(i-1)
    end do

    do i = 1, nth
        theta(i) = (i-1) * 2.0_rkind*pi / (nth-1)
    end do

    Ai = (/2._rkind, 3._rkind, 4._rkind, 14._rkind, 16._rkind/)

    regulopt = 0.0_rkind

    do i=1,nx
        ! Ion densities
        Ni(i,1) = 4e19_rkind*(1.0_rkind-xn(i)**2)**2  ! Deuterium (example)
        Ni(i,2) = 1e19_rkind*(1.0_rkind-xn(i)**2)**2  ! Tritium (example)
        Ni(i,3) = 1e14_rkind*(1.0_rkind-xn(i)**2)**2  ! Helium (example)
        Ni(i,4) = 5e9_rkind*(0.3_rkind + exp(6.0_rkind*(xn(i)-1.0_rkind))) ! Minor impurity 1 (edge peaked)
        Ni(i,5) = 2e10_rkind*(0.3_rkind + exp(6.0_rkind*(xn(i)-1.0_rkind))) ! Minor impurity 2 (edge peaked)
        Na(i)   = 6e12_rkind*(0.3_rkind + exp(6.0_rkind*(xn(i)-1.0_rkind))) ! Main impurity density (edge peaked)
        Ne(i)   = 6e19_rkind*(1.0_rkind-xn(i)**2)**2  ! Electron density (charge neutrality is not strictly enforced here for simplicity)

        ! Ion charges
        Zi(i,1)=1.0_rkind
        Zi(i,2)=1.0_rkind
        if (xn(i)>=0.6_rkind) then
            Zi(i,3)=1.0_rkind
        else
            Zi(i,3)=2.0_rkind
        end if
        Zi(i,4)= 6.0_rkind - 5.0_rkind*exp(-((xn(i)-1.0_rkind)/0.25_rkind)**2)
        Zi(i,5)= 7.0_rkind - 6.0_rkind*exp(-((xn(i)-1.0_rkind)/0.25_rkind)**2)
        Za(i)= 60.0_rkind - 55.0_rkind*exp(-((xn(i)-1.0_rkind)/0.25_rkind)**2) ! Main impurity charge

        ! Temperatures 
        Ti(i) = 1.5_rkind*(1.0_rkind-xn(i)**2)**2
        Te(i) = 3.0_rkind*(1.0_rkind-xn(i)**2)**2

       
        qmag(i) = 1.0_rkind + (4-1.0_rkind)*xn(i)**2
        FV(i) = R0*B0
    enddo

    do i = 1,nx  
        do j = 1, nth
            RV(i,j) = R0 * (1.0_rkind + xn(i)*invaspct*cos(theta(j))) 
            BV(i,j) = FV(i)**2 / (RV(i,j)**2) * (1.0_rkind + invaspct**2 /(RV(i,j)**2) * xn(i)**2 /(qmag(i)**2) ) 
            jacob(i,j) = xn(i) * R0 / (1.0_rkind - invaspct*cos(theta(j)))
        enddo
        epsi(i) = xn(i) * invaspct 
        anis(i) = (Te(i)/Ti(i) - 1.0_rkind) 
    enddo

    S1(:) = 0.0_rkind 
    S2(:) = Na(:) 

    do j = 1, nions
        S1 = S1 + Ni(:,j)*Ai(j) *mp 
        S2 = S2 + Ni(:,j)
    end do

    ! Impurity rotation frequency (Omega)
    Omega = Za * q_e * B0/(Aa *mp ) 
    
    ! Mach number squared (simplified)
    Machi2 = (ma*Na + S1)* R0**2 *Omega/(S2*2.0_rkind*Ti)

    ! Calculate dpsidx (derivative of poloidal flux with respect to radial coordinate)
    do i = 1, nx
        dpsidx(i) = ( 2.0_rkind * xn(i) * qmag(i) - xn(i)**2 * (2.0_rkind * (4-1.0_rkind) * xn(i)) )  / ( qmag(i)**2 )
    end do

    ! Calculate Asymmetries (AsymPhi and AsymN)
    do i=1,nx 
        term_ICRH(i) = anis(i)  / ( (Te(i)/Ti(i))) 
        term_rot(i)  = sqrt(Machi2(i)) * 2.0_rkind 
        
        delta_phi(i) = epsi(i)/(1.0_rkind + Za(i)*Te(i)/Ti(i)) * (term_ICRH(i) + term_rot(i))
        Deltam_phi(i) = -0.165_rkind*epsi(i)        
    enddo
      
    AsymPhi(:,1) = delta_phi
    AsymPhi(:,2) = Deltam_phi

    ! Density asymmetries (simplified relation to potential)
    do j = 1, nions
    !    AsymN(:,j,1) = - Zi(:,j) * (Te/Ti) * delta_phi
    !    AsymN(:,j,2) = - Zi(:,j) * (Te/Ti) * Deltam_phi
        AsymN(:,j,1) = 1
        AsymN(:,j,2) = 1
    enddo


    !===========================================================
    ! CALCULATE GRADIENTS (Centered finite difference)
    !===========================================================
    ! Grad Ni
    do j= 1, nions
        gradNi(1,j) = (Ni(2,j) - Ni(1,j)) / (xn(2) - xn(1) +1e-33)
        do i = 2, nx-1
            gradNi(i,j) = (Ni(i+1,j) - Ni(i-1,j)) / (xn(i+1) - xn(i-1)+1e-33)
        end do
        gradNi(nx,j) = (Ni(nx,j) - Ni(nx-1,j)) / (xn(nx) - xn(nx-1)+1e-33)
    enddo

    ! Grad Na
    gradNa(1) = (Na(2) - Na(1)) / (xn(2) - xn(1)+1e-33)
    do i = 2, nx-1
        gradNa(i) = (Na(i+1) - Na(i-1)) / (xn(i+1) - xn(i-1)+1e-33)
    end do
    gradNa(nx) = (Na(nx) - Na(nx-1)) / (xn(nx) - xn(nx-1)+1e-33)

    ! Grad Ti
    gradTi(1) = (Ti(2) - Ti(1)) / (xn(2) - xn(1)+1e-33)
    do i = 2, nx-1
        gradTi(i) = (Ti(i+1) - Ti(i-1)) / (xn(i+1) - xn(i-1)+1e-33)
    end do
    gradTi(nx) = (Ti(nx) - Ti(nx-1)) / (xn(nx) - xn(nx-1)+1e-33)


    !===========================================================
    ! SET OPTIONS (These are now module variables)
    !===========================================================
    pol_asym    = .true.
    rotation    = .true.
    full_geom   = .false.
    solution    = 3

    call FACIT(nx, nth, xn,nions, theta, Za, Aa, Zi, Ai, Te, Ti, Ne, Ni, Na, Machi2, &
                    gradTi, gradNi, gradNa, invaspct, B0, R0, qmag, dpsidx, FV, BV, RV, jacob, AsymPhi, AsymN, & 
                    pol_asym, full_geom,solution, rotation, regulopt, Flux_imp, Da, Da_M, Vconv,Vconv_M, nn, dmin, dmaj, & 
                    Da_PS_M, Da_BP_M, Da_CL_M, Ka_PS_M, Ka_BP_M, Ka_CL_M, Ha_PS_M, Ha_BP_M, Ha_CL_M,Va_PS_M,Va_BP_M,Va_CL_M, &
                    Da_PS, Da_BP, Da_CL, Ka_PS, Ka_BP, Ka_CL, Ha_PS, Ha_BP, Ha_CL,Va_PS,Va_BP,Va_CL) 

     
     
    open(20,file='facit_output.dat')
    do i=1,nx
        write(20,203) xn(i), Da(i), Vconv(i), Flux_imp(i), dmin(i), dmaj(i)
    enddo
    close(20)

    203 format(6(E14.7,1X))

    deallocate(theta,Te, Ti, Na, gradTi, Ne, gradNa)
    deallocate(Ni, gradNi)
    deallocate(qmag, xn, dpsidx, FV, Za,Machi2)
    deallocate(Ai)
    deallocate(Zi)
    deallocate(jacob, BV,RV)
    deallocate(AsymPhi)
    deallocate(AsymN)
    deallocate(regulopt)
    deallocate(Da, Vconv, dmin, dmaj, Flux_imp)
    deallocate(Da_M, Vconv_M)
    deallocate(nn)
    deallocate(Da_BP_M, Da_PS_M, Da_CL_M, Ka_BP_M, Ka_PS_M, Ka_CL_M, Ha_BP_M, Ha_PS_M, Ha_CL_M, Va_BP_M, Va_PS_M, Va_CL_M)
    deallocate(Da_BP, Da_PS, Da_CL, Ka_BP, Ka_PS, Ka_CL, Ha_BP, Ha_PS, Ha_CL, Va_BP, Va_PS, Va_CL)
    deallocate(term_ICRH, term_rot, delta_phi, Deltam_phi, anis, epsi)
    deallocate(S1, S2, Omega)

end program facit_try
