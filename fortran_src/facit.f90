! Software name : FACIT
! Authors : P. Maget and D. Fajardo, C. Angioni, P. Manas
! Copyright holders : Commissariat à l’Energie Atomique et aux Energies Alternatives (CEA), France, Max-Planck Institut für Plasmaphysik, Germany
! CEA and IPP authorize the use of the FACIT software under the CeCILL-C open source license https://cecill.info/licences/Licence_CeCILL-C_V1-en.html  
! The terms and conditions of the CeCILL-C license are deemed to be accepted upon downloading the software and/or exercising any of the rights granted under the CeCILL-C license.
! 
subroutine FACIT(nx, nth, xn,nions, theta, &                                ! grid parameters
                 Za, Aa, Zi, Ai, &                                    ! impurity and main ion charge and mass
                 Te, Ti, Ne, Ni, Na, Machi, &                         ! plasma profiles
                 gradTi, gradNi, gradNa, &                            ! gradients
                 invaspct, B0, R0, qmag, dpsidx, FV, BV, RV, jacob, & ! equilibrium
                 AsymPhi, AsymN, &                                    ! input asymmetries
                 pol_asym, full_geom, rotation, regulopt, &           ! different options in the model
                 Flux_imp, Da, Vconv, nn, dmin, dmaj, &               ! main outputs: flux, transport coefficients and poloidal asymmetry
                 Da_PS, Da_BP, Da_CL, Ka_PS, Ka_BP, Ka_CL, &          ! supplementary outputs: PS, BP, CL & centrifugal components
                 Ha_PS, Ha_BP, Ha_CL,Va_PS,Va_BP,Va_CL)               ! supplementary outputs: PS, BP, CL

!*******************************************************************************
! Self-consistent calculation of collisional impurity transport and poloidal
! distribution of the impurity density
!------------
! References:
! Plasma Phys. Control. Fusion (2020) 62 105001 (doi: 10.1088/1361-6587/aba7f9)
! Modifications : 
! * 2022/01/10 : 1. correction to rotation term (P. Maget, corrigendum submitted to PPCF)
! *              2. last three outputs are Va_PS,Va_BP,Va_CL (one more than before)
! *              3. Use of a metric qmag_metrics consistent with full equilibrium when dpsidx available
! * 2022/01/21 : 1. correction of Cgeo_G for full_geom=.t. (line 288)
! *              2. correction of GGG for full_geom=.t. (line 778)
!* 2025/07/04 : 1. qmag input is the good quantity, compute dpsidx from it when not given in input
!*                          2. Correct mu_ie expression to be as in D. Fajardo's paper
!*******************************************************************************
! INPUTS:
! -----------> description [unit] {variable type, shape}
! ------------
! nx --------> size of radial arrays [-] {int}
! nth -------> size of poloidal arrays [-] {int}
!nions --------> number of ion species [-] {int}
! xn --------> radial coordinate [-] {arr, nx}
! theta -----> poloidal coordinate [-] {arr, nth}
! Za --------> impurity charge number [-] {arr, nx}
! Aa --------> impurity mass number [-] {float}
! Zi --------> main ion charge number [-] {array, (nx, nions)}
! Ai --------> main ion mass number [-] {array, nions}
! Te --------> electron temperature [eV] {arr, nx}
! Ti --------> main ion temperature [eV] {arr, nx}
! Ne --------> electron density [1/m^3] {arr, nx}
! Ni --------> main ion denisty [1/m^3] {arr, (nx,nions)}
! Na --------> impurity density [1/m^3] {arr, nx}
! Machi -----> Mach number of main ion [-] {arr, nx} check
! gradTi ----> main ion temperature gradient [eV/-] {arr, nx}
! gradNi ----> main ion density gradient [1/m^3/-] {arr, (nx,nions)}
! gradNa ----> impurity density gradient [1/m^3/-] {arr, nx}
! invaspct --> inverse aspect ratio [-] {float}
! B0 --------> magnetic field at magnetic axis [T] {float}
! R0 --------> major radius at magnetic axis [m] {float}
! qmag ------> safety factor [-] {arr, nx}
! dpsidx ----> radial derivative of poloidal flux [V*s/-] {arr, nx}
! FV --------> poloidal current flux function [T*m] {arr, nx}
! BV --------> magnetic field [T] {arr, (nx,nth)}
! RV --------> major radius [m] {arr, (nx,nth)}
! jacob -----> Jacobian of the coordinate system [m/T] {arr, (nx,nth)}
! AsymPhi ---> poloidal asymmetry of electrostatic potential [-] {arr, (nx,2)}
! AsymN -----> poloidal asymmetry of main ion density [-] {arr, (nx,2)}
! pol_asym --> poloidally symmetric (False) or asymmetric (True) system {boolean}
! full_geom -> analytical (False) or iterative (True) geometry {boolean}
! rotation --> consider toroidal rotation {boolean}
! regulopt --> options for iterative calculations [-] {arr, 4}
!*******************************************************************************
! OUTPUTS:
! --------
! Flux_imp --> surface-averaged impurity flux [1/(m^2 s)] {arr, nx}
! Da --------> diffusion coefficient [m^2/s] {arr, nx}
! Vconv -----> convective velocity [m/s] {arr, nx}
! nn --------> poloidal asymmetry of the impurity density Na/<Na> [-] {arr, (nx,nth)}
! dmin ------> horizontal asymmetry of the impurity density [-] {arr, nx}
! dmaj ------> vertical asymmetry of the impurity density [-] {arr, nx}
! Da_* ------> Pfirsch-Schlüter, Banana-Plateau and classical components of
!              the diffusion coefficient [m^2/s] {arr, nx}
! Ka_* ------> components of the coefficient of the main ion density gradient [m^2/s] {arr, nx}
! Ha_* ------> components of the coefficient of the main ion temperature gradient [m^2/s] {arr, nx}
! Va_* ------> Pfirsch-Schlüter, Banana-Plateau and classical components of
!              the pinch velocity [m/s] {arr, nx}!*******************************************************************************

  use constants, only: rkind, mp, me, q_e, eps_pi_fac, sqrt2
  implicit none

  !------------------------------ declarations --------------------------------

  ! INPUTS
  integer :: nx, nth, nions
  real(rkind), dimension(nth) :: theta
  real(rkind), dimension(nx)  :: Te, Ti, gradTi, Ne, gradNa
  real(rkind), dimension(nx,nions)  :: Ni, gradNi
  real(rkind), dimension(nx)  :: qmag, xn, dpsidx, FV, Za
  real(rkind), dimension(nx)  :: Machi 
  real(rkind) :: B0, R0, invaspct, Aa
  real(rkind), dimension(nions)  :: Ai
  real(rkind), dimension(nx,nions)  :: Zi
  real(rkind), dimension(nx,nth) :: jacob, BV, RV, PhiV, NV
  real(rkind), dimension(nx,2) :: AsymPhi, AsymN
  real(rkind), dimension(4) :: regulopt
  logical :: pol_asym, full_geom, rotation

  ! OUTPUTS
  real(rkind), dimension(nx,nions) :: Da, Vconv, Flux_imp, dmin, dmaj
  real(rkind), dimension(nx) ::  Das, dmins, dmajs, Vconvs, Flux_imps
  real(rkind), dimension(nx, nth) :: nn
  real(rkind), dimension(nx,nions) :: Da_BP, Da_PS, Da_CL, Ka_BP, Ka_PS, Ka_CL, Ha_BP, Ha_PS, Ha_CL, Va_BP, Va_PS, Va_CL
  real(rkind), dimension(nx) :: Da_BPs, Da_PSs, Da_CLs, Ka_BPs, Ka_PSs, Ka_CLs, Ha_BPs, Ha_PSs, Ha_CLs, Va_BPs, Va_PSs, Va_CLs

  ! OTHER
  integer :: i,j, p, n, it, ix!, info, ierr, ierrmax
  real(rkind) :: amin, ma, ftrap, ki_Redl, C2
  real(rkind), dimension(nions)  :: mi
  real(rkind), dimension(nx) ::  grad_ln_na, Te15, Ti15, grad_ln_Ti
  real(rkind), dimension(nx,nions)  :: grad_ln_ni
  real(rkind), dimension(nx) :: epsk, eps15, eps2, ft, wca, dD2, deltaM
  real(rkind), dimension(nx,nions) :: ki, C0a, g  
  real(rkind), dimension(nx) :: f1, f2, f3, y1, y2, y3, y4, adps
  real(rkind), dimension(nx) :: LneeNRL, LneimpNRL
  real(rkind), dimension(nx,nions) :: LneiNRL
  real(rkind), dimension(nx) :: LnimpeNRL, LnimpimpNRL
  real(rkind), dimension(nx,nions) :: LniiNRL, LniimpNRL
  real(rkind), dimension(nx) :: Tauee, Taueimp
  real(rkind), dimension(nx,nions) :: Tauei
  real(rkind), dimension(nx,nions) :: Tauie, Tauii, Tauiimp
  real(rkind), dimension(nx) :: Tauimpe, Tauimpimp
  real(rkind), dimension(nx,nions) :: Tauimpi
  real(rkind), dimension(nx) :: wee, wimpimp, nuestar, nuimpstar, rhoLimp2
  real(rkind), dimension(nx,nions) :: wii, nuistar
  real(rkind), dimension(nx, nions) :: L11impi, nuswca, mu_ie, alpha
  real(rkind), dimension(nx) ::  B2avg, dNH, dNV, dminphia, dmajphia, Zeff
  real(rkind), dimension(nx,nions) :: UU, GG
  real(rkind), dimension(nx) :: Cgeo_G, Cgeo_U, Cgeo_Gcl
  real(rkind), dimension(nx, nth) :: b2
  real(rkind), dimension(nx) :: b2navg, nb2avg, nNVavg, b2NVavg
  real(rkind), dimension(nx) :: K11a, K12a, K22a
  real(rkind), dimension(nx,nions) :: K11i, K12i, K22i
  real(rkind), dimension(nx,nions) :: Vra_BP, Vra_PS, Vra_CL
  real(rkind), dimension(nx,nions) :: Ka, Ha

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI
  !(previous two lines here because "A specification statement cannot appear in the executable section")

  !---------------------------------------------------------------------------

  mi = Ai*mp ! main ion mass
  ma = Aa*mp ! impurity mass

  epsk  = xn*invaspct ! local inverse aspect ratio
  eps15 = (epsk + 1.0e-33)**1.5
  eps2  = (epsk + 1.0e-33)**2
  amin  = invaspct*R0 ! minor radius

  if (maxval(abs(dpsidx)).lt.1.0e-33_rkind) then
	dpsidx = amin**2*B0*xn/qmag
  endif

Zeff =Za**2*Na/(Ne)
do j = 1,nions
  Zeff += Zi(j)**2*Ni(:,j)/(Ne)
enddo

  do i = 1, nx 
    grad_ln_na(i) = gradNa(i)/(Na(i) + 1.e-33)
    grad_ln_Ti(i) = gradTi(i)/(Ti(i) + 1.e-33)
    ! trapped particle fraction
    ft(i) = ftrap(epsk(i))
    do j = 1, nions
    ! logarithmic gradients
      grad_ln_ni(i,j) = gradNi(i,j)/(Ni(i,j) + 1.e-33)
    enddo  
  enddo


  if (.not.rotation) then
    Machi = 0.0_rkind  !check with Patrick
  endif

  deltaM = 2*(Aa/Ai)*Machi**2*epsk ! rotation strength parameter


  ! Coulomb Logarithms (from NRL formulary)
  LneeNRL = 23.5 - 0.5*log(Ne/1e6) + 1.25*log(Te) - sqrt(1e-5 + (1.0/16.0)*(log(Te) - 2)**2)

  do p = 1, nx
    do j = 1, nions
      if ((Ti(p)*me/(mi(j))<Te(p)).and.(Te(p)<10*Zi(p)**2)) then
        LneiNRL(p,j) = 23 - 0.5*log(Zi**2*Ne(p)/1e6) + 1.5*log(Te(p))
      elseif ((Ti(p)*me/(mi(j))<10*Zi**2).and.(Te(p)>10*Zi(p)**2)) then
        LneiNRL(p,j) = 24 - 0.5*log(Ne(p)/1e6) + log(Te(p))
      else
        LneiNRL(p,j) = 30 - log(Zi(p)**2/Ai(p)*(Ni(p,j)/1e6)**0.5) + 1.5*log(Ti(p))
      endif
    enddo 

    if ((Ti(p)*me/(Aa*mp)<Te(p)).and.(Te(p)<10*Za(p)**2)) then
      LneimpNRL(p) = 23. - log(Za(p)*(Ne(p)/1e6)**0.5) + 1.5*log(Te(p))
    elseif ((Ti(p)*me/(Aa*mp)<10.*Za(p)**2).and.(Te(p)>10.*Za(p)**2)) then
      LneimpNRL(p) = 24. - 0.5*log(Ne(p)/1e6) + log(Te(p))
    else
      LneimpNRL(p) = 30. - log(Za(p)**2/Aa*(Na(p)/1e6)**0.5) + 1.5*log(Ti(p))
    endif

  enddo

  do j = 1,nions
    LniiNRL(:,j)     = 23. - log(Zi(j)*Zi(j)*sqrt(2*(Ni(:,j)/1e6)*Zi(j)**2)) + 1.5*log(Ti)
    LniimpNRL(;,j)   = 23. - log(Zi(j)*Za*sqrt((Ni(:,j)/1e6)*Zi(j)**2+(Na/1e6)*Za**2)) + 1.5*log(Ti)
  enddo 
  LnimpimpNRL = 23. - log(Za*Za*sqrt((Na/1e6)*Za**2+(Na/1e6)*Za**2)) + 1.5*log(Ti)


  ! Collision times (Braginskii)
  Ti15 = Ti**1.5
  Te15 = Te**1.5

  Tauee     = (eps_pi_fac*sqrt(me)*Te15)/(Ne*LneeNRL)
  do j=1,nions
    Tauei(:,j)     = (eps_pi_fac*sqrt(me)*Te15)/(Zi**2*Ni(:,j)*LneiNRL(:,j))
  enddo
  Taueimp   = (eps_pi_fac*sqrt(me)*Te15)/(Zi**2*Za**2*Na*LneimpNRL) !Why Zi

  do j =1,nions
    Tauie(:,j)     = (eps_pi_fac*sqrt(mi(j))*Ti15)/(Zi(j)**2*Ne*LneiNRL(:,j))
    Tauii(:,j)     = (eps_pi_fac*sqrt(mi(j))*Ti15)/(Zi(j)**4*Ni(:,j)*LniiNRL(:,j))
    Tauiimp(:,j)   = (eps_pi_fac*sqrt(mi(j))*Ti15)/(Zi(j)**2*Za**2*Na*LniimpNRL(:,j))
  enddo

  Tauimpe   = (eps_pi_fac*sqrt(ma)*Ti15)/(Za**2*Ne*LneimpNRL)
  do j =1, nions
    Tauimpi(:,j)   = (eps_pi_fac*sqrt(ma)*Ti15)/(Zi(j)**2*Za**2*Ni(:,j)*LniimpNRL(:,j))
  enddo
  Tauimpimp = (eps_pi_fac*sqrt(ma)*Ti15)/(Za**4*Na*LnimpimpNRL)


  ! impurity collision frequency
  do j = 1,nions
    L11impi(j) = 1.0/(sqrt(1.0 + Aa/Ai(j))*Tauimpi(:,j))
  enddo

  ! transit frequencies
  wee     = (2.0*q_e*Te/me)**0.5/(R0*qmag)
  do j = 1,nions
    wii(:,j)     = (2.0*q_e*Ti/mi(j))**0.5/(R0*qmag)
  enddo
  wimpimp = (2.0*q_e*Ti/ma)**0.5/(R0*qmag)

  ! collisionalities
  nuestar   = 1.0/((eps15 + 1.e-33)*wee*Tauee)
  do j = 1,nions
    nuistar(:,j)   = 1.0/((eps15 + 1.e-33)*wii(:,j)*Tauii(:,j))
  enddo
  nuimpstar = 1.0/((eps15 + 1.e-33)*wimpimp*Tauimpimp)

  wca    = q_e*Za*B0/ma ! impurity cyclotron frequency
  do j=1,nions
    nuswca = L11impi(:,j)/wca ! ratio of impurity collision frequency to cyclotron frequency
  enddo

  ! fitted factors
  call facs(nx, Za, ft, f1, f2, f3, y1, y2, y3, y4, adps)
  do j = 1,nions
    g(:,j) = nuistar(:,j)*eps15 ! collisionality parameter
    mu_ie(:,j) = (96.0*sqrt2/125.0)*(1/Zi(j)**2)*sqrt(me/mi(j))*(Ti15/Te15) ! ion-electron heat exchange term (Fülöp-Helander PoP '01)
    alpha(:,j) = Na*Za**2/(Ni(:,j)*Zi(j)**2) ! impurity strength parameter
  enddo

  do i = 1, nx
    do j=1,nions
      C0a(i,j) = C2(alpha(i,j), g(i,j), f1(i), f2(i), Aa, Ai(j))/(1 + f3(i)*mu_ie(i,j)*g(i,j)**2) ! coefficient of ion heat flux in impurity-ion friction
      ki(i,j)  = ki_Redl(nuistar(i,j), ft(i), Zeff(i)) ! neoclassical ion flow coefficient
    enddo
  enddo

  rhoLimp2 = (2.0*q_e*Ti/ma)/wca**2 ! Impurity Larmor radius (squared)

  ! thermodynamic gradients (for asymmetry calculations)
  do j =1,nions
    UU(:,j)  = -(Za/Zi(j))*(C0a(:,j) + ki(:,j))*grad_ln_Ti
    GG(:,j)  = grad_ln_na - (Za/Zi(j))*grad_ln_ni(:,j) + (1 + (Za/Zi(j))*(C0a(:,j) - 1))*grad_ln_Ti
  enddo


  !---------------------------------------------------------------------------
  !-----------------------  Poloidal asymmetry  ------------------------------
  !---------------------------------------------------------------------------


  if (full_geom) then

    ! Full geometry calculation

    call fluxavg(nx, nth, theta, BV**2, jacob, B2avg)

    do ix = 1, nx
      b2(ix,:) = BV(ix,:)**2/B2avg(ix)
    enddo


    if (pol_asym) then

      ! Poloidally asymmetric case


      do i = 1, nx
        PhiV(i,:) = AsymPhi(i,1)*cos(theta) + AsymPhi(i,2)*sin(theta)
       	NV(i,:) = 1. + AsymN(i,1)*cos(theta) + AsymN(i,2)*sin(theta)
      enddo

      call asymmetry_fg(nx, nth, theta, BV, RV, jacob, FV, dpsidx, Machi, L11impi, &
                        R0, Ai, Aa, Zi, Za, B2avg, UU, GG, PhiV, NV, Te, Ti, regulopt, &
                        dmin, dmaj, nn, dmins, dmajs)

    else

      ! Poloidally symmetric case

      dmins = 0.0_rkind
      dmajs = 0.0_rkind
      nn   = 1.0_rkind
      NV   = 1.0_rkind

    endif


    call fluxavg(nx, nth, theta, b2/nn, jacob, b2navg)
    call fluxavg(nx, nth, theta, nn/b2, jacob, nb2avg)
    call fluxavg(nx, nth, theta, nn/NV, jacob, nNVavg)
    call fluxavg(nx, nth, theta, b2/NV, jacob, b2NVavg)

    ! geometric coefficients in the equation for the flux
    Cgeo_G = nb2avg - 1._rkind/b2navg
    Cgeo_U = b2NVavg/b2navg - nNVavg

    Cgeo_Gcl = nb2avg

  else

    ! Simplified geometry

    B2avg = B0**2*(1.0 + 0.5*epsk**2)


    if (pol_asym) then

      !asymmetries of ES potential and main ion density
      dminphia = Za*(Te/Ti)*AsymPhi(:,1)
      dmajphia = Za*(Te/Ti)*AsymPhi(:,2)
      dNH = AsymN(:,1)
      dNV = AsymN(:,2)

      call asymmetry_an(nx, xn, UU, GG, epsk, invaspct, qmag, nuswca, deltaM, Ai, Aa, Zi, Za, &
                        dNH, dNV, dminphia, dmajphia, dmin, dmaj, dmins, dmajs)

      do i = 1, nx
        nn(i,:) = 1.0 + dmins(i)*cos(theta) + dmajs(i)*sin(theta)
      enddo


    else

      ! Poloidally symmetric case

      dminphia = 0.0_rkind
      dmajphia = 0.0_rkind
      dNH = 0.0_rkind
      dNV = 0.0_rkind

      dmins = 0.0_rkind
      dmajs = 0.0_rkind
      nn   = 1.0_rkind

    endif


    dD2 = 0.5*(dmins**2 + dmajs**2)

    ! geometric coefficients in the equation for the flux

    Cgeo_G = 2.0*epsk*dmins + 2.0*eps2 + dD2
    Cgeo_U = epsk*(dNH - dmins) - dD2 + 0.5*(dmins*dNH + dmajs*dNV)
 
    Cgeo_Gcl = 1.0 + epsk*dmins + 2*eps2  

  endif


  !---------------------------------------------------------------------------
  !---------------------------  Impurity flux  -------------------------------
  !---------------------------------------------------------------------------

  ! General form: Va = -D*grad_ln_na + (K*grad_ln_ni + H*grad_ln_Ti + Vrot)
  !                  = -D*grad_ln_na + Vconv


  ! Pfirsch-Schlüter flux

  !Da_PS   = adps*ma*L11impi*FV**2*q_e*Ti*Cgeo_G*amin**2/(Za**2*q_e**2*B2avg*(dpsidx**2 + 1.e-33))
  Da_PSs=0.0
  Ka_PSs=0.0
  Ha_PSs=0.0
  Va_PSs=0.0
  Vra_PSs=0.0
  do j=1,nions
    Da_PS(:,j)   = adps*qmag**2*rhoLimp2*L11impi*(Cgeo_G/(2.0*eps2))
    !Da_PS   = qmag**2*rhoLimp2*adps*L11impi*FV**2/(R0**2*B2avg)*(Cgeo_G/(2.0*eps2))
    Ka_PS(:,j)   = (Za/Zi(j))*Da_PS(:,j)
    Ha_PS(:,j)   = -((1.0 + (Za/Zi(j))*(C0a(:,j) - 1.0)) + (Cgeo_U/Cgeo_G)*(Za/Zi(j))*(C0a(:,j) + ki(:,j)))*Da_PS(:,j)
    Vra_PS(:,j)  = -Da_PS*grad_ln_na/amin + Ka_PS*grad_ln_ni/amin + Ha_PS*grad_ln_Ti/amin 
    Va_PS(:,j)  = Ka_PS(:,j)*grad_ln_ni(:,j)/amin + Ha_PS(:,j)*grad_ln_Ti/amin
    Da_PSs += Da_PS(:,j)
    Ka_PSs += Ka_PS(:,j)
    Ha_PSs += Ha_PS(:,j)
    Vra_PSs += Vra_PS(:,j)
    Va_PSs += Va_PS(:,j)
  enddo 


  ! Classical flux
  Da_CLs=0.0
  Ka_CLs=0.0
  Ha_CLs=0.0
  Va_CLs=0.0
  Vra_CLs=0.0
  do j=1,nions
    Da_CL(:,j)   = (2.0*eps2*Cgeo_Gcl/Cgeo_G)*Da_PS(:,j)/(2*qmag**2)
    Ka_CL(:,j)   = (Za/Zi(j))*Da_CL(:,j)
    Ha_CL(:,j)   = -(1.0 + (Za/Zi(j))*(C0a(:,j) - 1.0))*Da_CL(:,j)
    Vra_CL(:,j)  = -Da_CL(:,j)*grad_ln_na/amin + Ka_CL(:,j)*grad_ln_ni(:,j)/amin + Ha_CL(:,j)*grad_ln_Ti/amin
    Va_CL(:,j) = Ka_CL(:,j)*grad_ln_ni(:,j)/amin + Ha_CL(:,j)*grad_ln_Ti/amin
    Da_CLs += Da_CL(:,j)
    Ka_CLs += Ka_CL(:,j)
    Ha_CLs += Ha_CL(:,j)
    Vra_CLs += Vra_CL(:,j)
    Va_CLs += Va_CL(:,j)
  enddo

  ! Banana-Plateau flux

  call K_VISC(nx, Ni, Na, Ti, wii, wimpimp, Zi, Za, Ai, Aa, Tauii, &
              Tauimpi, Tauiimp, Tauimpimp, eps2, R0, qmag, ft, y1, y2, y3, y4, &
              K11a, K12a, K22a, K11i, K12i, K22i)


  Da_BPs=0.0
  Ka_BPs=0.0
  Ha_BPs=0.0
  Va_BPs=0.0
  Vra_BPs=0.0
  do j=1,nions
    Da_BP = 1.5*q_e*Ti*(1.0/(1.0/K11a(:,j) + 1.0/K11i(:,j)))/(Za**2*q_e**2*FV**2*Na)
    Ka_BP = (Za/Zi(j))*Da_BP(:,j)
    Ha_BP = ((Za/Zi(j))*(K12i(:,j)/K11i(:,j) - 1.5) - (K12a(:,j)/K11a(:,j) - 1.5))*Da_BP(:,j)

    Vra_BP(:,j) = -Da_BP(:,j)*grad_ln_na/amin + Ka_BP(:,j)*grad_ln_ni(:,j)/amin + Ha_BP(:,j)*grad_ln_Ti/amin
    Va_BP(:,j) = Ka_BP(:,j)*grad_ln_ni(:,j)/amin + Ha_BP(:,j)*grad_ln_Ti/amin
    Da_BPs += Da_BP(:,j)
    Ka_BPs += Ka_BP(:,j)
    Ha_BPs += Ha_BP(:,j)
    Vra_BPs += Vra_BP(:,j)
    Va_BPs += Va_BP(:,j)

  ! Total transport coefficients

  Da = Da_PS + Da_BP + Da_CL ! total diffusion coefficient matrix
  Ka = Ka_PS + Ka_BP + Ka_CL
  Ha = Ha_PS + Ha_BP + Ha_CL

  Das = Da_PSs + Da_BPs + Da_CLs ! total diffusion coefficient
  Kas = Ka_PSs + Ka_BPs + Ka_CLs
  Has = Ha_PSs + Ha_BPs + Ha_CLs

  !Vra = Vra_PS + Vra_BP + Vra_CL

  ! total convective velocity
  Vconvs =0.0_rkind
  Flux_imp = 0.0_rkind
  do j = 1,nions
    Vconv = Ka(:,j)*grad_ln_ni(:,j)/amin + Ha(:,j)*grad_ln_Ti/amin !check 
    Vconvs += Vconv(:,j)

  ! Total surface-averaged flux
    Flux_imp(:,j) = -Da(:,j)*gradNa + Na*Vconv(:,j)
    Flux_imps += Flux_imp(:,j)
  enddo

  return
end subroutine FACIT



!*******************************************************************************
!*********************** Complementary subroutines *****************************
!*******************************************************************************



subroutine fluxavg(nx,nth,thetay,AF,JJ,Aavg) !check

!*******************************************************************************
! Calculates the flux surface average (FSA) of a function AF(x,theta)
!*******************************************************************************
! INPUTS:
! -------
! nx -----> size of radial arrays [-] {int}
! nth ----> size of poloidal arrays [-] {int}
! thetay -> poloidal grid [-] {arr, nth}
! AF -----> function to average [-] {arr, (nx,nth)}
! JJ -----> Jacobian of the coordinate system [-] {arr, (nx,nth)}
!*******************************************************************************
! OUTPUT:
! ------
! Aavg ---> FSA of the AF function [-] {arr, nx}
!*******************************************************************************

  use constants, only: rkind
  implicit none

  real(rkind), dimension(nx,nth) :: AF, JJ
  real(rkind), dimension(nth) :: thetay
  real(rkind), dimension(nx) :: Aavg
  real(rkind) :: denom
  integer   :: nx, nth, ix, ith

  Aavg = 0._rkind

  do ix = 1,nx
    denom = 0._rkind

    do ith=2,nth
      denom = denom + 0.5*(thetay(ith)-thetay(ith-1))*(JJ(ix,ith) + JJ(ix,ith-1))
    enddo

    do ith = 2,nth
      Aavg(ix) = Aavg(ix)+ 0.5*(thetay(ith)-thetay(ith-1))*(AF(ix,ith)*JJ(ix,ith) + AF(ix,ith-1)*JJ(ix,ith-1))
    enddo

    if (denom.eq.0.) then
      Aavg(ix) = sum(AF(ix,:))/nth
    else
      Aavg(ix) = Aavg(ix)/denom
    endif

  enddo

  return
end subroutine fluxavg 




subroutine fluxavgscal(nth,thetay,AF,JJ,Aavg) !check

!*******************************************************************************
! Calculates the flux surface average (FSA) of a function AF(x0, theta) at a
! given radial point x0
!*******************************************************************************
! INPUTS:
! -------
! nth ----> size of poloidal arrays [-] {int}
! thetay -> poloidal grid [-] {arr, nth}
! AF -----> function to average [-] {arr, nth}
! JJ -----> Jacobian of the coordinate system at x0 [-] {arr, nth}
!*******************************************************************************
! OUTPUT:
! ------
! Aavg ---> FSA of the AF function [-] {float}
!*******************************************************************************

  use constants, only: rkind
  implicit none

  real(rkind), dimension(nth) :: AF, JJ
  real(rkind), dimension(nth) :: thetay
  real(rkind) :: denom, Aavg
  integer   :: nth, ith

  Aavg = 0._rkind
  denom = 0._rkind

  do ith=2,nth
    denom = denom + 0.5*(thetay(ith)-thetay(ith-1))*(JJ(ith) + JJ(ith-1))
  enddo

  do ith=2,nth
    Aavg = Aavg+ 0.5*(thetay(ith)-thetay(ith-1))*(AF(ith)*JJ(ith) + AF(ith-1)*JJ(ith-1))
  enddo

  if (abs(denom).gt.0.) then
    Aavg = Aavg/denom
  else
    Aavg = sum(AF)/nth
  endif

  return
end subroutine fluxavgscal



function ftrap(epsK)
!*******************************************************************************
! Trapped particle fraction as a function of the local inverse aspect ratio
!*******************************************************************************
! INPUT:
! ------
! epsK --> local inverse aspect ratio [-] {float}
!*******************************************************************************
! OUTPUT:
! ------
! ftrap -> trapped particle fraction [-] {float}
!*******************************************************************************
  use constants, only: rkind
  implicit none

  real(rkind) :: ftrap
  real(rkind), intent(inout) :: epsK

  ftrap = 1. - (1. - epsK)**2/(sqrt(1. - epsK**2)*(1+1.46*sqrt(epsK)))

  return
end function ftrap



subroutine asymmetry_an(nx, xn, UU, GG, epsk, invaspct, qmag, nuswca, deltaM, Ai, Aa, Zi, Za, &
                        dNH, dNV, dminphia, dmajphia, dmin, dmaj, dmins, dmajs)

!*******************************************************************************
! Poloidal asymmetry of the impurity density distribution, analytical
! calculation in circular geometry
!*******************************************************************************
! INPUTS:
! -------
! nx -------> size of radial arrays [-] {int}
! xn -------> radial grid [-] {arr, nx}
! UU -------> thermodynamic gradient U [-] {arr, nx}
! GG -------> thermodynamic gradient G [-] {arr, nx}
! epsK -----> local inverse aspect ratio [-] {arr, nx}
! invaspct -> inverse aspect ratio [-] {float}
! qmag -----> safety factor [-] {arr, nx}
! nuswca ---> ratio of impurity coll. freq. to cyclotron freq. [-] {arr, nx}
! deltaM ---> rotation strength parameter [-] {arr, nx}
! Ai -------> main ion mass number [-] {float}
! Aa -------> impurity mass number [-] {float}
! Zi -------> main ion charge number [-] {float}
! Za -------> impurity charge number [-] {arr, nx}
! dNH ------> horizontal asymmetry of main ion density [-] {arr, nx}
! dNV ------> vertical asymmetry of main ion density [-] {arr, nx}
! dminphia -> horizontal asymmetry of electrostatic potential [-] {arr, nx}
! dmajphia -> vertical asymmetry of electrostatic potential [-] {arr, nx}
!*******************************************************************************
! OUTPUTS:
! --------
! dmin -----> horizontal asymmetry of impurity density Matrix [-] {arr, nx}
! dmaj -----> vertical asymmetry of impurity density Matrix [-] {arr, nx}
! dmins -----> horizontal asymmetry of impurity density sum [-] {arr, nx}
! dmajs -----> vertical asymmetry of impurity density sum [-] {arr, nx}
!*******************************************************************************

  use constants, only: rkind
  implicit none

  integer :: nx, nions
  real(rkind), dimension(nx) :: epsK, qmag, deltaM, Za
  real(rkind), dimension(nx,nions) :: UU, GG, nuswca, dNH, dNV, Zi
  real(rkind), dimension(nx) :: dminphia, dmajphia, xn
  real(rkind) :: Aa, invaspct
  real(rkind), dimension(nions) :: Ai
  real(rkind), dimension(nx,nions) :: dmin, dmaj
  real(rkind), dimension(nx,nions) :: RR, UG, Ae, AGe, CD0, HH, QQ, FF, KK
   real(rkind), dimension(nx) :: S1, S2, S3
  real(rkind), dimension(nx) :: CD, CDV, RD, DD, num, cosa, sina, dmin, dmaj



  !RR = 0.5*(1-(Ai*Za)/(Aa*Zi))/(epsk + 1.e-33)
  !RR = 0.5*(1-Ai*Za/(Aa*Zi))/(xn+1.e-33)
  RR = 0.0_rkind

  UG = 1 + UU/GG
  !Ae = nuswca*qmag**2/(epsk + 1.e-33)
  do j= 1,nions
    Ae(:,j) = nuswca(:,j)*qmag**2/invaspct
  enddo
  AGe = Ae*GG
  do j = 1,nions
    CD0 = -(epsk + 1.e-33)/UG(:,j)
  enddo
  !HH = 1.0 + deltaM*CD0*RR/GG
  !HH = 1.0_rkind
  !QQ = CD0*(dNV/(epsk + 1.e-33))*UU/GG
  do j=1,nions
    QQ(:,j) = CD0(:,j)*(dNV(:,j)/(epsk + 1.e-33))*(UG(:,j)-1.0)
  
  !FF = CD0*(1-0.5*dNH/(epsk + 1.e-33)*UU/GG - deltaM*(RR/GG)/(epsk + 1.e-33))
  !FF = CD0*(1-0.5*dNH*(UG-1.0)/(epsk + 1.e-33) - deltaM*(RR/GG)/(epsk + 1.e-33))
    FF(:,j) = CD0(:,j)*(1-0.5*dNH(:,j)*(UG(:,j)-1.0)/(epsk + 1.e-33) )
  enddo
  KK = 1.0_rkind

  S1 = 0.0_rkind
  S2 = 0.0_rkind
  S3 = 0.0_rkind

  do j = 1, nions
    S1 += A(:,j)/CD0(:,j)
    S2 += A(:,j)/CD0(:,j) * QQ(:,j) 
    S3 += A(:,j)/CD0(:,j) * FF(:,j) 
  enddo

  CD = FF - 0.5*(dminphia - deltaM + S2)
  CDV = -0.5*(dmajphia + S3)
  RD = sqrt(CD**2 +CDV**2)
  DD(:,j) = RD(:,j)**2 + AGe**2*(RD(:,j)/CD0(:,j))**2

  num(:,j)  = ((AGe(:,j)/CD0(:,j))**2 - 1)*(FF(:,j)/(CD0(:,j)) + 0.5*(dminphia-deltaM)/CD0(:,j)) + &
          (AGe(:,j)/CD0(:,j))*(0.5*dNV(:,j)*(UG(:,j)-1.0)/(epsk + 1.e-33) - 0.5*dmajphia/CD0(:,j))
  cosa(:,j) = RD(:,j)*CD0(:,j)*num(:,j)/DD(:,j)

  num(:,j)  = 2*AGe(:,j)*(FF(:,j)/CD0(:,j) + 0.5*(dminphia-deltaM)/CD0(:,j))+((AGe(:,j)/CD0(:,j))**2-1)*&
          (0.5*dmajphia - 0.5*dNV(:,j)*CD0(:,j)*(UG(:,j)-1.0)/((epsk + 1.e-33)))
  sina(:,j) = RD(:,j)*num(:,j)/DD(:,j)



dmin += dmin
dmaj += dmajs


end subroutine asymmetry_an

!New subroutine for matrix calculation. 

subroutine asymmetry_an_mat(nx, xn,nions, UU, GG, epsk, J, F, nuswca, deltaM, Aa, R0,B0 Za, &
                        dNH, dNV, dminphia, dmajphia, dmin_M, dmaj_M)

!*******************************************************************************
! Poloidal asymmetry of the impurity density distribution, analytical matrix solution
!*******************************************************************************
! INPUTS:
! -------
! nx -------> size of radial arrays [-] {int}
! nions ------> numer of ion species [-] {int}
! xn -------> radial grid [-] {arr, nx}
! UU -------> thermodynamic gradient U [-] {arr, (nx,nions)}
! GG -------> thermodynamic gradient G [-] {arr, (nx,nions)}
! epsK -----> local inverse aspect ratio [-] {arr, nx}
! nuswca ---> ratio of impurity coll. freq. to cyclotron freq. [-] {arr, nx}
! deltaM ---> rotation strength parameter [-] {arr, nx}
! Aa -------> impurity mass number [-] {float}
! Za -------> impurity charge number [-] {arr, nx}
! dNH ------> horizontal asymmetry of main ion density [-] {arr, (nx,nions)}
! dNV ------> vertical asymmetry of main ion density [-] {arr, (nx,nions)}
! dminphia -> horizontal asymmetry of electrostatic potential [-] {arr, nx}
! dmajphia -> vertical asymmetry of electrostatic potential [-] {arr, nx}
!*******************************************************************************
! OUTPUTS:
! --------
! dmin_M -----> horizontal asymmetry of impurity density Matrix solution [-] {arr, nx}
! dmaj_M -----> vertical asymmetry of impurity density Matrix solution [-] {arr, nx}
!*******************************************************************************

  use constants, only: rkind
  implicit none

  integer :: nx, nions
  real(rkind), dimension(nx) :: epsK, deltaM, Za, J,F !check for F and J
  real(rkind), dimension(nx,nions) :: UU, GG, nuswca, dNH, dNV, A_M
  real(rkind), dimension(nx) :: dminphia, dmajphia, xn, detinv, S2, S3, S1
  real(rkind) :: Aa, invaspct, ma, R0, B0
  real(rkind), dimension(nx) :: dmin_M, dmaj_M
  
  ma= Aa*mp

  S1 = 0.0_rkind
  S2 = 0.0_rkind
  S3 = 0.0_rkind
  detinv = 0.0_rkind

  do j = 1,nions
    A_M(:,j) = J*F*ma*nuswca(:,j)/(Za*q_e)
    detinv += 1/(1+A_M(:,j)**2 *(GG(:,j)+UU(:,j))**2) !determinant of the matrix 
    S1 += A_M(:,j)*(GG(:,j)+UU(:,j))
    S2 +=  A_M(:,j)*dNV(:,j)*UU(:,j)
    S3 += 2*A_M(:,j)*epsK*GG(:,j) - dNH(:,j)*A(:,j)*UU(:,j)
  enddo

  !Calculation of the deltas
  dmin_M = detinv*(deltaM-dminphia+S2 - S1*(-dmajphia+S3))
  dmaj_M = detinv*(S1*(deltaM-dminphia+S2) -dmajphia+S3)


end subroutine asymmetry_an_Mat



subroutine asymmetry_fg(nx, nth, theta, BV, RV, jacob, FV, dpsidx, Machi, L11impi, &
                        R0, Ai, Aa, Zi, Za, B2avg, UU, GG, PhiV, NV, Te, Ti, regulopt, &
                        dmin, dmaj, nn)
!*******************************************************************************
! Poloidal asymmetry of the impurity density distribution, iterative
! calculation in full geometry
!*******************************************************************************
! INPUTS:
! -------
! nx -------> size of radial arrays [-] {int}
! nth ------> size of poloidal arrays [-] {int}
! theta ----> poloidal grid [-] {arr, nth}
! BV -------> magnetic field [T] {arr, (nx,nth)}
! RV -------> major radius [m] {arr, (nx,nth)}
! jacob ----> Jacobian of the coordinate system [m/T] {arr, (nx,nth)}
! FV -------> poloidal current flux function [T*m] {arr, nx}
! dpsidx ---> radial derivative of poloidal flux [V*s/-] {arr, nx}
! Machi ----> Mach number of main ion [-] {arr, nx}
! L11impi --> impurity - main ion collision frequency [1/s] {arr, nx}
! R0 -------> major radius at magnetic axis [m] {float}
! Ai -------> main ion mass number [-] {float}
! Aa -------> impurity mass number [-] {float}
! Zi -------> main ion charge number [-] {float}
! Za -------> impurity charge number [-] {arr, nx}
! B2avg ----> FSA of magnetic field squared <BV**2> [T^2] {arr, nx}
! UU -------> thermodynamic gradient U [-] {arr, nx}
! GG -------> thermodynamic gradient G [-] {arr, nx}
! PhiV -----> poloidal asymmetry of electrostatic potential [-] {arr, (nx,nth)}
! NV -------> poloidal asymmetry of main ion density [-] {arr, (nx,nth)}
! Te -------> electron temperature [eV] {arr, nx}
! Ti -------> main ion temperature [eV] {arr, nx}
! regulopt -> options for iterative calculations [-] {arr, 4}
!*******************************************************************************
! OUTPUTS:
! --------
! dmin -----> horizontal asymmetry of impurity density [-] {arr, nx}
! dmaj -----> vertical asymmetry of impurity density [-] {arr, nx}
! nn -------> poloidal asymmetry of the impurity density Na/<Na> [-] {arr, (nx,nth)}
!*******************************************************************************

  use constants, only: rkind, mp, q_e, pi
  implicit none

  integer :: nx, nth, nions
  real(rkind), dimension(nth) :: theta
  real(rkind), dimension(nx) :: FV, dpsidx, Machi, Za, Te, Ti, B2avg
  real(rkind), dimension(nx,nions) :: UU, GG, L11impi
  real(rkind) :: R0, Aa
  real(rkind), dimension(nions) :: Zi, Ai
  real(rkind), dimension(nx, nth) :: BV, RV, jacob, PhiV, NV
  real(rkind), dimension(4) :: regulopt

  real(rkind), dimension(nx,nions) :: dmin, dmaj
  real(rkind), dimension(nx, nth) :: nn


  real(rkind), dimension(nx,nions) :: asym_error, Factrot0
  real(rkind), dimension(1:nth+1,1:nth+1) :: AAA, LL, CC, DD
  real(rkind), dimension(1:nth+1) :: BB, nnya
  real(rkind), dimension(1:nth,nions) :: Apsi, FFF, GGG, HHH
  real(rkind), dimension(1:nth) :: dum, nnx, nny, b2, nnp, TermGG, TermUU
  real(rkind), dimension(size(AAA,1)) :: work  ! work array for LAPACK
  real(rkind)  :: b2snavg, nsb2avg, b2sNNavg, FactV
  !real(rkind)  :: TermUUneo, TermRRneo, TermRRneoa, TermRRneob, TermGGcl
  !real(rkind)  :: TermGGneo, TermUUneoa, TermUUneob, TermRRcla, TermRRclb, TermRRcl
  real(rkind)  :: err, prog, regulweight, dtheta, Erreur, dumscal
  integer, dimension(size(AAA,1)) :: ipiv   ! pivot indices
  integer :: n, info, ix, it, ierr, ierrmax



  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  !allocate(AA(nth+1,nth+1),LL(nth+1,nth+1),CC(nth+1,nth+1),DD(nth+1,nth+1),BB(nth))

  err         = regulopt(1)
  prog        = regulopt(2)
  regulweight = regulopt(3)
  ierrmax     = int(regulopt(4))


  ! Poloidal asymmetry
  
  nn = 1.0_rkind

  do ix = 1,nx
    do j= 1,nions

 	    Factrot0(ix,j) = (Aa/Ai(j))*(Machi(ix)/R0)**2
    enddo

    do it=1,nth
        b2(it) = BV(ix,it)**2/B2avg(ix)
        do j = 1,nions
          Apsi(it,j) = jacob(ix,it)*FV(ix)*(Aa*mp)*L11impi(ix,nions)/(q_e*Za(ix))/((dpsidx(ix))**2+1.e-33)
        enddo
    enddo

   	call fluxavgscal(nth, theta , b2/NV(ix,:), jacob(ix,:), b2sNNavg)

  	nnp = 2
  	nnx = nn(ix,:)
  	ierr = 1
  	Erreur = 2*err

  	do while ( Erreur>err .and. ierr<ierrmax)

      AAA = 0._rkind
  		LL  = 0._rkind
  		BB  = 0._rkind

  		call fluxavgscal(nth, theta, b2/nnx, jacob(ix,:), b2snavg)

  		TermGG = 1. - b2/nnx/b2snavg
  		TermUU = b2/NV(ix,:) - b2sNNavg*b2/nnx/b2snavg
      do j= 1, nions
        FFF  = Apsi(:,j)*( GG(ix,j) + b2/NV(ix,:)*UU(ix,j) )
        GGG = -Za(ix)*Te(ix)/Ti(ix)*(PhiV(ix,:)-PhiV(ix,1)) + Factrot0(ix,j)*(RV(ix,:)**2-RV(ix,1)**2)
        HHH = Apsi(:,j)*b2/b2snavg*(GG(ix,j) + b2sNNavg*UU(ix,j) )
      enddo

      ! AA*n = BB  CHECK
      ! Remplissage matrices AA, BB & LL

  		do it = 2,nth-1

  			dtheta = 0.5*(theta(it+1) - theta(it-1))
  			AAA(it,it-1) = -0.5/dtheta
  			AAA(it,it) = -FFF(it)-0.5/dtheta*(GGG(it+1)-GGG(it-1))
  			AAA(it,it+1) = 0.5/dtheta

  			LL(it, it-1) = 1.
  			LL(it, it) = -2.
  			LL(it, it+1) = 1.

  			BB(it) = -HHH(it)

  		enddo

  		it = 1
  		dtheta = 0.5*(theta(2)-theta(nth) + 2*PI)
  		AAA(it,it)   = -FFF(it) - (0.5/dtheta)*(GGG(it+1) - GGG(nth-1))
  		AAA(it,it+1) = 0.5/dtheta
  		AAA(it,nth)  = -0.5/dtheta
  		LL(it, it)   = -2.
  		LL(it, it+1) = 1.
  		LL(it, nth)  = 1.
  		BB(it) = -HHH(it)


  		it = nth
  		dtheta = 0.5*(theta(1)-theta(nth-1) + 2*PI)
  		AAA(it,it-1) = -0.5/dtheta
  		AAA(it,it)   = -FFF(it)-0.5/dtheta*(GGG(1)-GGG(it-1))
  		AAA(it,it+1) = 0.5/dtheta
  		LL(it, it-1) = 1.
  		LL(it, it)   = -2.
  		LL(it, it+1) = 1.
  		BB(it) = -HHH(it)

  		it = nth+1
  		dtheta = 0.5*(theta(2)-theta(nth) + 2*PI)
      		AAA(it,it) = -1.
  		AAA(it,1)  = 1
  		BB(it) = 0.

  		! CC = transpose(AA)*AA+regulweight*transpose(LL)*LL
  		CC = transpose(AAA)
  		CC = matmul(CC,AAA)
  		DD = transpose(LL)
  		DD = matmul(DD,LL)
  		CC = CC + regulweight*DD

  		n  = size(CC,1)
  		DD = CC
  		call DGETRF(n, n, DD, n, ipiv, info)

   		if (info /= 0) then
			   stop 'Matrix is numerically singular!'
  		end if

  		! DGETRI computes the inverse of a matrix using the LU factorization
  		! computed by DGETRF.
   		call DGETRI(n, DD, n, ipiv, work, n, info)

      if (info /= 0) then
			   stop 'Matrix inversion failed!'
  		end if

      CC = matmul(DD,transpose(AAA))
  		nnya = matmul(CC,BB)
      !nnya = 1.

  		call fluxavgscal(nth, theta, nnya(1:nth), jacob(ix,:), dumscal)

  		nny = nnya(1:nth)/dumscal
  		nny = max(nny, 1.e-5)

  		nnp = nnx
  		nnx = prog*nnp + (1-prog)*nny

  		do it = 2, nth-1
			   dum(it) = (log(nnx(it+1))-log(nnx(it-1))-GGG(it+1)+GGG(it-1))/(theta(it+1)-theta(it-1))
  		enddo

  		it = 1
  		dum(it) = (log(nnx(it+1))-log(nnx(it))-GGG(it+1)+GGG(it))/(theta(it+1)-theta(it))

  		it = nth
  		dum(it) = (log(nnx(it))-log(nnx(it-1))-GGG(it)+GGG(it-1))/(theta(it)-theta(it-1))

  		Erreur = maxval(abs(dum-FFF+HHH/nnx))

  		ierr = ierr+1

    end do

  	asym_error(ix) = Erreur

  	nn(ix,:) = nnx
  	dmin(ix) = 2.*sum((nnx-1.)*cos(theta))/float(nth)
  	dmaj(ix) = 2.*sum((nnx-1.)*sin(theta))/float(nth)

  enddo

end subroutine asymmetry_fg







function ki_Redl(nuistar, ft, Zeff)

!*******************************************************************************
! Neoclassical ion flow coefficient, fitted w.r.t. NEO in Redl PoP (2021)
!*******************************************************************************
! INPUTS:
! -------
! nuistar -> main ion collisionality [-]
! ft ------> trapped particle fraction [-]
! Zeff ----> effective charge [-]
!*******************************************************************************
! OUTPUT:
! ------
! ki ------> main ion flow coefficient [-]
!*******************************************************************************

  use constants, only: rkind
  implicit none

  real(rkind), intent(in) :: nuistar, ft, Zeff
  real(rkind) :: alpha0
  real(rkind) :: ki_Redl

  alpha0 = -(0.62 + 0.055*(Zeff - 1.0))*(1.0 - ft)/((0.53 + 0.17*(Zeff - 1.0))*&
            (1.0 - (0.31 - 0.065*(Zeff - 1.0))*ft - 0.25*ft**2))

  ki_Redl = ((alpha0 + 0.7*Zeff*SQRT(ft*nuistar))/(1.0 + 0.18*SQRT(nuistar))-&
              0.002*nuistar**2*ft**6)/(1.0 + 0.004*nuistar**2*ft**6)

  return
end function ki_Redl




function C2(alpha, g, f1, f2, Aimp, Ai)

!*******************************************************************************
! new parametrized version of the C2 function
! see Hirshman-Sigmar NF (1981), page 61 (1138)
!*******************************************************************************
! INPUTS:
! -------
! alpha -> impurity strength parameter [-]
! g -----> main ion collisionality parameter [-]
! f1 ----> f1(Zimp) factor [-]
! f2 ----> f2(Zimp) factor [-]
!*******************************************************************************
! OUTPUT:
! ------
! C2 ----> term in impurity-ion friction C0a [-]
!*******************************************************************************

  use constants, only: rkind
  implicit none

  real(rkind), intent(in) :: alpha, g, f1, f2, Aimp, Ai
  real(rkind) :: C2

  C2 = 1.5/(1.0 + (Ai/Aimp)*f1) - (0.29 + 0.68*alpha)/(0.59 + alpha + (1.34 + f2)/g**2)

  return
end function C2



subroutine facs(nx, Zimp, ft, f1, f2, f3, y1, y2, y3, y4, adps)

!*******************************************************************************
! set of factors fitted with respect to NEO
!*******************************************************************************
! INPUTS:
! -------
! nx ---> size of radial arrays [-]
! Zimp -> impurity charge [-]
! ft ---> trapped particle fraction [-]
!*******************************************************************************
! OUTPUTS:
! --------
! f1 --> low collisionality saturation of PS C0a coefficient [-]
! f2 --> correction to low-Z impurity-ion friction in C0a coefficient [-]
! f3 --> factor in i-e heat exchange contribution to C0a [-]
! y1 --> factor of the (1,1) banana viscosity coefficient of the impurity [-]
! y2 --> factor of the (1,2) plateau viscosity coefficient of the impurity [-]
! y3 --> factor of the (1,2) Pfirsch-Schlüter viscosity coefficient of the impurity [-]
! y4 --> factor of the (1,2) banana viscosity coefficient of the main ion [-]
!*******************************************************************************

  use constants, only: rkind
  implicit none


  integer :: nx
  !real(dp) :: Zimp
  real(rkind), dimension(nx) :: Zimp
  real(rkind), dimension(nx) :: ft

  !real(dp) :: f1, f2, f3
  real(rkind), dimension(nx) :: f1, f2, f3
  real(rkind), dimension(nx) :: adps, y2, y3
  real(rkind), dimension(nx) :: w11, w12, w13, w14, y1
  real(rkind), dimension(nx) :: wp1, wp2, wp3, wp4, yp, y4


  ! f's:
  f1 =  (-6.83808564e5 + 2.46855534e6*Zimp)/( 1 + 6.04692708e5*Zimp**1.61470425)
  f2 = (88.28389935 + 10.50852772*Zimp)/( 1 + 0.2157175*Zimp**2.57338463)
  f3 = (-3.15171054e6 + 1.92908543e6*Zimp)/(1+5.26716920e6*Zimp**8.33610108e-01)


  ! y's:
  w11 = 1.23805214e-05*Zimp**3 - 1.03611576e-03*Zimp**2 + 1.85221287e-02*Zimp + 1.29758029
  w12 = -5.84996779e-05*Zimp**3 + 4.79623045e-03*Zimp**2 -1.01924030e-01*Zimp - 5.81816797
  w13 = 5.98152997e-05*Zimp**3 - 4.84255587e-03*Zimp**2 + 1.03585208e-01*Zimp + 5.71139227
  w14 = -1.33253954e-05*Zimp**3 + 1.04416741e-03*Zimp**2 - 1.94973421e-02*Zimp - 1.10061492

  y1 = w11*ft**2 + w12*ft + w13*ft**0.5 + w14

  wp1 = ( 0.11603574 + 0.47297835*Zimp**0.94456671)/( 1 + 0.1245426*Zimp**1.20221441)
  wp2 = ( 0.6327602 - 2.92116611*Zimp**1.06310182)/( 1 + 0.30469549*Zimp**1.16495283)
  wp3 = ( -0.69318477 + 2.85619511*Zimp**1.11765611)/( 1 + 0.34127361*Zimp**1.1952383)
  wp4 = ( 1.80217558 - 1.0080554*Zimp**0.96345934)/( 1 + 0.51677339*Zimp**1.11131896)

  yp = wp1*ft**2 + wp2*ft + wp3*ft**0.5 + wp4
  y4 = yp/y1

  y2 = 21.31*ft**3 - 21.88*ft**2 + 7.316*ft + 0.6264
!  y3 = ((4.2857e5 - 4.4978e5*Zimp)/(1 - 1.3557e5*Zimp**1.9))*((-9.09204093e6 + 8.15802759e6*Zimp)/(1.0 + 3.25263641e6*Zimp**1.07640221))
  y3 = 8.85/Zimp**2.98 - 7.96/Zimp**1.82 - 9.27/Zimp**1.98 + 8.34/Zimp**0.82

!  adps = (-110378.3491 + 753838.926571*Zimp**1.06107833841)/(1+ 1007273.22737*Zimp)
  adps = -0.109/Zimp + 0.743*Zimp**0.06


end subroutine facs


subroutine K_VISC(nx, ni, nimp, Ti, wii, wimpimp, Zi, Zimp, Ai, Aimp, Tauii, &
                  Tauimpi, Tauiimp, Tauimpimp, eps2, R0, qmag, ft, y1, y2, y3, y4, &
                  K11a, K12a, K22a, K11i, K12i, K22i)

!*******************************************************************************
! Calculates the positive definite matrix of neoclassical viscosity coefficients
! in a fully analytical way by solving for the coefficients in the individual
! banana, plateau and Pfirsch-Schlüter collisionality regimes and using a
! rational approximation to interpolate between them, as well as evaluating
! kinetic integrals for a Maxwellian distribution
!*******************************************************************************
! INPUTS:
! -------
! nx --------> size of radial arrays [-]
! ni --------> main ion density [1/m^3]
! nimp ------> impurity density [1/m^3]
! Ti --------> main ion temperature [eV]
! wii -------> ion transit frequency [1/s]
! wimpimp ---> impurity transit frequecy [1/s]
! Zi --------> main ion charge [-]
! Zimp ------> impurity charge [-]
! Ai --------> main ion mass [-]
! Aimp ------> impurity mass [-]
! Tauii -----> ion-ion collision time [s]
! Tauimpi ---> impurity-ion collision time [s]
! Tauiimp ---> ion-impurity collision time [s]
! Tauimpimp -> impurity impurity collision time [s]
! eps2 ------> local inverse aspect ratio squared [-]
! R0 --------> major radius at magnetic axis [m]
! qmag ------> safety factor [-]
! ft --------> trapped particle fraction [-]
! y1 --------> fitted factor from facs subroutine [-]
! y2 --------> fitted factor from facs subroutine [-]
! y3 --------> fitted factor from facs subroutine [-]
! y4 --------> fitted factor from facs subroutine [-]
!*******************************************************************************
! OUTPUTS:
! --------
! K11a ------> (1,1) impurity viscosity coefficient [kg/(m*s)]
! K12a ------> (1,2)=(2,1) impurity viscosity coefficient [kg/(m*s)]
! K22a ------> (2,2) impurity viscosity coefficient [kg/(m*s)]
! K11i ------> (1,1) main ion viscosity coefficient [kg/(m*s)]
! K12i ------> (1,2)=(2,1) main ion viscosity coefficient [kg/(m*s)]
! K22i ------> (2,2) main ion viscosity coefficient [kg/(m*s)]
!*******************************************************************************

  use constants, only: rkind, mp, q_e, pi, sqrt2
  implicit none

  integer :: nx,nions
  real(rkind), dimension(nx,nions) :: ni, wii, Tauii, Tauimpi, Tauiimp
  real(rkind), dimension(nx) :: nimp, Ti, wimpimp, Tauimpimp, eps2, qmag, ft
  real(rkind), dimension(nions) :: Zi, Ai, mi
  real(rkind) ::  Aimp, R0, mimp
  !real(dp) :: Zimp
  real(rkind), dimension(nx) :: Zimp
  real(rkind), dimension(nx,nions) :: K11a, K12a, K22a, K11i, K12i, K22i
  !real(dp) :: f1, f2, f3
  real(rkind), dimension(nx) :: y1, y2, y3, y4

  real(rkind), dimension(nx) :: fac_a_P, K11aP, K12aP, K22aP
  real(rkind), dimension(nx,nions) :: fac_i_P, K11iP, K12iP, K22iP
  real(rkind) :: r00, r01, r11
  real(rkind), dimension(nions) :: xai, xia, x2ai, x2ia, xfac_ai, xfac_ia
  real(rkind) :: qaa00, qaa01, qaa11, qai00, qia00, qai01, qia01, qai11, qia11
  real(rkind), dimension(nx,nions) :: fac_qai_PS, fac_qia_PS
  real(rkind), dimension(nx,nions) :: qa00, qi00, qa01, qi01, qa11, qi11, Qa, Qi
  real(rkind), dimension(nx,nions) :: la11, la12, la22, li11, li12, li22
  real(rkind), dimension(nx) :: fac_imp_PS 
  real(rkind), dimension(nx,nions) :: fac_ion_PS
  real(rkind), dimension(nx,nions) :: K11aPS, K12aPS, K22aPS, K11iPS, K12iPS, K22iPS
  real(rkind), dimension(nx) :: fac_B
  real(rkind), dimension(nx,nions) :: nuDai, nuD2ai, nuD4ai
  real(rkind), dimension(nx,nions) :: nuDia, nuD2ia, nuD4ia
  real(rkind), dimension(nx) :: nuDaa, nuD2aa, nuD4aa
  real(rkind), dimension(nx,nions) :: nuDii, nuD2ii, nuD4ii
  real(rkind), dimension(nx,nions) :: K11aB, K12aB, K22aB, K11iB, K12iB, K22iB


  mimp = Aimp*mp
  mi = Ai*mp


  ! Plateau regime

  fac_a_P = nimp*(q_e*Ti)*SQRT(pi)/(3.0*wimpimp)
  do j = i,nions
  fac_i_P(:,j) = ni(:,j)*(q_e*Ti)*SQRT(pi)/(3.0*wii(:,j))
  enddo

  K11aP = fac_a_P*2.0
  K12aP = fac_a_P*2.0*3.0
  K22aP = fac_a_P*2.0*3.0*4.0

  K11iP = fac_i_P*2.0
  K12iP = fac_i_P*2.0*3.0
  K22iP = fac_i_P*2.0*3.0*4.0

  ! Pfirsch-Schlüter regime

  r00 = 1.0/sqrt2
  r01 = 1.5/sqrt2
  r11 = 3.75/sqrt2

  xai = SQRT(Aimp/Ai)
  xia = 1.0/xai

  x2ai = xai**2
  x2ia = xia**2

  xfac_ai = (1+x2ai)**0.5
  xfac_ia = (1+x2ia)**0.5

  qaa00 = 8.0/2**1.5
  qaa01 = 15.0/2**2.5
  qaa11 = 132.5/2**3.5

  qai00 = (3.0+5.0*x2ai)/xfac_ai**3
  qia00 = (3.0+5.0*x2ia)/xfac_ia**3

  qai01 = 1.5*(3.0+7.0*x2ai)/xfac_ai**5
  qia01 = 1.5*(3.0+7.0*x2ia)/xfac_ia**5

  qai11 = (35.0*x2ai**3 + 38.5*x2ai**2 + 46.25*x2ai + 12.75)/xfac_ai**7
  qia11 = (35.0*x2ia**3 + 38.5*x2ia**2 + 46.25*x2ia + 12.75)/xfac_ia**7

  do j= 1,nions
    fac_qai_PS(:,j) = (ni(:,j)*Zi(j)**2/(nimp*Zimp**2))
    fac_qia_PS(:,j) = (nimp*Zimp**2/(ni(:,j)*Zi(j)**2))
  enddo

  do j=1,nions
    qa00(:,j) = fac_qai_PS(:,j)*qai00 +  qaa00 - r00
    qi00(:,j) = fac_qia_PS(:,j)*qia00 + qaa00 - r00

    qa01(:,j) = fac_qai_PS(:,j)*qai01 + qaa01 - r01
    qi01(:,j) = fac_qia_PS(:,j)*qia01 + qaa01 - r01

    qa11(:,j) = fac_qai_PS(:,j)*qai11 + qaa11 - r11
    qi11(:,j) = fac_qia_PS(:,j)*qia11 + qaa11 - r11
  enddo

  Qa = 0.4*(qa00*qa11-qa01*qa01)
  Qi = 0.4*(qi00*qi11-qi01*qi01)


  la11 = qa11/Qa
  la12 = 3.5*(qa11+qa01)/Qa
  la22 = 12.25*(qa11+qa00+2*qa01)/Qa

  li11 = qi11/Qi
  li12 = 3.5*(qi11 + qi01)/Qi
  li22 = 12.25*(qi11+qi00 + 2.0*qi01)/Qi

  fac_imp_PS = nimp*(q_e*Ti)*Tauimpimp
  do j= 1,nions
    fac_ion_PS(:,j) = ni(:,j)*(q_e*Ti)*Tauii(:,j)
  enddo

  do j=1,nions
    K11aPS(:,j) = fac_imp_PS*la11(:,j)
    K12aPS(:,j) = fac_imp_PS*la12(:,j)
    K22aPS(:,j) = fac_imp_PS*la22(:,j)
  enddo

  K11iPS = fac_ion_PS*li11
  K12iPS = fac_ion_PS*li12
  K22iPS = fac_ion_PS*li22



  ! Banana regime

  fac_B = (ft/(1.0 - ft))*(2.0*R0**2*qmag**2/(3.0*eps2))

  ! Maxwellian integrals

  do j = 1,nions
    nuDai(:,j)  = (xfac_ai + x2ai*LOG(xai/(1 + xfac_ai)))/Tauimpi(:,j)
    nuD2ai(:,j) = 1.0/(xfac_ai*Tauimpi(:,j))
    nuD4ai(:,j) = 2.0*(1.0 + 1.25*x2ai)/(xfac_ai**3*Tauimpi(:,j))

    nuDia(:,j)  = (xfac_ia + x2ia*LOG(xia/(1 + xfac_ia)))/Tauiimp(:,j)
    nuD2ia(:,j) = 1.0/(xfac_ia*Tauiimp(:,j))
    nuD4ia(:,j) = 2.0*(1.0 + 1.25*x2ia)/(xfac_ia**3*Tauiimp(:,j))
  enddo

  nuDaa  = (sqrt2 + LOG(1/(1 + sqrt2)))/Tauimpimp
  nuD2aa = 1.0/(sqrt2*Tauimpimp)
  nuD4aa = 4.5/(2.0**1.5*Tauimpimp)

  nuDii  = (sqrt2 + LOG(1/(1 + sqrt2)))/Tauii
  nuD2ii = 1.0/(sqrt2*Tauii)
  nuD4ii = 4.5/(2.0**1.5*Tauii)

  ! Banana regime viscosity coefficient

  do j=1,nions
    K11aB(:,j) = fac_B*nimp*mimp*(nuDai(:,j) + nuDaa)
    K12aB(:,j) = fac_B*nimp*mimp*(nuD2ai(:,j) + nuD2aa)
    K22aB(:,j) = fac_B*nimp*mimp*(nuD4ai(:,j) + nuD4aa)

    K11iB(:,j) = fac_B*ni(:,j)*mi(j)*(nuDia(:,j) + nuDii(:,j))
    K12iB(:,j) = fac_B*ni(:,j)*mi(j)*(nuD2ia(:,j) + nuD2ii(:,j))
    K22iB(:,j) = fac_B*ni(:,j)*mi(j)*(nuD4ia(:,j) + nuD4ii(:,j))
  enddo

  ! total viscosity coefficients, rational approximation for interpolation
  do j= 1,nions 
  K11a(:,j) = y1*K11aB(:,j)/((1.0 + y1*K11aB(:,j)/(K11aP))*(1.0 + K11aP/(K11aPS(:,j))))
  K12a(:,j) = K12aB(:,j)/((1.0 + K12aB(:,j)/(y2*K12aP))*(1.0 + y2*K12aP/(y3*K12aPS(:,j))))
  K22a(:,j) = K22aB(:,j)/((1.0 + K22aB(:,j)/(K22aP))*(1.0 + K22aP/(K22aPS(:,j))))

  K11i(:,j) = K11iB(:,j)/((1.0 + K11iB(:,j)/(K11iP(:,j)))*(1.0 + K11iP(:,j)/(K11iPS(:,j))))
  K12i(:,j) = y4*K12iB(:,j)/((1.0 + y4*K12iB/(K12iP(:,j)))*(1.0 + K12iP(:,j)/(K12iPS(:,j))))
  K22i(:,j) = K22iB(:,j)/((1.0 + K22iB(:,j)/(K22iP(:,j)))*(1.0 + K22iP(:,j)/(K22iPS(:,j))))



  return
end subroutine K_VISC




