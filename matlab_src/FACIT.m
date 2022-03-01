% Collisional impurity flux
% --------------------------------------------------------------------------------
% out = FACIT(in)
%Te,Ti,Ni,Nimp,qmag,B0,R0,x,dpsidx,jacob,BV,RV,FV,invaspct,Ai,Zi,AZV,PhiV,NV,Machi,param)
% 
% Inputs: in.
% Te, Ti 	: electron and ion temperature profile (eV)
% Ni, Nimp	: ion and impurity density profiles
% gradTi, gradNi, gradNimp : radial derivative of Ti, Ni and Nimp profiles with respect to x = r/a
% Alternatively :
% grad_ln_Ti, grad_ln_Ni, grad_ln_Nimp : radial derivative of ln(Ti), ln(Ni) and ln(Nimp) profiles with respect to x = r/a
% qmag 		: safety factor
% B0, R0	: magnetic field and major radius (scalars)
% x		: r/a
% invaspct	: a/R0
% Ai, Zi	: mass number and charge of the main ion
% Aimp, Zimp(x)	: mass number and charge of the impurity species
% PhiV		: Phi(theta) : electrostatic potential. 
% BV	 	: B(x,theta) : magnetic induction (T)
% RV	 	: R(x,theta) : major radius (m)
% FV	 	: F(x) : R*B_phi 
% NV	 	: N(x,theta) : main ion density. 
% Machi 	: toroidal ion Mach number = mi*V^2/(2*Ti)
% param		: [tolerance, mix, regulweight, iter max]
% theta    	: poloidal angle used for the metrics
% dpsidx 	: derivarive of poloidal flux with respect to coordinate x;
% jacob 	: jacobian of the coordinate system (x,theta,phi)
% ionregime     : O = banana regime, 1 = extended model by Daniel Farajo & Clemente Angioni
% geom		: O = analytical model, 1 = realistic flux surface geometry
% pol_asym      : 0 = no poloidal asymmetries, 1 = with poloidal asymmetries
% 
% out. : 
% DPSa 	: Pirsh-Schluter diffusion coefficient
% Vra 	: radial velocity
% Vra0 	: radial velocity without asymmetry
% dmin 	: horizontal asymmetry
% dmaj 	: vertical asymmetry
% Da 	: diffusion coefficient 
% Va	: pinch velocity
% Da0 	: diffusion coefficient without asymmetry
% Va0	: pinch velocity without asymmetry
% ki	: neoclassical coefficient for poloidal ion velocity, equal to -1.017 in banana regime
% nn 	: poloidal distribution of the impurity
% --------------------------------------------------------------------------------

function out = FACIT(in);

% Inputs
Ai = in.Ai;
Zi = in.Zi;
Aimp = in.Aimp;
Zimp = in.Zimp;
Te = in.Te;
Ti = in.Ti;
Ne = in.Ne;
Ni = in.Ni;
Nimp = in.Nimp;
x = in.x;

if isfield(in,'gradTi')==1
	grad_ln_Ti = in.gradTi./Ti;
else
	grad_ln_Ti = in.grad_ln_Ti;
end
if isfield(in,'gradNi')==1
	grad_ln_ni = in.gradNi./Ni;
else
	grad_ln_ni = in.grad_ln_Ni;
end
if isfield(in,'gradNimp')==1
	grad_ln_na = in.gradNimp./Nimp;
else
	grad_ln_na = in.grad_ln_Nimp;
end
if isfield(in,'pol_asym')==1
	pol_asym = in.pol_asym;
else
	pol_asym = 1;
end
qmag = in.qmag;
Machi = in.Machi;
invaspct = in.invaspct;
B0 = in.B0;
R0 = in.R0;
ionregime = in.ionregime;
geom = in.geom;
if geom > 0
	dpsidx = in.dpsidx;
	jacob = in.jacob;
	BV = in.BV;
	RV = in.RV;
	FV = in.FV;
	NV = in.NV;
	PhiV = in.PhiV;
	theta = in.theta';
	param = in.param;
	if length(param)==0
		err = 1e-2; % Convergence on parallel momentum equation
		prog = 0.5;
		regulweight = 1.e-5;
		ierrmax = 1e2;
	else
		err = param(1); % Convergence on parallel momentum equation
		prog = param(2);
		regulweight = param(3);
		ierrmax = param(4);
	end
else
	AsymN = in.AsymN;
	AsymPhi = in.AsymPhi;
end
if length(Zimp)==1
	Zimp = Zimp*ones(size(x));
end
amin  = invaspct*R0;
if isfield(in,'dpsidx')==1
	if isfield(in,'FV')==1
		FV = in.FV;
	else
		FV = R0*B0*ones(size(in.x));
	end
	qmag_metrics = amin^2*FV.*in.x./(R0.*in.dpsidx);
else
	qmag_metrics = qmag;
end

interpmeth = 'spline';

eps0 = 8.8542e-12;
mu0 = 4*pi*1e-7;
me = 9.1096e-31;
mp = 1.6726e-27;
qe = 1.602e-19;

mi = Ai*mp;

if geom == 0
	epsilon = x*invaspct;
	B2avg = B0^2.*(1+0.5*epsilon.^2);
else
	for ix=1:length(x)
		epsilon(ix,1) = (max(BV(ix,:))-min(BV(ix,:)))./(max(BV(ix,:))+min(BV(ix,:)));
	end
	B2avg = fluxavg(BV.^2,jacob);
end

FTRAP = sqrt(2*epsilon);

mimp = Aimp*mp;

if geom == 0
	if length(AsymPhi(:,1))>1
	  dphia = Zimp.*Te./Ti.*AsymPhi(:,1);
	  Dphia = Zimp.*Te./Ti.*AsymPhi(:,2);
	else
	  dphia = Zimp.*Te./Ti.*AsymPhi(1);
	  Dphia = Zimp.*Te./Ti.*AsymPhi(2);
	end
	if length(AsymN(:,1))>1
		dNH =  AsymN(:,1);
		dNV =  AsymN(:,2);
	else
		dNH =  AsymN(1);
		dNV =  AsymN(2);
	end
	deltaM = 2*mimp/mi*Machi.^2.*epsilon;
end

%************************
% FRICTION NEOCLASSIQUE *
%************************
% COLLISIONS
%***********
% Logarithme coulombien (NRL)
LneeNRL=23.5-0.5*log(Ne./1e6)+1.25*log(Te)-sqrt(1e-5+(1/16)*(log(Te)-2).^2);
for p=1:length(Ne)
   if (Ti(p)*me./(mi)<Te(p) & Te(p)<10*Zi^2)
      LneiNRL(p,1)=23-0.5*log(Zi^2*Ne(p)/1e6)+1.5*log(Te(p));
   elseif (Ti(p)*me./(mi)<10*Zi^2 & Te(p)>10*Zi^2)
      LneiNRL(p,1)=24-0.5*log(Ne(p)./1e6)+log(Te(p));
   else
      LneiNRL(p,1)=30-log(Zi^2/Ai*(Ni(p)/1e6)^0.5)+1.5*log(Ti(p));
   end
   if (Ti(p)*me./(Aimp*mp)<Te(p) & Te(p)<10*Zimp(p)^2)
      LneimpNRL(p,1)=23-log(Zimp(p)*(Ne(p)/1e6)^0.5)+1.5*log(Te(p));
   elseif (Ti(p)*me./(Aimp*mp)<10*Zimp(p)^2 & Te(p)>10*Zimp(p)^2)
      LneimpNRL(p,1)=24-0.5*log(Ne(p)./1e6)+log(Te(p));
   else
      LneimpNRL(p,1)=30-log(Zimp(p)^2/Aimp*(Nimp(p)/1e6)^0.5)+1.5*log(Ti(p));
   end
end
LniimpNRL = 23-log(Zi*Zimp.*sqrt((Ni./1e6).*Zi^2+(Nimp./1e6).*Zimp.^2))+1.5*log(Ti);
LnimpimpNRL = 23-log(Zimp.*Zimp.*sqrt((Nimp./1e6).*Zimp.^2+(Nimp./1e6).*Zimp.^2))+1.5*log(Ti);
LniiNRL=23-log(Zi*Zi*sqrt(2.*(Ni./1e6)*Zi^2))+1.5*log(Ti);

% Temps intercollision (Braginskii)
Tauee=(3*eps0^2*me^2*(2*pi*qe/me*Te).^1.5)./(qe^4*Ne.*LneeNRL);	
Tauei=(3*eps0^2*me^2*(2*pi*qe/me*Te).^1.5)./(qe^4*Zi^2*Ni.*LneiNRL);	
Tauie=(3*eps0^2*mi^2*(2*pi*qe/mi*Ti).^1.5)./(qe^4*Zi^2*Ne.*LneiNRL);	
Tauii=(3*eps0^2*mi^2*(2*pi*qe/mi*Ti).^1.5)./(qe^4*Zi^4*Ni.*LniiNRL);	
Taueimp=(3*eps0^2*me^2*(2*pi*qe/me*Te).^1.5)./(qe^4*Zimp.^2.*Nimp.*LneimpNRL);		
Tauimpe=(3*eps0^2*mimp^2*(2*pi*qe/mimp*Ti).^1.5)./(qe^4*Zimp.^2.*Ne.*LneimpNRL);
Tauiimp=(3*eps0^2*mi^2*(2*pi*qe/mi*Ti).^1.5)./(qe^4*Zi^2*Zimp.^2.*Nimp.*LniimpNRL);
Tauimpi=(3*eps0^2*mimp^2*(2*pi*qe/mimp*Ti).^1.5)./(qe^4*Zi^2*Zimp.^2.*Ni.*LniimpNRL);
Tauimpimp=(3*eps0^2*mimp^2*(2*pi*qe/mimp*Ti).^1.5)./(qe^4*Zimp.^4.*Nimp.*LnimpimpNRL);

% FORCES DE FRICTIONS 
%********************
xei=(me/mi*Ti./Te).^0.5;		xie=1./xei;
xeimp=(me/mimp*Ti./Te).^0.5;		ximpe=1./xeimp;
xiimp=(mi/mimp)^0.5*ones(size(xei));	ximpi=1./xiimp;

% FORMALISME DE HOULBERG ( signe - hors diagonale )
% Pour N11, on prend sqrt(Ta/Tb) pour que la condition 
% de symetrie lijab=ljiba soit verifiee.
%#####################################################
Mooaa=-(2)./(2)^1.5;
Mo1aa=1.5*(2)./(2).^2.5;
M1oaa=Mo1aa;
M11aa=-(13/4+4+15/2)./(2).^2.5;
Nooaa=-Mooaa;
No1aa=-Mo1aa;	

N1oaa=-M1oaa;
N11aa=27/4./(2).^2.5;

Mooei=-(1+me/mi)./(1+xei.^2).^1.5;
Mo1ei=1.5*(1+me/mi)./(1+xei.^2).^2.5;
M1oei=Mo1ei;
M11ei=-(13/4+4*xei.^2+15/2*xei.^4)./(1+xei.^2).^2.5;
Nooei=-Mooei;
No1ei=-xei.^2.*Mo1ei;
N1oei=-M1oei;
N11ei=27/4*sqrt(Te./Ti).*xei.^2./(1+xei.^2).^2.5;

Mooeimp=-(1+me/mimp)./(1+xeimp.^2).^1.5;
Mo1eimp=1.5*(1+me/mimp)./(1+xeimp.^2).^2.5;
M1oeimp=Mo1eimp;
M11eimp=-(13/4+4*xeimp.^2+15/2*xeimp.^4)./(1+xeimp.^2).^2.5;
Nooeimp=-Mooeimp;
No1eimp=-xeimp.^2.*Mo1eimp;
N1oeimp=-M1oeimp;
N11eimp=27/4*sqrt(Te./Ti).*xeimp.^2./(1+xeimp.^2).^2.5;

Mooie=-(1+mi/me)./(1+xie.^2).^1.5;
Mo1ie=1.5*(1+mi/me)./(1+xie.^2).^2.5;
M1oie=Mo1ie;
M11ie=-(13/4+4*xie.^2+15/2*xie.^4)./(1+xie.^2).^2.5;
Nooie=-Mooie;
No1ie=-xie.^2.*Mo1ie;
N1oie=-M1oie;
N11ie=27/4*sqrt(Ti./Te).*xie.^2./(1+xie.^2).^2.5;

Mooiimp=-(1+mi/mimp)./(1+xiimp.^2).^1.5;
Mo1iimp=1.5*(1+mi/mimp)./(1+xiimp.^2).^2.5;
M1oiimp=Mo1iimp;
M11iimp=-(13/4+4*xiimp.^2+15/2*xiimp.^4)./(1+xiimp.^2).^2.5;
Nooiimp=-Mooiimp;
No1iimp=-xiimp.^2.*Mo1iimp;
N1oiimp=-M1oiimp;
N11iimp=27/4*xiimp.^2./(1+xiimp.^2).^2.5;

Mooimpe=-(1+mimp/me)./(1+ximpe.^2).^1.5;
Mo1impe=1.5*(1+mimp/me)./(1+ximpe.^2).^2.5;
M1oimpe=Mo1impe;
M11impe=-(13/4+4*ximpe.^2+15/2*ximpe.^4)./(1+ximpe.^2).^2.5;
Nooimpe=-Mooimpe;
No1impe=-ximpe.^2.*Mo1impe;
N1oimpe=-M1oimpe;
N11impe=27/4*sqrt(Ti./Te).*ximpe.^2./(1+ximpe.^2).^2.5;

Mooimpi=-(1+mimp/mi)./(1+ximpi.^2).^1.5;
Mo1impi=1.5*(1+mimp/mi)./(1+ximpi.^2).^2.5;
M1oimpi=Mo1impi;
M11impi=-(13/4+4*ximpi.^2+15/2*ximpi.^4)./(1+ximpi.^2).^2.5;
Nooimpi=-Mooimpi;
No1impi=-ximpi.^2.*Mo1impi;
N1oimpi=-M1oimpi;
N11impi=27/4*ximpi.^2./(1+ximpi.^2).^2.5;

% l(ij,ab) / na.ma 
%#################
Looee=Mooaa./Tauee+Mooei./Tauei+Mooeimp./Taueimp + Nooaa./Tauee;
Lo1ee=Mo1aa./Tauee+Mo1ei./Tauei+Mo1eimp./Taueimp + No1aa./Tauee;
L1oee=M1oaa./Tauee+M1oei./Tauei+M1oeimp./Taueimp + N1oaa./Tauee;
L11ee=M11aa./Tauee+M11ei./Tauei+M11eimp./Taueimp + N11aa./Tauee;
Looii=Mooaa./Tauii+Mooie./Tauie+Mooiimp./Tauiimp + Nooaa./Tauii;
Lo1ii=Mo1aa./Tauii+Mo1ie./Tauie+Mo1iimp./Tauiimp + No1aa./Tauii;
L1oii=M1oaa./Tauii+M1oie./Tauie+M1oiimp./Tauiimp + N1oaa./Tauii;
L11ii=M11aa./Tauii+M11ie./Tauie+M11iimp./Tauiimp + N11aa./Tauii;
Looimpimp=Mooaa./Tauimpimp+Mooimpe./Tauimpe+Mooimpi./Tauimpi + ...
	Nooaa./Tauimpimp;
Lo1impimp=Mo1aa./Tauimpimp+Mo1impe./Tauimpe+Mo1impi./Tauimpi + ...
	No1aa./Tauimpimp;
L1oimpimp=M1oaa./Tauimpimp+M1oimpe./Tauimpe+M1oimpi./Tauimpi + ...
	N1oaa./Tauimpimp;
L11impimp=M11aa./Tauimpimp+M11impe./Tauimpe+M11impi./Tauimpi + ...
	N11aa./Tauimpimp;
Looei=Nooei./Tauei;
Lo1ei=No1ei./Tauei;
L1oei=N1oei./Tauei;
L11ei=N11ei./Tauei;
Looeimp=Nooeimp./Taueimp;
Lo1eimp=No1eimp./Taueimp;
L1oeimp=N1oeimp./Taueimp;
L11eimp=N11eimp./Taueimp;
	Looie=Nooie./Tauie;
	Lo1ie=No1ie./Tauie;
	L1oie=N1oie./Tauie;
	L11ie=N11ie./Tauie;
	Looiimp=Nooiimp./Tauiimp;
	Lo1iimp=No1iimp./Tauiimp;
	L1oiimp=N1oiimp./Tauiimp;
	L11iimp=N11iimp./Tauiimp;
Looimpe=Nooimpe./Tauimpe;
Lo1impe=No1impe./Tauimpe;
L1oimpe=N1oimpe./Tauimpe;
L11impe=N11impe./Tauimpe;
Looimpi=Nooimpi./Tauimpi;
Lo1impi=No1impi./Tauimpi;
L1oimpi=N1oimpi./Tauimpi;
L11impi=N11impi./Tauimpi;

% REGIME DE COLLISIONNALITE
%--------------------------
wee=(2*qe*Te./me).^0.5./(R0*qmag);
wii=(2*qe*Ti./mi).^0.5./(R0*qmag);
wimpimp=(2*qe*Ti./mimp).^0.5./(R0*qmag);
nuestar = 1./((epsilon+eps).^1.5.*wee.*Tauee);	% Cf Hirshman
nuistar = 1./((epsilon+eps).^1.5.*wii.*Tauii);	
nuimpstar = 1./((epsilon+eps).^1.5.*wimpimp.*Tauimpimp);
nuimpinorm = 1./((epsilon+eps).^1.5.*wii.*Tauimpi);

% VISCOSITES (Cf. fonction viscon.m de V. Basiuk & Article Kessel Nucl.Fus.1994 p.1221)
%****************************************************************************************
for i=1:length(Ne)
	tauxx=[Tauee(i) Tauii(i) Tauimpimp(i)]';
	m=[me mi mimp]';
	n=[Ne(i) Ni(i) Nimp(i)]';
	T=1e-3*[Te(i) Ti(i) Ti(i)]';
	Vt=(2*qe*[Te(i)./me Ti(i)./mi Ti(i)./mimp]').^0.5;
	nu=[nuestar(i) nuistar(i) nuimpstar(i)]';
	Z=[-1 Zi Zimp(i)]';
	nesp=3;
	epsK = epsilon(i);
	ft=FTRAP(i);
	[ni1(i,1) ni2(i,1) ni3(i,1)]=f_viscon(tauxx,m,n,T,Vt,1,2,nu,Z,nesp,epsK,ft);
end
muooi=ni1.*Tauii./(mi*Ni); 
muo1i=ni2.*Tauii./(mi*Ni);
mu1oi=muo1i; 
mu11i=ni3.*Tauii./(mi*Ni);

if ionregime == 1	% Daniel & Clemente extended ion regime formulas
	C0a = 1.5/(1+Ai/Aimp);
	ki = muo1i./muooi;
else
	C0a = 1.5/(1+Ai/Aimp);
	ki = muo1i./muooi;
end 

% l(ij,ab) / na.ma 
%#################
wca = qe*Zimp.*B0./mimp;
wci = qe*Zi*B0/mi;
nuswca = Looimpi./wca;
DPSa = 2*qmag_metrics.^2.*Looimpi./wca.^2.*B0.^2./B2avg.*qe.*Ti./(mimp);

UU = -Zimp./Zi.*(C0a+ki).*grad_ln_Ti;
GG = grad_ln_na - Zimp./Zi.*grad_ln_ni + (1+Zimp./Zi.*(C0a-1)).*grad_ln_Ti ;
GG0 = GG - grad_ln_na;

if geom == 0
		GG0_FH = -Zimp./Zi.*grad_ln_ni + Zimp./Zi.*(C0a-1).*grad_ln_Ti;
		
		UG = 1+UU./GG;
		Ae = nuswca.*FV./(R0.*B0).*qmag_metrics.^2./invaspct;
		AGe = Ae.*GG;
		AGe0 = Ae.*GG0;
		CDx0 = -epsilon./UG;
		QQ = CDx0.*(dNV./epsilon).*UU./GG;
		FF = CDx0.*(1-0.5*dNH./epsilon.*UU./GG);

		CDx = FF -0.5*(dphia-deltaM);
		CDVx = -0.5*(Dphia + QQ);

		RDx = sqrt((FF+0.5*(dphia-deltaM)).^2+0.25.*(Dphia-QQ).^2);
		DD = RDx.^2+AGe.^2.*(RDx./CDx0).^2;

		num = ((AGe./CDx0).^2-1).*(FF./(CDx0)+0.5*(dphia-deltaM)./CDx0)+AGe./CDx0.*(0.5*dNV.*UU./(epsilon.*GG)-0.5*Dphia./CDx0);
		cosa = RDx.*CDx0.*num./DD;

		num = 2.*AGe.*(FF./CDx0+0.5*(dphia-deltaM)./CDx0)+((AGe./CDx0).^2-1).*(0.5*Dphia-0.5*dNV.*CDx0.*UU./(epsilon.*GG));
		sina = RDx.*num./DD;

		alpha = acos(cosa);
		inda = find( sina < 0.);
		if length(inda)>0
			alpha(inda) = -alpha(inda);
		end

		dminx = CDx + RDx.*cosa;
		dmajx = CDVx + RDx.*sina;
		CD = CDx;
		CDV = CDVx;
		RD = RDx;
		RDV = RDx;

		facG = 1+dminx./epsilon+0.25*(dminx.^2+dmajx.^2)./epsilon.^2;
		facU = 0.5*(dminx-dNH)./epsilon+0.25*(dminx.^2+dmajx.^2-dminx.*dNH-dmajx.*dNV)./epsilon.^2;

		Vra_neo = -DPSa./amin.*( facG.*GG + facU.*UU )  ;

		Vra_cl = -DPSa./amin.*1./(2.*qmag.^2).*(1 + epsilon.*dminx + 2.*epsilon.^2).*GG  ;

		Vra = Vra_neo + Vra_cl ;

		dmin = dminx;
		dmaj = dmajx;
		Da_neo = DPSa.*facG;
		Da_cl = DPSa./(2.*qmag.^2).*(1 + epsilon.*dminx + 2.*epsilon.^2);
		Da = DPSa.*( facG + 1./(2.*qmag.^2).*(1 + epsilon.*dminx + 2.*epsilon.^2));
		Va_neo = - DPSa./amin.*( facG.*GG0 + facU.*UU ) ;
		Va_cl = - DPSa./amin.*1./(2.*qmag.^2).*(1 + epsilon.*dminx + 2.*epsilon.^2).*GG0;
		Va = - DPSa./amin.*( facG.*GG0 + facU.*UU ) ...
			- DPSa./amin.*1./(2.*qmag.^2).*(1 + epsilon.*dminx + 2.*epsilon.^2).*GG0;

		Vra0_neo = -DPSa.*GG./amin ;
		Vra0_cl = -DPSa./amin.*1./(2.*qmag.^2).*( GG ) ;
		Vra0 = Vra0_neo + Vra0_cl ;
		Da0 = DPSa.*( 1. + 1./(2.*qmag.^2).*(1 + 2.*epsilon.^2));
		Va0 = - DPSa./amin.*GG0 ...
			- DPSa./amin.*1./(2.*qmag.^2).*GG0 ...
			+  V_rot0_neo + V_rot0_cl;
		if pol_asym == 0
			Vra_neo = -DPSa.*GG./amin ;
			Vra_cl = -DPSa./amin.*1./(2.*qmag.^2).*( GG ) ;
			Vra = Vra_neo + Vra_cl ;
			Da = DPSa.*( 1. + 1./(2.*qmag.^2).*(1 + 2.*epsilon.^2));
			Va = - DPSa./amin.*GG0 ...
			- DPSa./amin.*1./(2.*qmag.^2).*GG0;
			dmin  = 0.*qmag;
			dmaj  = 0.*qmag;
		end

% Other ouputs of interest
		GGa = GG;
		GG0a = GG0;
		UUa = UU;

	out = struct('rhon',x, ...
	'DPSa',DPSa, ...
	'nuestar',nuestar, ...
	'nuistar',nuistar, ...
	'nuimpstar',nuimpstar, ...
	'nuimpinorm',nuimpinorm, ...
	'Vra',Vra, ...
	'Vra0',Vra0, ...
	'Vra_neo',Vra_neo, ...
	'Vra0_neo',Vra0_neo, ...
	'Vra_cl',Vra_cl, ...
	'Vra0_cl',Vra0_cl, ...
	'deltaM',deltaM, ...
	'dmin',dmin, ...
	'dmaj',dmaj, ...
	'Da',Da, ...
	'Da_neo',Da_neo, ...
	'Da_cl',Da_cl, ...
	'Va_neo',Va_neo, ...
	'Va_cl',Va_cl, ...
	'Va_neo_rot',Va_neo_rot, ...
	'Va_cl_rot',Va_cl_rot, ...
	'Va',Va, ...
	'Va_rot',Va_rot, ...
	'Da0',Da0, ...
	'Va0',Va0, ...
	'ki',ki, ...
	'C0a',C0a, ...
	'GGa',GGa, ...
	'GG0a',GG0a, ...
	'AGe0',AGe0, ...
	'UUa',UUa, ...
	'CD',CD, ...
	'CDV',CDV, ...
	'RD',RD, ...
	'RDV',RDV, ...
	'alpha',alpha);

else
% Poloidal asymmetry
% Convergence criterion: <b2/n>
%thetax = linspace(0,2*pi,length(BV(1,:))+1);
%theta = thetax(1:length(thetax)-1);
thetalong = [theta-2*pi,theta,theta+2*pi];
for it = 1:length(BV(1,:))
%	gradR2(:,it) = gradient(RV(:,it).^2,x+1.e-33);
	gradR2(:,it) = 0.*RV(:,it);
end

if length(Machi)==1
	for ix = 1:length(x)
		Factrot0(ix,1) = (mimp/mi).*(Machi./R0).^2;
		Factrot(ix,1) = (mimp/mi).*(Machi./R0).^2.*(1.-(Ai/Aimp).*Zimp(ix)./Zi);
	end
else
	for ix = 1:length(x)
		Factrot0(ix,1) = (mimp/mi).*(Machi(ix)./R0).^2;
		Factrot(ix,1) = (mimp/mi).*(Machi(ix)./R0).^2.*(1.-(Ai/Aimp).*Zimp(ix)./Zi);
	end
end

nn = ones(size(BV));
for ix = 1:length(x)
if pol_asym == 1
	b2 = BV(ix,:).^2./B2avg(ix,1);
	Apsi = jacob(ix,:).*FV(ix).*mimp*Looimpi(ix)./(qe*Zimp(ix))./((dpsidx(ix)).^2+1.e-33);
	ApsiX(ix,:) = Apsi;

	b2sNNavg = fluxavg(b2./NV(ix,:),jacob(ix,:));
%	gradR2avg = fluxavg(gradR2(ix,:),jacob(ix,:));
	gradR2avg = 0.;

	nnp = 2;
	nnx = nn(ix,:);
	ierr = 1;
	progx = prog;
	Erreur = 2*err;
	max_grad = 3;
	%dtheta = theta(2)-theta(1);
	matopt = 0;
	AA = zeros(length(theta)+1,length(theta)+1);
	LL = zeros(length(theta)+1,length(theta)+1);
	BB = zeros(length(theta),1);
	while ( Erreur>err & ierr<ierrmax)

		b2snavg = fluxavg(b2./nnx,jacob(ix,:));
		TermGG = 1. - b2./nnx./b2snavg;
		TermUU = b2./NV(ix,:) - b2sNNavg.*b2./nnx./b2snavg;
		TermRR = b2./nnx./b2snavg.*gradR2avg - gradR2(ix,:);
	
		FFF  = Apsi.*( GG(ix) + b2./NV(ix,:).*UU(ix) - Factrot(ix).*gradR2(ix,:));
		GGG = -Zimp(ix).*Te(ix)./Ti(ix).*(PhiV(ix,:)-PhiV(ix,1)) + Factrot0(ix).*(RV(ix,:).^2-RV(ix,1).^2);
		HHH = Apsi.*b2./b2snavg.*(GG(ix) + b2sNNavg.*UU(ix) - Factrot(ix).*gradR2avg);
		dGGG = interp1(thetalong,gradient([GGG,GGG,GGG],thetalong),theta,interpmeth);
		% AA*n = BB
		% Remplissage matrices AA, BB & LL
		for ii = 2:length(theta)-1
			dtheta = 0.5*(theta(ii+1) - theta(ii-1));
			AA(ii,ii-1) = -0.5/dtheta;
			AA(ii,ii) = -FFF(ii)-0.5/dtheta*(GGG(ii+1)-GGG(ii-1));
			AA(ii,ii+1) = 0.5/dtheta;

			LL(ii, ii-1) = 1.;
			LL(ii, ii) = -2.;
			LL(ii, ii+1) = 1.;

			BB(ii,1) = -HHH(ii);
		end
		ii = 1;
		dtheta = 0.5*(theta(2)-theta(end) +2*pi);
		AA(ii,ii) = -FFF(ii)-0.5/dtheta*(GGG(ii+1)-GGG(length(theta)-1));
		AA(ii,ii+1) = 0.5/dtheta;
		AA(ii,length(theta)) = -0.5/dtheta;
		LL(ii, ii) = -2.;
		LL(ii, ii+1) = 1.;
		LL(ii, length(theta)) = 1.;
		BB(ii,1) = -HHH(ii);

		ii = length(theta);
		dtheta = 0.5*(theta(1)-theta(end-1) + 2*pi);
		AA(ii,ii-1) = -0.5/dtheta;
		AA(ii,ii) = -FFF(ii)-0.5/dtheta*(GGG(1)-GGG(ii-1));
		AA(ii,ii+1) = 0.5/dtheta;
		LL(ii, ii-1) = 1.;
		LL(ii, ii) = -2.;
		LL(ii, ii+1) = 1.;
		BB(ii,1) = -HHH(ii);

		ii = length(theta)+1;
		dtheta = 0.5*(theta(2)-theta(end) + 2*pi);
		AA(ii,ii) = -1.;
		AA(ii,1) = 1;
		BB(ii,1) = 0.;

		CC = transpose(AA)*AA+regulweight*transpose(LL)*LL;
		nnya = inv(CC)*transpose(AA)*BB;
		nny = nnya(1:length(theta))'/fluxavg(nnya(1:length(theta))',jacob(ix,:));
		nny = max(nny, 1e-5);

		nnp = nnx;
		nnx = progx.*nnp + (1-progx).*nny;
		
		Erreur = max(abs(gradient(log(nnx)-GGG,theta)-FFF+HHH./nnx));
		Error(ix,ierr) = Erreur;

		ierr = ierr+1;
		%disp(['r/a=',num2str(x(ix),3),' - ierr=',int2str(ierr),' - err=',num2str(Erreur)])
		%plot(theta,nnx,theta,nnp,'--')
		%keyboard
		%pause
	end
	asym_error(ix,1) = Erreur;

	nn(ix,:) = nnx;
	dmin(ix,1) = 2.*mean((nn(ix,:)-1).*cos(theta));
	dmaj(ix,1) = 2.*mean((nn(ix,:)-1).*sin(theta));
else
	dmin(ix,1) = 0.*theta;
	dmaj(ix,1) = 0.*theta;
end
	% Radial flux
	FactV = amin./(dpsidx(ix).^2+1.e-33).*mimp*Looimpi(ix).*qe.*Ti(ix)./(qe*Zimp(ix)).^2.*FV(ix).^2./B2avg(ix,1);

	b2snavg = fluxavg(b2./nn(ix,:),jacob(ix,:));
	nsb2avg = fluxavg(nn(ix,:)./b2,jacob(ix,:));
	TermGGneo = 1./b2snavg - nsb2avg;
	TermUUneoa = fluxavg(b2./NV(ix,:),jacob(ix,:))./b2snavg;
	TermUUneob = fluxavg(nn(ix,:)./NV(ix,:),jacob(ix,:));
	TermUUneo = TermUUneoa - TermUUneob;

	Vra_neo(ix,1) = FactV.*( TermGGneo.*GG(ix) + TermUUneo.*UU(ix)  );
	Va_neo(ix,1) = FactV.*( TermGGneo.*GG0(ix) + TermUUneo.*UU(ix)  );
	Da_neo(ix,1) = -amin*FactV.*TermGGneo;

	TermGGcl = B2avg(ix,1)./FV(ix).^2.*fluxavg(nn(ix,:).*RV(ix,:).^2,jacob(ix,:))-nsb2avg;
	Vra_cl(ix,1) = -FactV.*( TermGGcl.*GG(ix)  );
	Va_cl(ix,1) = -FactV.*( TermGGcl.*GG0(ix)  );
	Da_cl(ix,1) = amin.*FactV.*TermGGcl;

	Vra(ix,1) = Vra_neo(ix,1) + Vra_cl(ix,1);
	Va(ix,1) = Va_neo(ix,1) + Va_cl(ix,1);
	Da(ix,1) = Da_neo(ix,1) + Da_cl(ix,1);

	% For test
	%test(ix,1) = TermGGneo;
	%disp([' TermGGcl = ',num2str(TermGGcl),' - nsb2avg = ',num2str(nsb2avg)])
	% End test
end

	out = struct('nn',nn, ...
	'rhon',x, ...
	'dmin',dmin, ...
	'dmaj',dmaj, ...
	'Vra_neo',Vra_neo, ...
	'Va_neo',Va_neo, ...
	'Da_neo',Da_neo, ...
	'Vra_cl',Vra_cl, ...
	'Va_cl',Va_cl, ...
	'Da_cl',Da_cl, ...
	'Vra',Vra, ...
	'Va',Va, ...
	'Da',Da, ...
	'asym_error',asym_error, ...
	'Error',Error);
end

return

%====================================================================================
%   FUNCTION F_VISCON
function [nui1,nui2,nui3]=f_viscon(tau,m,n,T,vt,k,i,nustar,Z,nesp,epsK,ft)
% F_VISCON : coefficient de viscosité de l'espèce i pour la couche k
%
% Repris de viscon.m (V. Basiuk)
%
% Syntaxe  :  [nui1,nui2,nui3]=F_viscon(tau,m,n,T,vt,k,i,nustar,Z,nesp,espK,ft);
%
% masse : m, densite n, Temperature T
%
% k: numero de la couche
%
% Programmes appeles :
%   chand.m
%
% Auteur : V. Basiuk
% Date   : 28/11/2000
% 

np        = 201;
x         = linspace(0.0001,4,np);
dx        = x(2)-x(1);
fc        = 1-ft;
nuaDv     = 0;
nuasv     = 0;
nuapv     = 0;

for j=1:nesp
  xj        = x*vt(i,k)/vt(j,k);
  fact1     = n(j,k)/n(i,k)*Z(j)^2/Z(i)^2;
  F1        = (erf(xj)-chand(xj))./x.^3;
  fact2     = 3*sqrt(pi)/4/tau(i,k)*fact1;
  nuaDv     = nuaDv+fact2.*F1;
  nuasv     = nuasv+fact2*2*(T(i,k)/T(j,k) + (vt(i,k)/vt(j,k))^2)*chand(xj)./x;
  nuapv     = nuapv+fact2.*2*chand(xj)./x.^3;
end
nuaEv     = 2*nuasv - 2*nuaDv - nuapv;
nuatv     = 3*nuaDv + nuaEv;
F2        = 1 + 2.48*nustar(i,k)*nuaDv*tau(i,k)./x;
F3        = 1 + 1.9634*nustar(i,k)*nuatv*tau(i,k)*epsK^1.5./x;
nuatotv   = nuaDv./F2./F3;
int11     = nuatotv*tau(i,k).*x.^4..*exp(-x.^2);
int12     = nuatotv*tau(i,k).*x.^6..*exp(-x.^2);
int22     = nuatotv*tau(i,k).*x.^8..*exp(-x.^2);
coef      = ones(size(x));
pair      = 2:2:np;
coef(pair)= 2*ones(size(1:(np-1)/2));
coef(1)   = 0.5;
coef(np)  = 0.5;
Av11(k)   = coef*int11'*2/3*dx;
Av12(k)   = coef*int12'*2/3*dx;
Av22(k)   = coef*int22'*2/3*dx;
fact3     = ft(k)/fc(k)*n(i,k)*m(i)/tau(i,k)*8/3/sqrt(pi);
Ki11      = fact3*Av11(k);
Ki21      = fact3*Av12(k);
Ki22      = fact3*Av22(k);
nui1      = Ki11;
nui2      = Ki21 - 2.5*Ki11;
nui3      = Ki22 - 5*Ki21 + 25/4*Ki11;
return

%--------------------------------------------
function y=chand(x);
% CHAND : Fonction de Chandrasekar
%
% Syntaxe : y=chand(x)
% 
% Auteur : V. Basiuk

derfx=2./sqrt(pi).*exp(-x.^2);
y=(erf(x)-x.*derfx)./(2*x.^2);
return

%--------------------------------------------
function Aavg = fluxavg(AF,JJ);
thetay = linspace(0,1,length(AF(1,:)));
for ix=1:length(AF(:,1))
	denom = trapz(thetay,JJ(ix,:));
	if abs(denom)>0
		Aavg(ix,1) = trapz(thetay,AF(ix,:).*JJ(ix,:))./trapz(thetay,JJ(ix,:));
	else
		Aavg(ix,1) = mean(AF(ix,:));
	end
end
return

