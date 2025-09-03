# -*- coding: utf-8 -*-
''' This routine runs FACIT code from IMAS data via input.dat 
    Inputs: shot, time, [vacuum, symmetric]
    Author : P. Maget (patrick.maget@cea.fr)
    Date   : 01 / 08 / 2022
'''
# Standard python modules
from __future__ import (unicode_literals, absolute_import, \
                        print_function, division)
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
#import matplotlib.ticker as tick
import numpy as np
import pywed as pw
import imas
from ypstruct import *	# to have a structured output
from scipy.io import savemat
import os
try:
    from imas import imasdef
except ImportError as err:
    from imas.backends.imas_core.imas_interface import imasdef

# facit options
optss=0
t1 = 0.
t2 = 1e4
x1 = 0.41
x2 = 0.49
geo = 0
Zeff = 0
optphinat=0
optasym=0
rotation = 1
ionregime = 1

# IMAS IDS
shot = 58248
user = 'PM156629'
machine= 'west'
run = 1
occ = 0

print(' ')
print('shot           =', shot)
print('machine        =', machine)
print('user           =', user)
print('run            =', run)
print('occ            =', occ)

# Constants
mu0 = 4*np.pi*1e-7
qe = 1.602e-19
eps0=8.8542e-12
me=9.1096e-31
mp=1.6726e-27

def F_Ionisation_stat_W(Te_in):
	te_ref = [0,   30,40,  50,  60,  70,  100, 150, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000, 2300, 2700, 3000, 3500, 4000, 5000, 6000, 7000, 10000, 12000, 15000, 20000, 25000, 30000, 40000, 100000]
	z_ref =  [0,  6.77,  7.95,  9.13,  10.25,  11.46, 13.47, 15.70, 17.71, 20.18, 21.62, 22.81, 23.75, 25.33, 26.67, 30.47, 33.93, 36.08, 40.39, 42.84, 44.43, 45.30, 46.72, 48.20, 50.00, 54.58, 57.02, 59.42, 62.02, 63.64, 64.99, 66.83, 74]

	zave  = np.interp(Te_in,te_ref,z_ref)
	zave  = np.amin([74*np.ones(len(Te_in)),np.amax([0*zave,zave],axis=0)],axis=0)
	return  zave

#%Check FACIT interface is present in current directory
FDIR = '/Home/PM156629/facit/fortran_src'
dum = os.path.exists('FACIT_interface')
if dum!=1:
	cmd = '! cp '+FDIR+'/FACIT_interface .'
	os.system(cmd)

# Default charge and mass of light impurity when Zeff is imposed (>1)
Az = 14
Zz = 7

# Parameters for natural asymmetry
bb = -1.65
phif0 = 0.1

# Open shot and specific run of machine
backend_id = imasdef.HDF5_BACKEND
idx = imas.DBEntry(backend_id,machine,shot,run,user, data_version='3')
idx.open()
equil = idx.get('equilibrium',occ)
cp = idx.get('core_profiles',0)
idx.close()

itermax = 10
regulweight = 1e0
param = [1.e-3, 0.5, regulweight, itermax]

maxiterWconv = 40
errmax = 1.e-2
mix = 0.2

nth = 64
thetax = np.linspace(0,2*np.pi,nth+1)
theta = thetax[0:nth]

i_run = np.squeeze(np.argwhere((cp.time>t1) & (cp.time<t2)))
ntime = len(cp.time[i_run])
time_out = cp.time[i_run]

nx = len(cp.profiles_1d[0].grid.psi)
n_ions = len(cp.profiles_1d[0].ion)
psi = np.zeros((ntime,nx))
dpsidx = np.zeros((ntime,nx))
psin = np.zeros((ntime,nx))
rhon = np.zeros((ntime,nx))
qpsi = np.zeros((ntime,nx))
Ne = np.zeros((ntime,nx))
Ni = np.zeros((nx,n_ions))
Ni_all = np.zeros(nx)
Te = np.zeros((ntime,nx))
Zi = np.zeros((nx,n_ions))
Ti = np.zeros((ntime,nx))
qmin = np.zeros(ntime)
Ip = np.zeros(ntime)
Te0 = np.zeros(ntime)
Ti0= np.zeros(ntime)
Ne0 = np.zeros(ntime)
F = np.zeros(nx)
Aion = np.zeros(n_ions)
Zion = np.zeros(n_ions)
nion0 = np.zeros(n_ions)

k = i_run[int(ntime/2)]

timecp = time_out[i_run[k]]
print('Time in core_profiles IDS =', timecp)
psi[k,:] = cp.profiles_1d[i_run[k]].grid.psi
rhon_in = cp.profiles_1d[i_run[k]].grid.rho_tor_norm
rhon[k,:] = rhon_in
		
kimas = np.min(np.argwhere(equil.time>=timecp))
x_imas = equil.time_slice[kimas].profiles_1d.rho_tor_norm
q_in = equil.time_slice[kimas].profiles_1d.q
amin = equil.time_slice[kimas].boundary.minor_radius
Rgeo = equil.time_slice[kimas].boundary.geometric_axis.r
Rmag = equil.time_slice[kimas].global_quantities.magnetic_axis.r
invaspct = amin/Rgeo
EPSILON = amin/Rgeo*x_imas
Bmag = equil.time_slice[kimas].global_quantities.magnetic_axis.b_field_phi
Bgeo = (Rmag*Bmag/Rgeo)
if geo==1:
	dum = len(equil.time_slice[kimas].profiles_2d)
	if dum>=2:
		theta_imas = equil.time_slice[kimas].profiles_2d[1].theta
		nth_imas = len(theta_imas[0,:])
		jacob_imas = equil.time_slice[kimas].coordinate_system.jacobian
		if len(jacob_imas)==0:
			jacob_imas = np.zeros((len(x_imas),nth_imas))
			for ix in range(len(x_imas)):
				jacob_imas[ix,:] =  amin**2*Rgeo*x_imas[ix]
		phi_imas = equil.time_slice[kimas].profiles_1d.phi
		psi_imas = equil.time_slice[kimas].profiles_1d.psi
		bphi_imas = equil.time_slice[kimas].profiles_2d[1].b_field_phi
		br_imas = equil.time_slice[kimas].profiles_2d[1].b_field_r
		bz_imas = equil.time_slice[kimas].profiles_2d[1].b_field_z
		R_imas = equil.time_slice[kimas].profiles_2d[1].r
		Z_imas = equil.time_slice[kimas].profiles_2d[1].z
		F_imas = bphi_imas[:,0]*R_imas[:,0]
		B_imas = (bphi_imas**2 + br_imas**2 + bz_imas**2)**0.5
		jacob_dum = np.zeros((len(x_imas),len(theta)))
		B_dum = np.zeros((len(x_imas),len(theta)))
		R_dum = np.zeros((len(x_imas),len(theta)))
		Z_dum = np.zeros((len(x_imas),len(theta)))
		for ix in range(len(x_imas)):
			if ix == 0:
				theta_vect = np.concatenate((theta_imas[ix+1,1:nth_imas-1]-2*np.pi,theta_imas[ix+1,1:nth_imas-1],theta_imas[ix+1,2:nth_imas]+2*np.pi))
			else:
				theta_vect = np.concatenate((theta_imas[ix,1:nth_imas-1]-2*np.pi,theta_imas[ix,1:nth_imas-1],theta_imas[ix,2:nth_imas]+2*np.pi))
			jacob_vect = np.concatenate((jacob_imas[ix,1:nth_imas-1],jacob_imas[ix,1:nth_imas-1],jacob_imas[ix,2:nth_imas]))
			B_vect = np.concatenate((B_imas[ix,1:nth_imas-1],B_imas[ix,1:nth_imas-1],B_imas[ix,2:nth_imas]))
			R_vect = np.concatenate((R_imas[ix,1:nth_imas-1],R_imas[ix,1:nth_imas-1],R_imas[ix,2:nth_imas]))
			Z_vect = np.concatenate((Z_imas[ix,1:nth_imas-1],Z_imas[ix,1:nth_imas-1],Z_imas[ix,2:nth_imas]))
			isort = np.argsort(theta_vect)
			jacob_dum[ix,:] = np.interp(theta,theta_vect[isort],jacob_vect[isort])
			B_dum[ix,:] = np.interp(theta,theta_vect[isort],B_vect[isort])
			R_dum[ix,:] = np.interp(theta,theta_vect[isort],R_vect[isort])
			Z_dum[ix,:] = np.interp(theta,theta_vect[isort],Z_vect[isort])
		
		F = np.interp(rhon_in,x_imas,np.abs(F_imas))
		psi_dum = np.interp(rhon_in,x_imas,psi_imas)
		dpsidx = np.gradient(-psi_dum/(2*np.pi),rhon_in)
		dpsidx_imas = np.interp(x_imas,rhon_in,dpsidx)
		jacob = np.zeros(len(rhon_in),len(theta))
		B = np.zeros((len(rhon_in),len(theta)))
		R = np.zeros((len(rhon_in),len(theta)))
		Z = np.zeros((len(rhon_in),len(theta)))
		for itheta in range(len(theta)):
			jacob[:,itheta]  = 2*np.pi*np.interp(rhon_in,x_imas,jacob_dum[:,itheta] *dpsidx_imas)
			B[:,itheta]  = np.interp(rhon_in,x_imas,B_dum[:,itheta] )
			R[:,itheta]  = np.interp(rhon_in,x_imas,R_dum[:,itheta] )
			Z[:,itheta]  = np.interp(rhon_in,x_imas,Z_dum[:,itheta] )
	else:
		print('Geometry option geo=1 not possible with the current input IDS')
		geo = 0
		dpsidx = amin**2*Bgeo*rhon_in/q_in
		nth = 64
		thetax = np.linspace(0,2*np.pi,nth+1)
		theta = thetax[0:nth]
		R = np.zeros((nx,nth))
		Z = np.zeros((nx,nth))
		B = np.zeros((nx,nth))
		jacob = np.zeros((nx,nth))
		for ix in range(len(rhon_in)):
			jacob[ix,:] = amin**2*Rgeo*rhon_in[ix]
			F[ix] = Rgeo*Bgeo
			R[ix,:] = Rgeo*(1+EPSILON[ix]*np.cos(theta))
			B[ix,:] = Bgeo/(1+EPSILON[ix]*np.cos(theta))*(1+(EPSILON[ix]**2/q_in[ix]**2))**0.5
else:
	dpsidx = amin**2*Bgeo*rhon_in/q_in
	F = np.zeros(nx)
	R = np.zeros((nx,nth))
	B = np.zeros((nx,nth))
	jacob = np.zeros((nx,nth))
	for ix in range(len(rhon_in)):
		R[ix,:] = Rgeo* (1.+ EPSILON[ix]*np.cos(theta))
		B[ix,:] = Rgeo*Bgeo/R[ix,:]*(1  + (EPSILON[ix]/q_in[ix])**2)**0.5
		jacob[ix,:] = amin**2*rhon_in[ix]*R[ix,:]
		F[ix] = Rgeo*Bgeo

ne_in = cp.profiles_1d[i_run[k]].electrons.density
Te_in     = cp.profiles_1d[i_run[k]].electrons.temperature
	
Ti_in = cp.profiles_1d[i_run[k]].t_i_average
ni_in = cp.profiles_1d[i_run[k]].n_i_thermal_total
zeff_in = cp.profiles_1d[i_run[k]].zeff

Aieff = 0.
Zieffx = 0.*rhon_in
neff = 0.*rhon_in
for ii in range(n_ions):
	Ni[:,ii] = cp.profiles_1d[i_run[k]].ion[ii].density
	Zi[:,ii] = cp.profiles_1d[i_run[k]].ion[ii].z_ion_1d
	Aion[ii] = cp.profiles_1d[i_run[k]].ion[ii].element[0].a
	Zion[ii] = cp.profiles_1d[i_run[k]].ion[ii].element[0].z_n
	nion0[ii] = Ni[0,ii]
	Zieffx = Zieffx + Zion[ii]**2*Ni[:,ii]
	Aieff = Aieff + Aion[ii]*Ni[:,ii] 
	neff = neff + Zion[ii]*Ni[:,ii] 
Zieff = np.mean(Zieffx/neff)
neff = neff/Zieff
Aieff = Zieff*np.mean(Aieff/ne_in)
mieff = mp*Aieff
imain = np.argwhere(nion0==max(nion0))[0][0]

factPhi = EPSILON/(1+zeff_in*Te_in/Ti_in)

if len(cp.profiles_1d[i_run[k]].rotation_frequency_tor_sonic)>0:
	vtor_in = cp.profiles_1d[i_run[k]].rotation_frequency_tor_sonic*(2*np.pi*Rgeo)
else:
	vtor_in = 0*rhon_in
Machi = (0.5*mp*Aieff*vtor_in**2/(qe*Ti_in))**0.5

if max(Aion)<100:
	flagW = 0
	cW = 1e-5
	print('Add tungsten with concentration cW0 = '+str(cW))
else:
	flagW = 1
if Zeff>1:
	Zieff = Zeff;
	Aieff = 2*Zieff;
Zimp = F_Ionisation_stat_W(Te_in)
if flagW == 1:
	Aimp = max(Aion)
	iw = np.argwhere(Aion==Aimp)[0][0]
	nimp_in = cp.profiles_1d[i_run[k]].ion[iw].density
else:
	nimp_in = cW*ne_in
	Aimp = 184

AsymPhi = np.zeros((nx,2))
AsymN = np.zeros((nx,2))
PhiV = np.zeros((nx,nth))
NV = np.zeros((nx,nth))
		
TperpsTpar = np.ones(len(rhon_in))
fH = 0.0
bC = 1.

if optphinat == 1:
	wci = qe*Zieff*Bgeo/mieff
	LniiNRL = 23-np.log(Zieff**2*(2.*(ni_in/1e6)**0.5*Zieff**2))+1.5*np.log(Ti_in)
	Tauii = (3*eps0**2*mieff**2*(2*np.pi*qe/mieff*Ti_in)**1.5)/(qe**4*Zieff**4*ni_in*LniiNRL)	
	grad_ln_Ti = np.gradient(Ti_in,rhon_in)/Ti_in
	AsymPhinat = bb/(Zieff+Ti_in/Te_in)*q_in**2*rhon_in/((amin/Rgeo*rhon_in)**2.5)*grad_ln_Ti/(Tauii*wci)*phif0
else:
	AsymPhinat = 0*rhon_in

AsymPhi[:,0] = factPhi*(fH*(TperpsTpar-1)*bC/(bC+TperpsTpar*(1-bC)) + 2.*Machi**2)
AsymPhi[:,1] = AsymPhinat
AsymN[:,0] = - Zieff*Te_in/Ti_in*AsymPhi[:,0] + 2.*EPSILON*Machi**2
AsymN[:,1] = - Zieff*Te_in/Ti_in*AsymPhi[:,1]
for ix in range(len(rhon_in)):
	PhiV[ix,:] = AsymPhi[ix,0]*np.cos(theta)+AsymPhi[ix,1]*np.sin(theta)
	NV[ix,:] = 1+AsymN[ix,0]*np.cos(theta)+AsymN[ix,1]*np.sin(theta)

if Zeff<1:
	ind_impurities = np.argwhere((Aion<Aimp) & (nion0>0))
else:
	ind_impurities = imain

indx = np.squeeze(np.argwhere((rhon_in>x1) & (rhon_in<x2)))

for za in range(len(ind_impurities)):
	if Zeff<1:
		Ai = cp.profiles_1d[i_run[k]].ion[ind_impurities[za][0]].element[0].a
		Zi = cp.profiles_1d[i_run[k]].ion[ind_impurities[za][0]].element[0].z_n
		ni_in = cp.profiles_1d[i_run[k]].ion[ind_impurities[za][0]].density
	else:
		ni_in = ne_in/Zeff
		Ai = 2*Zeff
		Zi = Zeff
	x_Nimp = nimp_in[indx]
	dum = np.gradient(nimp_in,rhon_in)
	x_gradNimp = dum[indx]
	x_Zimp = Zimp[indx]

	x_Ni = ni_in[indx]
	dum = np.gradient(ni_in,rhon_in)
	x_gradNi = dum[indx]
	x_AsymN = AsymN[indx,:]
	x_AsymPhi= AsymPhi[indx,:]
	x_NV = np.zeros((len(indx),len(theta)))
	x_PhiV = np.zeros((len(indx),len(theta)))
	for ith in range(len(theta)):
		x_NV[:,ith] = np.squeeze(NV[indx,ith])
		x_PhiV[:,ith] = np.squeeze(PhiV[indx,ith])
	x_x = rhon_in[indx]
	x_dpsidx = dpsidx[indx]
	x_FV = F[indx]
	if geo == 1:
		x_BV = B[indx,:]
		x_RV = R[indx,:]
		x_jacob = jacob[indx,:]
	x_Te  = Te_in[indx]
	x_Ti =	Ti_in[indx]
	x_gradTi = np.gradient(x_Ti,x_x)
	x_Ne =ne_in[indx]
	x_qmag = np.abs(q_in[indx])
	x_Machi = Machi[indx]
	invaspct = amin/Rgeo
	details = 0

	VS = np.zeros(max(18,nth))
	VS[0:15] = [nx, nth, Aimp,Ai,Zi,Bgeo,Rgeo,invaspct,optasym,geo,rotation,param[0],param[1],param[2],param[3]]
	VV = [VS]
	for ii in range(len(indx)):
		VS = np.zeros(max(18,nth))
		VS[0:18] = [ x_x[ii],x_Zimp[ii],x_Te[ii],x_Ti[ii],x_Ne[ii],x_Ni[ii],x_Nimp[ii],x_qmag[ii],x_Machi[ii],x_gradNi[ii],x_gradTi[ii],x_gradNimp[ii],x_FV[ii],x_dpsidx[ii],x_AsymPhi[ii,0],x_AsymPhi[ii,1],x_AsymN[ii,0],x_AsymN[ii,1]]
		VV.append(VS)
	if geo>0:
		VS = np.zeros(max(18,nth))
		for ii in range(nth):
			VS[ii] = theta[ii]
		VV.append(VS)
		for ii in range(len(indx)):
			VS = np.zeros(max(18,nth))
			VS[0:nth] = x_BV[ii,:]
			VV.append(VS)
		for ii in range(len(indx)):
			VS = np.zeros(max(18,nth))
			VS[0:nth] = x_RV[ii,:]
			VV.append(VS)
		for ii in range(len(indx)):
			VS = np.zeros(max(18,nth))
			VS[0:nth] = x_jacob[ii,:]
			VV.append(VS)
		for ii in range(len(indx)):
			VS = np.zeros(max(18,nth))
			VS[0:nth] = x_PhiV[ii,:]
			VV.append(VS)
		for ii in range(len(indx)):
			VS = np.zeros(max(18,nth))
			VS[0:nth] = x_NV[ii,:]
			VV.append(VS)
	WW = np.array(VV)
	np.savetxt('facit_input.dat',WW,fmt='%14.7e', delimiter='  ')

	cmd = './FACIT_interface'
	os.system(cmd)

	np.loadtxt('facit_output.dat', delimiter='  ')



