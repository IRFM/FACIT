#Example of a new interface to run the code.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import h5py




#Part 1 Prepare the inputs:

#numbers
    #integers
nx = 1000 #size of radial arrays
nth = 1000 #size of poloidal arrays
nions = 5  #number of ion species
Aa = 184

    #reals
invaspct = 0.2
B0= 3.7
R0 = 2.5
qa=4


#Vectors


xn = np.arange(0, 1, nx)
theta = np.arange(0, 2*np.pi, nth)
Za = np.zeros(nx)
Ai = np.zeros(nions)
Te = np.zeros(nx)
Ti = np.zeros(nx)
Ne = np.zeros(nx)
Na = np.zeros(nx)
Machi = np.zeros(nx)
gradTi = np.zeros(nx)
gradNa = np.zeros(nx)
qmag = np.zeros(nx)
dpsidx = np.zeros(nx)
FV = np.zeros(nx)
BV = np.zeros(nx)
RV = np.zeros(nx)
regulopt = np.zeros(4)


# Matrix.

Zi = np.zeros((nx, nions))
Ni = np.zeros((nx, nions))
gradNi = np.zeros((nx, nions))
jacob = np.zeros((nx, nth))
AsymPhi = np.zeros((nx, nth))

#Tensor.

AsymN =  np.zeros((nx, nth, nions))


#logical

pol_asym = True
rotation = True
solution = 3
full_geom = True

#extra 
ni0 = [4e19, 1e19, 1e14, 5e9, 2e10] #select inicial density for each case


#Profiles

Ai = [2,3,4,14,16]

Ni[:,0] = ni0[0] * (1 - xn**2)**2 #Deuterium
Ni[:,1] = ni0[1]* (1 - xn**2)**2 #Helium or T
Ni[:,2] = ni0[2]* (1 - xn**2)**2 #Helium or T
Ni[:,3] = ni0[3]* (0.3 + np.exp(6*(xn-1))) #light impurities
Ni[:,4] = ni0[4]* (0.3 + np.exp(6*(xn-1))) #light impurities

Na = 6e12 * (0.3 + np.exp(6*(xn-1)))

Ne = 6e19 * (1 - xn**2)**2

Zi[:,0] = 1 #Deuterium
Zi[:,1] = 1 # T
#definir ifs
if xn>=0.6:
    Zi[:,2] = 1 #Helium
else:
    Zi[:,2] = 2

Zi[:,3] = 6 - 5 * np.exp(-((xn - 1) / 0.25)**2)
Zi[:,4] = 7 - 6 * np.exp(-((xn - 1) / 0.25)**2)



Zimp = 60 - 55 * np.exp(-((xn - 1) / 0.3)**2)


qmag =1 + (qa-1)*xn**2

Btheta = 0.5 * (xn / qmag)
Btheta[0] = 0.0 


Ti = 1.5*(1-xn**2)**2

Te = 3 *(1-xn**2)**2

#Gradients

gradTi = np.grad(Ti,xn)
gradNa = np.grad(Na,xn)

for i in range(len(xn)):
    gradNi[:,i]= np.grad(Ni[:,i],xn)

for i in range(len(FV)):
    FV[i]= R0*B0

#BV=
#RV=
#jacob= 
#AsymPhi=
#AsymN= 


#Part 2 Comunicate with fortran (We can use f2py)

import facit_code

#should be enough with the f2py but to be determined.



#Part 3 Save the info.


class FACITResult:
    def __init__(self, outputs):

        (
            self.Flux_imp, self.Da, self.Da_M, self.Vconv,
            self.nn, self.dmin, self.dmaj,

            self.Da_PS_M, self.Da_BP_M, self.Da_CL_M,
            self.Ka_PS_M, self.Ka_BP_M, self.Ka_CL_M,
            self.Ha_PS_M, self.Ha_BP_M, self.Ha_CL_M,
            self.Va_PS_M, self.Va_BP_M, self.Va_CL_M,

            self.Da_PS, self.Da_BP, self.Da_CL,
            self.Ka_PS, self.Ka_BP, self.Ka_CL,
            self.Ha_PS, self.Ha_BP, self.Ha_CL,
            self.Va_PS, self.Va_BP, self.Va_CL
        ) = outputs

        # Dictionary per category.
        self.Da_comp = {
            "PS": self.Da_PS,
            "BP": self.Da_BP,
            "CL": self.Da_CL,
            "PS_M": self.Da_PS_M,
            "BP_M": self.Da_BP_M,
            "CL_M": self.Da_CL_M,
        }

        self.Ka_comp = {
            "PS": self.Ka_PS,
            "BP": self.Ka_BP,
            "CL": self.Ka_CL,
            "PS_M": self.Ka_PS_M,
            "BP_M": self.Ka_BP_M,
            "CL_M": self.Ka_CL_M,
        }

        self.Ha_comp = {
            "PS": self.Ha_PS,
            "BP": self.Ha_BP,
            "CL": self.Ha_CL,
            "PS_M": self.Ha_PS_M,
            "BP_M": self.Ha_BP_M,
            "CL_M": self.Ha_CL_M,
        }

        self.Va_comp = {
            "PS": self.Va_PS,
            "BP": self.Va_BP,
            "CL": self.Va_CL,
            "PS_M": self.Va_PS_M,
            "BP_M": self.Va_BP_M,
            "CL_M": self.Va_CL_M,
        }

raw_outputs = facit_code.facit( nx, nth, xn, nions, theta, Za, Aa, Zi, Ai,
    Te, Ti, Ne, Ni, Na, Machi, gradTi, gradNi, gradNa, invaspct, B0, R0, qmag, dpsidx,
    FV, BV, RV, jacob, AsymPhi, AsymN, pol_asym, full_geom, solution, rotation, regulopt)

result = FACITResult(raw_outputs)

def to_hdf5(self, filename="facit_output.h5"):
    with h5py.File(filename, "w") as f:
        f["Flux_imp"] = self.Flux_imp
        f["Da"] = self.Da
        f["Da_M"] = self.Da_M
        f["Vconv"] = self.Vconv
        f["nn"] = self.nn
        f["dmin"] = self.dmin
        f["dmaj"] = self.dmaj

        grp_Da = f.create_group("Da_comp")
        for k, v in self.Da_comp.items():
            grp_Da[k] = v
        
        grp_Ka = f.create_group("Ka_comp")
        for k, v in self.Ka_comp.items():
            grp_Ka[k] = v
        
        grp_Ha = f.create_group("Ha_comp")
        for k, v in self.Ha_comp.items():
            grp_Ha[k] = v
        
        grp_Va = f.create_group("Va_comp")
        for k, v in self.Va_comp.items():
            grp_Va[k] = v

#Part 3 Plots.



def plot_profile(x, y, xlabel="", ylabel="", title=""):
    plt.figure(figsize=(6,4))
    sns.lineplot(x=x, y=y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plot_profile( x=xn, y=result.Da, xlabel="r", ylabel="Da", title="Particle Diffusivity Da(x)")
plot_profile( x=result.dmin, y=result.dmaj, xlabel="delta h", ylabel="Delta V", title="Assymetry")

