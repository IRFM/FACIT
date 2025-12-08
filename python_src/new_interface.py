import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
import facitpy

# ------------------------------
# Part 1: Inputs and profiles
# ------------------------------

nx = 100
nth = 100
nions = 5
Aa = 184

mp = 1.67262192369e-27
q_e = 1.602176634e-19
invaspct = 0.2
B0 = 3.7
R0 = 2.5
qa = 4
ma = Aa*mp
cor = 1e-33

xn = np.linspace(0,1,nx)
theta = np.linspace(0,2*np.pi,nth)

Za = np.zeros(nx)
Ai = np.array([2,3,4,14,16])
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
regulopt = np.array([0,0,0,0],dtype=np.float64)

BV = np.zeros((nx,nth))
RV = np.zeros((nx,nth))
Zi = np.zeros((nx,nions))
Ni = np.zeros((nx,nions))
gradNi = np.zeros((nx,nions))
jacob = np.zeros((nx,nth))
AsymPhi = np.zeros((nx,2))
AsymN = np.zeros((nx,nions,2))

pol_asym = True
rotation = True
solution = 3
full_geom = True

ni0 = [4e19,1e19,1e14,5e9,2e10]

Ai = np.array([2,3,4,14,16])

Ni[:,0] = ni0[0]*(1-xn**2)**2
Ni[:,1] = ni0[1]*(1-xn**2)**2
Ni[:,2] = ni0[2]*(1-xn**2)**2
Ni[:,3] = ni0[3]*(0.3+np.exp(6*(xn-1)))
Ni[:,4] = ni0[4]*(0.3+np.exp(6*(xn-1)))

Na = 6e12*(0.3+np.exp(6*(xn-1)))
Ne = 6e19*(1-xn**2)**2

Zi[:,0] = 1
Zi[:,1] = 1

for i in range(nx):
    Zi[i,2] = 1 if xn[i]>=0.6 else 2
    FV[i] = R0*B0

Zi[:,3] = 6 - 5*np.exp(-((xn-1)/0.25)**2)
Zi[:,4] = 7 - 6*np.exp(-((xn-1)/0.25)**2)

Za = 60 - 55*np.exp(-((xn-1)/0.3)**2)

qmag = 1 + (qa-1)*xn**2
Ti = 1.5*(1-xn**2)**2
Te = 3*(1-xn**2)**2

for j in range(nth):
    RV[:,j] = xn*invaspct*np.cos(theta[j])
    BV[:,j] = FV**2/(RV[:,j]**2+cor)*(1+invaspct**2/(RV[:,j]**2+cor)*xn**2/(qmag**2+cor))
    jacob[:,j] = xn / (B0*(1-invaspct*np.cos(theta[j])) + cor)

S1 = 0
S2 = Na
for j in range(nions):
    S1 += Ni[:,j]*Ai[j]*mp
    S2 += Ni[:,j]

meff = S1/S2
eps = xn/R0
anis = Te/(Ti+cor)-1
Omega = Za * q_e * B0/(Aa*mp)
Machi = (ma*Na+S1)/(S2*2*Ti + cor)*R0*Omega

for i in range(nx):
    dpsidx[i] = (2*xn[i]*qmag[i] - xn[i]**2*(2*(4-1)*xn[i]))/(qmag[i]**2+cor)

gradTi = np.gradient(Ti,xn)
gradNa = np.gradient(Na,xn)
for i in range(nions):
    gradNi[:,i] = np.gradient(Ni[:,i],xn)

AsymPhi[:] = 0
AsymN[:] = 1

# ------------------------------
# Part 2: Normalize + Fortran call
# ------------------------------

def arrF(a): return np.asfortranarray(np.asarray(a, dtype=np.float64))

raw_outputs = facitpy.facit_mod.facit(
    arrF(xn), arrF(theta), arrF(Za), float(Aa), arrF(Zi), arrF(Ai),
    arrF(Te), arrF(Ti), arrF(Ne), arrF(Ni), arrF(Na), arrF(Machi),
    arrF(gradTi), arrF(gradNi), arrF(gradNa),
    float(invaspct), float(B0), float(R0), arrF(qmag), arrF(dpsidx),
    arrF(FV), arrF(BV), arrF(RV), arrF(jacob),
    arrF(AsymPhi), arrF(AsymN),
    int(pol_asym), int(full_geom), int(solution), int(rotation),
    arrF(regulopt)
)

# ------------------------------
# Part 3: Output handling
# ------------------------------

class FACITResult:
    def __init__(self, out):
        (
            self.Flux_imp, self.Da, self.Da_M,
            self.Vconv, self.Vconv_M, self.nn,
            self.dmin, self.dmaj,
            self.Da_PS_M, self.Da_BP_M, self.Da_CL_M,
            self.Ka_PS_M, self.Ka_BP_M, self.Ka_CL_M,
            self.Ha_PS_M, self.Ha_BP_M, self.Ha_CL_M,
            self.Va_PS_M, self.Va_BP_M, self.Va_CL_M,
            self.Da_PS, self.Da_BP, self.Da_CL,
            self.Ka_PS, self.Ka_BP, self.Ka_CL,
            self.Ha_PS, self.Ha_BP, self.Ha_CL,
            self.Va_PS, self.Va_BP, self.Va_CL
        ) = out

        self.Da_comp = {"PS":self.Da_PS,"BP":self.Da_BP,"CL":self.Da_CL,
                        "PS_M":self.Da_PS_M,"BP_M":self.Da_BP_M,"CL_M":self.Da_CL_M}

result = FACITResult(raw_outputs)

# ------------------------------
# Save
# ------------------------------
with h5py.File("results_facit.h5","w") as f:
    f["Da"] = result.Da
    f["Flux_imp"] = result.Flux_imp
print("saved in results_facit.h5")

# ------------------------------
# Plots
# ------------------------------
plt.figure()
plt.plot(xn, result.Da)
plt.title("Da(x)")
plt.grid()
plt.show()

