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

# radial and poloidal grids
xn = np.linspace(0, 0.99, nx)
theta = np.linspace(0, 2*np.pi, nth, endpoint=True)

# allocate arrays
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

term_iCRH = np.zeros(nx)
term_rot = np.zeros(nx)
delta_phi  = np.zeros(nx)
Deltam_phi  = np.zeros(nx)
S1 = np.zeros(nx)

pol_asym = True
rotation = True
solution = 3
full_geom = True

ni0 = [4e19,1e19,1e14,5e9,2e10]

# ------------------------------
# Densities
# ------------------------------
Ni[:,0] = ni0[0]*(1-xn**2)**2
Ni[:,1] = ni0[1]*(1-xn**2)**2
Ni[:,2] = ni0[2]*(1-xn**2)**2
Ni[:,3] = ni0[3]*(0.3+np.exp(6*(xn-1)))
Ni[:,4] = ni0[4]*(0.3+np.exp(6*(xn-1)))

Na = 6e12*(0.3+np.exp(4*(xn-1)))
Ne = 6e19*(1-xn**2)**2

Zi[:,0] = 1
Zi[:,1] = 1

for i in range(nx):
    Zi[i,2] = 1 if xn[i]>=0.6 else 2
    FV[i] = R0*B0

Zi[:,3] = np.clip(6 - 5*np.exp(-((xn-1)/0.25)**2), 1.0, 20.0)
Zi[:,4] = np.clip(7 - 6*np.exp(-((xn-1)/0.25)**2), 1.0, 20.0)
Za = np.clip(60 - 55*np.exp(-((xn-1)/0.25)**2), 1.0, 20.0)

qmag = 1 + (qa-1)*xn**2

# ------------------------------
# Temperaturas
# ------------------------------
T_floor = 0.05
Ti = 1.5*(1-xn**2)**2 + T_floor
Te = 3.0*(1-xn**2)**2 + 2.0*T_floor

# ------------------------------
# Geometry: RV, BV, jacob
# ------------------------------
min_denom = 1e-3
min_RV = 1e-2
min_qmag = 1e-6

for j in range(nth):
    denom = np.clip(1.0 - invaspct*np.cos(theta[j]), min_denom, None)

    RV[:, j] = R0*(1.0 + xn*invaspct*np.cos(theta[j]))
    RV[:, j] = np.clip(RV[:, j], min_RV, None)

    RV2 = RV[:, j]**2
    safe_RV2 = np.maximum(RV2, min_RV**2)
    safe_qmag2 = np.maximum(qmag**2, min_qmag**2)

    BV[:, j] = FV**2/safe_RV2 * (1.0 + invaspct**2/safe_RV2 * xn**2/safe_qmag2)

    jacob[:, j] = xn*R0/denom

# ------------------------------
# S1, S2, meff, anis, Omega, Machi
# ------------------------------
S1 = 0
S2 = Na.copy()
for j in range(nions):
    S1 += Ni[:,j]*Ai[j]*mp
    S2 += Ni[:,j]

meff = S1/np.maximum(S2,1e-12)
eps = xn*invaspct
anis = Te/Ti - 1.0
Omega = Za*q_e*B0/(Aa*mp)

with np.errstate(divide='ignore', invalid='ignore'):
    Machi_raw = (ma*Na + S1)*R0**2 * Omega/(S2*2.0*Ti)

Machi_raw = np.nan_to_num(Machi_raw)
Machi = np.clip(Machi_raw, 0.0, 0.8)

# ------------------------------
# dpsidx
# ------------------------------
for i in range(nx):
    dpsidx[i] = (2*xn[i]*qmag[i] - xn[i]**2 * (2*(4-1)*xn[i]))/max(qmag[i]**2,1e-12)
dpsidx[0] = dpsidx[1]

# ------------------------------
# Gradients
# ------------------------------
for i in range(nx):
    if i == 0:
        dx = xn[i+1] - xn[i]
        gradTi[i] = (Ti[i+1] - Ti[i]) / dx
        gradNa[i] = (Na[i+1] - Na[i]) / dx
        for j in range(nions):
            gradNi[i,j] = (Ni[i+1,j] - Ni[i,j]) / dx

    elif i == nx-1:
        dx = xn[i] - xn[i-1]
        gradTi[i] = (Ti[i] - Ti[i-1]) / dx
        gradNa[i] = (Na[i] - Na[i-1]) / dx
        for j in range(nions):
            gradNi[i,j] = (Ni[i,j] - Ni[i-1,j]) / dx

    else:
        dx = xn[i+1] - xn[i-1]
        gradTi[i] = (Ti[i+1] - Ti[i-1]) / dx
        gradNa[i] = (Na[i+1] - Na[i-1]) / dx
        for j in range(nions):
            gradNi[i,j] = (Ni[i+1,j] - Ni[i-1,j]) / dx

# ------------------------------
# Asymmetry
# ------------------------------
for i in range(nx):
    ratio_TeTi = max(Te[i]/Ti[i],1e-8)
    term_iCRH[i] = anis[i]/ratio_TeTi
    term_rot[i] = 2.0*np.sqrt(max(Machi[i],0.0))

    denom_phi = max(1.0 + Za[i]*Te[i]/Ti[i],1e-6)
    delta_phi[i] = eps[i]/denom_phi*(term_iCRH[i]+term_rot[i])
    Deltam_phi[i] = -0.165*eps[i]

AsymPhi[:,0] = delta_phi
AsymPhi[:,1] = Deltam_phi

for j in range(nions):
    AsymN[:,j,0] = -Zi[:,j]*(Te/Ti)*delta_phi
    AsymN[:,j,1] = -Zi[:,j]*(Te/Ti)*Deltam_phi

# ------------------------------
# Part 2: FACIT call
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
# Output structured
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

result = FACITResult(raw_outputs)

# ------------------------------
# Save EVERYTHING in HDF5
# ------------------------------
with h5py.File("resultado_facit_completo.h5","w") as f:

    # Inputs
    inputs = f.create_group("inputs")
    inputs["xn"] = xn
    inputs["theta"] = theta
    inputs["Za"] = Za
    inputs["Zi"] = Zi
    inputs["Ai"] = Ai
    inputs["Te"] = Te
    inputs["Ti"] = Ti
    inputs["Ne"] = Ne
    inputs["Ni"] = Ni
    inputs["Na"] = Na
    inputs["Machi"] = Machi
    inputs["gradTi"] = gradTi
    inputs["gradNi"] = gradNi
    inputs["gradNa"] = gradNa
    inputs["qmag"] = qmag
    inputs["dpsidx"] = dpsidx
    inputs["FV"] = FV
    inputs["BV"] = BV
    inputs["RV"] = RV
    inputs["jacob"] = jacob
    inputs["AsymPhi"] = AsymPhi
    inputs["AsymN"] = AsymN

    # Params
    params = f.create_group("params")
    params["Aa"] = Aa
    params["R0"] = R0
    params["B0"] = B0
    params["qa"] = qa
    params["invaspct"] = invaspct

    # Outputs
    outputs = f.create_group("outputs")
    for i, arr in enumerate(raw_outputs):
        outputs[f"out_{i}"] = np.asarray(arr)

    outputs["Da"] = result.Da
    outputs["Flux_imp"] = result.Flux_imp

print("✔ Todo guardado en resultado_facit_completo.h5 🚀")

# ------------------------------
# Plot
# ------------------------------
plt.figure()
plt.plot(xn, result.Da)
plt.title("Da(x)")
plt.grid()
plt.show()

