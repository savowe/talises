import talisestools as tt 
import numpy as np
from matplotlib import pyplot as plt

## Load binary data

data_psi1 = tt.readall(1)
data_psi2 = tt.readall(2)

# Get information for plots
psi1 = data_psi1["wavefunction"]
psi2 = data_psi2["wavefunction"]
t = data_psi1["t"]
xMin = data_psi1["xMin"]
xMax = data_psi1["xMax"]
nDimX = data_psi1["nDimX"]
x = np.linspace(xMin, xMax, nDimX)
max_density = np.max(np.abs(psi1)**2)

## Create animated gif
import gif
@gif.frame
def plot(psi1, psi2, x, t, max_density):
    den1_real = np.abs(np.real(psi1))**2
    den2_real = np.abs(np.real(psi2))**2
    den1_imag = np.abs(np.imag(psi1))**2
    den2_imag = np.abs(np.imag(psi2))**2
    den1 = den1_real + den1_imag
    den2 = den2_real + den2_imag
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.set_title(r"t = {:10.0f} $\mu$s".format(t))
    ax1.plot(x, den2, label = r"$|\Psi|^2$")
    ax1.plot(x, den2_real, label = r"|real$(\Psi )|^2$")
    ax1.plot(x, den2_imag, label = r"|imag$(\Psi )|^2$")
    ax1.set_ylabel(r"$|\Psi|^2$ [m$^{-1}$] excited")
    ax2.plot(x, den1, label = r"$|\Psi|^2$")
    ax2.plot(x, den1_real, label = r"|real$(\Psi )|^2$")
    ax2.plot(x, den1_imag, label = r"|imag$(\Psi )|^2$")
    ax2.set_ylabel(r"$|\Psi|^2$ [m$^{-1}$] ground")
    ax2.set_xlabel(r"position [m]")
    ax1.set_xlim(np.min(x), np.max(x))
    ax1.set_ylim(0, max_density)
    ax2.set_ylim(0, max_density)
    ax1.legend()
    ax2.legend()
    ax1.grid()
    ax2.grid()
    plt.tight_layout()

frames = []
for i in range(0,len(t)):
    frame = plot(psi1[:,i], psi2[:,i], x, t[i], max_density)
    frames.append(frame)
    print("Generated plot "+str(i)+"/"+str(len(t)-1))

gif.save(frames, "eval.gif", duration=100)