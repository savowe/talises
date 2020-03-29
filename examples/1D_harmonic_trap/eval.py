import talisestools as tt 
import numpy as np
from matplotlib import pyplot as plt

## Load binary data
# Load first state
psi_1_1 = tt.readbin("Seq_1_1.bin")
psi_2_1 = tt.readbin("Seq_2_1.bin")
psi_3_1 = tt.readbin("Seq_3_1.bin")
# Load second state
psi_1_2 = tt.readbin("Seq_1_2.bin")
psi_2_2 = tt.readbin("Seq_2_2.bin")
psi_3_2 = tt.readbin("Seq_3_2.bin")

# concatenate seperated wavefunction and time arrays
psi_1 = np.concatenate((psi_1_1["wavefunction"],psi_2_1["wavefunction"],psi_3_1["wavefunction"]), axis=1)
psi_2 = np.concatenate((psi_1_2["wavefunction"],psi_2_2["wavefunction"],psi_3_2["wavefunction"]), axis=1)
t = np.concatenate((psi_1_1["t"],psi_2_1["t"],psi_3_1["t"]))

# Get additional information for plots
xMin = psi_1_1["xMin"]
xMax = psi_1_1["xMax"]
nDimX = psi_1_1["nDimX"]

x = np.linspace(xMin, xMax, nDimX)
# Calculate probabillity densities
den_1 = np.abs(psi_1)**2
den_2 = np.abs(psi_2)**2
max_density = np.max(den_2)

## Create animated gif
import gif
@gif.frame
def plot(den1, den2, x, t, max_density):
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(x, den2, label = "t="+"{:10.3f}".format(t)+r' $\mu$s')
    ax1.set_ylabel(r"$|\Psi|^2$ [m$^{-1}$] excited")
    ax2.plot(x, den1)
    ax2.set_ylabel(r"$|\Psi|^2$ [m$^{-1}$] ground")
    ax2.set_xlabel(r"position [m]")
    ax1.set_xlim(np.min(x), np.max(x))
    ax1.set_ylim(0, max_density)
    ax2.set_ylim(0, max_density)
    ax1.grid()
    ax2.grid()
    ax1.legend()
    plt.tight_layout()

frames = []
for i in range(0,len(t)):
    frame = plot(den_1[:,i], den_2[:,i], x, t[i], max_density)
    frames.append(frame)
    print("Generated plot "+str(i)+"/"+str(len(t)-1))

gif.save(frames, "eval.gif", duration=200)