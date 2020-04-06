import talisestools as tt 
import numpy as np
from matplotlib import pyplot as plt

## Load binary data

data_psi1 = tt.readbin("Seq_1_1.bin")
data_psi2 = tt.readbin("Seq_1_2.bin")
data_psi3 = tt.readbin("Seq_1_3.bin")

# Get information for plots
psi1 = data_psi1["wavefunction"]
psi2 = data_psi2["wavefunction"]
psi3 = data_psi3["wavefunction"]
t = data_psi1["t"]

# Integration over x for population number
N1 = np.sum(np.power(np.abs(psi1),2), axis=0)*data_psi1["dx"]
N2 = np.sum(np.power(np.abs(psi2),2), axis=0)*data_psi2["dx"]
N3 = np.sum(np.power(np.abs(psi3),2), axis=0)*data_psi3["dx"]

fig, ax = plt.subplots(dpi=150)
ax.plot(t, N1, label = r"$|g\rangle$")
ax.plot(t, N2, label = r"$|e\rangle$")
ax.plot(t, N3, label = r"$|i\rangle$")
ax.set_ylabel("N")
ax.set_xlabel(r"t [$\mu$s]")
plt.legend()
plt.grid()
fig.savefig("eval.png")