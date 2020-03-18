from readbin import readbin
import os,sys
import numpy as np
from matplotlib import pyplot as plt

def plotbin(filename):
    data = readbin(filename)
    if len(data["wavefunction"].shape)==1:
        x = np.arange(data["xMin"], data["xMax"], data["dx"])
        psi2 = np.abs(np.power(data["wavefunction"],2))
        plt.ylabel(r"$|\Psi |^2$")
        plt.xlabel("position")
        plt.plot(x, psi2, label="t="+str(data["t"]))
        plt.xlim(data["xMin"], data["xMax"])
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(filename+".png")
        print(filename+".png has been created.")
        plt.close()

    if len(data["wavefunction"].shape)==2:
        x = np.arange(data["xMin"], data["xMax"], data["dx"])
        y = np.arange(data["yMin"], data["yMax"], data["dy"])
        psi2= np.abs(np.power(data["wavefunction"],2))
        plt.title(r"$|\Psi |^2$ at t="+str(data["t"]))
        plt.ylabel("y")
        plt.xlabel("x")
        plt.pcolormesh(x, y, psi2)
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(filename+".png")
        print(filename+".png has been created.")
        plt.close()

if __name__ == '__main__':
    plotbin(sys.argv[1])

