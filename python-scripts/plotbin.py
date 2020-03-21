from readbin import readbin
from sys import argv
from numpy import abs, arange, power
from matplotlib import pyplot as plt

def plotbin(filename):
    data = readbin(filename)
    if len(data["wavefunction"].shape)==1:
        x = arange(data["xMin"], data["xMax"], data["dx"])
        psi2 = abs(power(data["wavefunction"],2))
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
        x = arange(data["xMin"], data["xMax"], data["dx"])
        y = arange(data["yMin"], data["yMax"], data["dy"])
        psi2= abs(power(data["wavefunction"],2))
        plt.title(r"$|\Psi |^2$ at t="+str(data["t"]))
        plt.ylabel("y")
        plt.xlabel("x")
        plt.pcolormesh(x, y, psi2)
        plt.grid()
        plt.tight_layout()
        plt.savefig(filename+".png")
        print(filename+".png has been created.")
        plt.close()

if __name__ == '__main__':
    plotbin(argv[1])

