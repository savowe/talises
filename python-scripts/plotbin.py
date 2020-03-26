from readbin import readbin
import os,sys
import numpy as np
from matplotlib import pyplot as plt

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="Name of binary file")
parser.add_argument("--gif", help="create gif of files")
args = parser.parse_args()

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
        plt.tight_layout()
        plt.savefig(filename+".png")
        print(filename+".png has been created.")
        plt.close()

if __name__ == '__main__':

    dir_path = os.path.dirname(os.path.realpath(args.filename))


    if args.gif:
        import gif
        file_lists = []
        directory = args.filename

        ### File Search ###
        print("Searchin for files in "+directory)
        # define function for time sorting of file
        def filename_to_float(file):
            return float(file[0:-6])
        # Seach for files with certain ending and append to file_list
        for i in range(1, int(args.gif)+1):
            file_list = [f for f in os.listdir(directory) if f.endswith(str(i)+'.bin')]
            file_list.sort(key = filename_to_float)
            file_lists.append(file_list)

        # Define plot for animation
        @gif.frame
        def plot(filelist, i_t):
            N = len(file_lists)
            fig, ax = plt.subplots(N)
            #plt.figure(figsize=(5, 3), dpi=100)
            data = readbin(directory+file_lists[0][0])
            max_density = 0
            for i in range(0,N):
                data = readbin(directory+file_lists[i][i_t])
                x = np.arange(data["xMin"], data["xMax"], data["dx"])
                psi2 = np.abs(np.power(data["wavefunction"],2))
                ax[i].plot(x, psi2, label="t="+str(data["t"]))
            plt.ylim(0,max_density)
            plt.tight_layout()
            plt.legend()

        frames = []
        for j in range(0, len(file_lists[0])):
            frame = plot(file_lists, j)
            frames.append(frame)
        gif.save(frames, directory+"talises.gif", duration=200)

        sys.exit()
    else:
        plotbin(sys.argv[1])