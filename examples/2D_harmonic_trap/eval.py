import talisestools as tt 
import os,sys
import numpy as np

#input
directory = os.path.abspath(os.getcwd())
n_int_state = 1
### File Search ###
print("Searchin for files in "+directory)
# define function for time sorting of file
def filename_to_float(file):
    if file[0:4] == "Seq_":
        return float file[0:3]
    return float(file[0:-6])

# Read and append data
file_list = []
file_list = [f for f in os.listdir(directory) if f.endswith(str(n_int_state)+'.bin')]
data = []
for i in range(0, len(file_list)):
    data.append(tt.readbin(file_list[i]))

# concatenate psi and time
psi = data[0]["wavefunction"]
t = data[0]["t"]
for i in range(1, len(file_list)):
    t = np.concatenate((t, data[i]["t"]))
    psi = np.concatenate((psi, data[1]["wavefunction"]), axis=data[0]["nDims"])

# sort after time
sorted_inds = t.argsort()
sorted_t = np.empty(shape=t.shape)
sorted_psi = np.empty(shape=psi.shape,dtype=np.complex_)
if data[0]["nDims"] == 1:
    for i in range(0, len(t)):
        sorted_psi[:,i] = psi[:,sorted_inds[i]]
        sorted_t[i] = t[sorted_inds[i]]

if data[0]["nDims"] == 2:
    for i in range(0, len(t)):
        sorted_psi[:,:,i] = psi[:,:,sorted_inds[i]]
        sorted_t[i] = t[sorted_inds[i]]

if data[0]["nDims"] == 3:
    for i in range(0, len(t)):
        sorted_psi[:,:,:,i] = psi[:,:,:,sorted_inds[i]]
        sorted_t[i] = t[sorted_inds[i]]
