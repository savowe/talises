# CODENAME raman_LNSE (still in development)

raman_LNSE is an easy-to-use C++ implementation of the Split-Step Fourier Method for numeric calculation of a wave function's time-propagation under nonlinear Schrödinger equations.
It's code is based on the [ATUS-2 package](https://github.com/GPNUM/atus2) programm developed at ZARM (Center of Applied Space Technology and Microgravity, University of Bremen).
### Features
- Calculation of a wavefunction's time propagation under a (non)linear Schrödinger equation: ![](https://latex.codecogs.com/gif.latex?i%5Chbar%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%20%5CPsi%20%3D%20%5CBig%5B%5Chat%7BV%7D%28%5CPsi%2C%5Cvec%7Br%7D%2Ct%29&plus;%5Cfrac%7B%5Chat%7Bp%7D%5E2%7D%7B2m%7D%5CBig%5D%5CPsi%28t%29)
- ![](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%5CPsi) includes arbitrary internal (e.g.: N-level system) and external (e.g.: 3D center-of-mass motion) degrees of freedom
- time and position dependent potentials ![](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20V%28%5Cvec%7Br%7D%2Ct%29)
- simple implementation of own Hamiltonians
- speed of C++
-------------
## Example: Light-Pulse Atom-Interferometry using stimulated Raman transitions
Light-pulse atom interferometry is a versatile tool for modern metrology, and can be used to measure (gravitational) acceleration, rotation, the fine-structure and gravitaitonal constant, gravitational waves, test Lorentz invariance, general relativity, and dark energy theories.  
In this example we combine important aspects of the raman_LNSE program and show how to simulate realistic light-pulse atom interferometry experiments including the description of internal as well as external degrees of freedom of a wave function subjected to a Hamitlonian of the form  
![](https://latex.codecogs.com/gif.latex?%5Chat%7BH%7D%20%3D%20%5Chat%7BV%7D%28%5Cvec%7Br%7D%2Ct%29&plus;%5Cfrac%7B%5Chat%7Bp%7D%5E2%7D%7B2m%7D.)
#### Seperabillity of internal and external dynamics
A composite system of interacting particles (like an atom) can be conveniently described using center of mass *R* and relative position *r* coordinates ([see Shore 2011, ch. 17.2](#sources)).
The wave function can thus be written as  
![](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cpsi%28%5Cvec%7BR%7D%2C%5Cvec%7Br%7D%29%20%3D%20%5Cpsi_%7BC%7D%28%5Cvec%7BR%7D%29%20%5Cpsi_%7Bnlm%7D%28%5Cvec%7Br%7D%29.)  
Here *C* denotes the wave function describing a center of mass motion and the indices *n,l,m* the part of the wave function describing the interaction of nuclei and electrons which give rise to some atomic orbital structure.  
The interaction of atom and light will be modelled using the dipole approximation (eletric dipole transitions are the dominant transitions between fine-structure states). The interaction of a charged particle with an eletric field is described with the transformation ![](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%5Cvec%7Bp%7D%20%5Crightarrow%20%5Cvec%7Bp%7D%20-%20%5Cfrac%7Be%7D%7Bc%7D%20%5Cvec%7BA%7D%28%5Cvec%7Br%7D%2Ct%29). In the dipole approximation we take the first order term of the vector potential's Taylor expansion at the center of mass coordinate *R*.
Thus, to first order, we infer that the laser's eletric field is spatially constant across the orbital wave function. This results in the famous eletric dipole term: ![](https://latex.codecogs.com/gif.latex?V%28%5Cvec%7Br%7D%2Ct%29%20%5Capprox%20-d%5Ccdot%20E%28%5Cvec%7BR%7D%2Ct%29)

### TBC
The set of basis vectors describing the internal strucure of the system is denoted as $|n\rangle $. The level structure of an atom is usualy very complicated and includes potentially infinetely many levels, however, for many physical situtations it is sufficent to approximate it as a N-level system. In the case of Raman transitions, we can approximate the atom as a three level system with basis vectors ![](https://latex.codecogs.com/gif.latex?%7Cg%5Crangle),![](https://latex.codecogs.com/gif.latex?%7Ce%5Crangle) and ![](https://latex.codecogs.com/gif.latex?%7Ci%5Crangle).  
Our
This means that our Hilbert space will be will be constituent of 



----------------------
### Sources
Shore, B. (2011). Manipulating Quantum Structures Using Laser Pulses. Cambridge: Cambridge University Press
