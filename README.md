# CODENAME raman_LNSE (still in development)

raman_LNSE is an easy-to-use C++ implementation of the Split-Step Fourier Method for numeric calculation of a wave function's time-propagation under nonlinear Schrödinger equations.
It's code is based on the [ATUS-2 package](https://github.com/GPNUM/atus2) programm developed at ZARM (Center of Applied Space Technology and Microgravity, University of Bremen).  
[Read more in the TALISES documentation.](https://sascha.vowe.eu/talises-doc/)
### Features
- Calculation of a wavefunction's time propagation under a (non)linear Schrödinger equation: ![](https://latex.codecogs.com/gif.latex?i%5Chbar%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%20%5CPsi%20%3D%20%5CBig%5B%5Chat%7BV%7D%28%5CPsi%2C%5Cvec%7Br%7D%2Ct%29&plus;%5Cfrac%7B%5Chat%7Bp%7D%5E2%7D%7B2m%7D%5CBig%5D%5CPsi%28t%29)
- ![](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%5CPsi) includes arbitrary internal (e.g.: N-level system) and external (e.g.: 3D center-of-mass motion) degrees of freedom
- time and position dependent potentials ![](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20V%28%5Cvec%7Br%7D%2Ct%29)
- simple implementation of own Hamiltonians
- speed of C++