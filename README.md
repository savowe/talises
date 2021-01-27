# TALISES

**TALISES** (This Ain't a LInear Schrödinger Equation Solver) is an easy-to-use C++ implementation of the Split-Step Fourier Method, for numeric calculation of a wave function's time-propagation under the Schrödinger equation.  

### Features
- Calculation of a wavefunction's time propagation under a (non)linear Schrödinger equation: ![](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20i%5Chbar%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%20%5CPsi%20%28%5Cvec%7Br%7D%2C%20t%29%20%3D%20%5CBig%5BV%28%5CPsi%2C%5Cvec%7Br%7D%2C%20t%29%20&plus;%20%5Cfrac%7B%5Chbar%5E2%7D%7B2m%7D%5Cnabla%5E2%20%5CBig%5D%20%5CPsi%20%28%5Cvec%7Br%7D%2C%20t%29)
- the wave-function ![](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%5CPsi) may include an arbitrary number of internal and external degrees of freedom
- simple implementation of Hamiltonians
- speed of C++, the FFTW and GSL libaries and multithreading

[Find more information and exemplary simulations in the documentation.](https://savowe.github.io/talises-doc/)

# Installing TALISES
TALISES source files can be obtained on [GitHub](https://github.com/savowe/talises) either by downloading it there or by directly cloning it with
```
git clone https://github.com/savowe/talises.git
```

TALISES uses the following libaries:

- [FFTW](http://www.fftw.org/) for Discrete-Fourier-Transforms
- [GSL(CBLAS)](https://www.gnu.org/software/gsl/) for numeric diagonalization and matrix exponentiation
- [muparser](https://beltoforion.de/article.php?a=muparser) for parsing of mathematical formulas

You need this three libaries in order to use TALISES. If they are already installed on your machine, you only need the libary paths in your environment variables (e.g LD_LIBRARY_PATH) whilst you compile the program by running `cmake .`, `make clean` and `make`.  
TALISES also comes with an installation script which is written in Common Lisp, that will download the required libaries, compile them and create [environment-module](http://modules.sourceforge.net/) files. The use of environment-modules is advised. 

### Prerequisites
Necessities are compilers (gcc, g++), build automation tools (make, cmake) and the [Boost C++ Libraries](https://www.boost.org/doc/libs/).
You can install these with
```text
sudo apt install build-essential cmake libboost-all-dev
```
If you want to use the install script that comes with TALISES you need a Common Lisp compiler (we will use [SBCL](http://www.sbcl.org/)) and environment-modules. They can be installed with
```text
sudo apt install sbcl curl
```
On newer versions of Ubuntu you can install environment-modules via
```text
sudo apt install tclsh environment-modules
```
After installation of environment-modules you may need to execute `add.modules` and reboot the system.
If this does not work we recommend installing [it from source.](http://modules.sourceforge.net/)

### Running the installation script

If you downloaded the git-repository and are within the directory, you can start the installation with
```text
sbcl --script install.lisp
```
The script will check for some dependencies e.g. compilers.
If all tests pass you will be greeted by the installer.
Here you have severall options to chose from.

````text
  _________    __    _________ ___________
 /_  __/   |  / /   /  _/ ___// ____/ ___/
  / / / /| | / /    / / \__ \/ __/  \__ \ 
 / / / ___ |/ /____/ / ___/ / /___ ___/ / 
/_/ /_/  |_/_____/___//____/_____//____/  

Options:
-h, --help            Print this help text
--no-fetch            Skip Download
--no-modules          Skip Module file installation
--no-install          Skip Installation
-j THREADS            Number of make threads
--build-dir DIR       Build directory
--install-dir DIR     Installation directory
--module-dir DIR      Module directory

What do you want to install?
Press number of each package to be installed and then press ENTER:
0 - gsl              (2.5)
1 - muparser         (2.2.6.1)
2 - fftw             (3.3.8)
3 - talises          (git)
a - all
q - Abort Installation.
````
The build directory will contain the source files, the installation directory the compiled binaries and the module directory the module files which will allow you to quickly switch between different environments.  
We recommend to install the libaries individually (first gsl, then muparser, etc.), so if something goes wrong you can easily identify the problems.  
When GSL, muparser and FFTW are installed you need to make sure that your environment variables point to the libary directories. Otherwise TALISES can not be compiled. Type `module avail` to see a list of the environment files you have. If environment-modules can find the folder you saved the modulefiles in you should see
````text
----- /home/username/local/modules/modulefiles/ ------
fftw-3.3.8       gsl-2.5          muparser-2.2.6.1
````
Then you can load those configuration files with `module load fftw-3.3.8 muparser-2.2.6.1 gsl-2.5` and check with `echo $LD_LIBRARY_PATH` whether the libaries' directories are in it.  
Now you can install TALISES either via the installation script or by running `cmake .` followed by `make clean` and `make`.  
If everything went right you will find the compiled binaries in the installation directory you set.

[Find more information and exemplary simulations in the documentation.](https://savowe.github.io/talises-doc/)
