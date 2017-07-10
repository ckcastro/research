README

##Synopsis

This repository contains all the source and output files required to simulate the propagation dynamics of light in a discrete array of N optical fibers in Kerr media described by the one-dimensional disordered discrete Nonlinear Schroedinger Equation. The solver is based in a ABC4 symplectic integrator method -- O(h^4) accurate. 

This repository contains:

Source files:

   Makefile          driver_ddnls_nonlinear_pert_ic.cpp
   abc4.cpp          matrix.cpp
   abc4.hpp          matrix.hpp
   ddnls_submit.sh   plot_appearance.m
   rhs.hpp           README
 
Output files
   
   NormDensitySpaceShortTime.pdf 	
   NormDensitySpaceTime.pdf 	
   energyrelerror_abc4.pdf 	
   finalNormDensity_abc4.pdf 	
   normrelerror_abc4.pdf 	

##Usage 

How To Run "driver.exe":
	Open Terminal. Change the current working directory to your local project. Type "Makefile". 
        Type "./driver.exe"


##Contributors

Claudia Castro-Castro

Daniel Reynolds


##Reference: 
     Skokos, Ch, E. Gerlach, J. D. Bodyfelt, G. Papamikos, and S. Eggl. 
     "High order three part split symplectic integrators: Efficient techniques for the long  
     time simulation of the disordered discrete nonlinear Schrödinger equation." 
     Physics Letters A 378, no. 26 (2014): 1809-1815.

