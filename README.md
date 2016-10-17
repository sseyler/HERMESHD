# PERSEUS Hydrodynamics

PERSEUSHydro is a Discontinuous Galerkin (DG) hydrodynamics code that solves the 3D compressible Euler equations.

The original PERSEUS algorithm (Physics of the Extended-mhd Relaxation System using an Efficient Upwind Scheme) is a finite volume method for solving the extended magnetohydrodynamics (XMHD) equations (by Matt Martin and Charles Seyler), which was further extended to the discontinuous Galerkin method (by Xuan Zhao, Nat Hamlin, and Charles Seyler). PERSEUSHydro is thus an extension of the PERSEUS 3D DG XMHD algorithm specialized to 3D compressible hydrodynamics.

There are 5 independent variables: density (rh), velocity (vx, vy, vz), and energy (en).
