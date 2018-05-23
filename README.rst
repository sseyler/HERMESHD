======================
HERMESHD
======================

HERMESHD is a Discontinuous Galerkin (DG) numerical code written in Fortran and Python for simulating 3D compressible hydrodynamics.

:Authors:      Sean L. Seyler (FH and modularized code), Charles Seyler (original code),
:Organization: Arizona State University, Cornell University
:Contact:      slseyler@asu.edu, ces7@cornell.edu
:Year:         2017
:License:      MIT License
:Copyright:    © 2017 Sean L. Seyler

----------------
Using this work
----------------

This work is currently distributed under the terms of the MIT license, though this may change in the near feature. If you choose to use any part of this work or these ideas, I kindly request that you cite `my talk presented at APS March Meeting 2018`_ [3]_—at least until we make our manuscript (in preparation) available! Many thanks for your consideration!

**Important note**: The original hydrodynamics code, derived from the PERSEUS XMHD code, was written by Charles E. Seyler and Matthew R. Martin [1]_; HERMESHD is based on the PERSEUS hydrodynamics derivative and will be adapted to fluctuating hydrodynamics (FH) simulations by Sean Seyler. See Background below!

-----------
Background
-----------

The original PERSEUS (Physics of the Extended-mhd Relaxation System using an Efficient Upwind Scheme) numerical code is a 3D finite volume (FV) method developed by Matt Martin and Charles Seyler in 2011 for solving the extended magnetohydrodynamics (XMHD) equations [1]_. PERSEUS finds applications in simulating High Energy Density (HED) plasmas, such as dense z-pinches, where a wide range of dynamical length scales and densities are encountered. In 2014, Xuan Zhao, Nat Hamlin, and Charles Seyler developed a discontinuous Galerkin method extending the original PERSEUS FV algorithm, further improving its accuracy and computational efficiency in HED applications [2]_. HERMESHD (Hyperbolic Equations and Relaxation Model for Extended Systems of HydroDynamics) is a specialized extension of the PERSEUS XMHD DG algorithm that solves the compressible Euler, Navier-Stokes, and higher-moment equations in three dimensions. The code is actively being extended by Sean Seyler (as part of his Blue Waters Graduate Fellowship Project) to (1) model fluctuating hydrodynamics, and (2) form a basis for hybrid atomistic-continuum simulation (with a view toward biological macromolecular simulation).


---------
Overview
---------

* Variables (5, 10, or 13 moment-equations):

  * 5 independent field variables from conventional hydrodynamics

    * density (``rh``)
    * velocity (``vx``, ``vy``, ``vz``)
    * energy (``en``)

  * 5 variables describing viscous stresses and 3 for heat flux

    * stress (``pxx``, ``pyy``, ``pzz``, ``pxy``, ``pxz``, ``pyz``) – 5 independent for symmetric, traceless stress
    * heat flux (``qx``, ``qy``, ``qz``)

* Units—a value of unity for each variable or parameter corresponds to the following dimensional units

  * length (``L0``)
  * time (``t0``)
  * number density (``n0``)
  * velocity (``v0``)
  * temperature (``te0``)



----------------
Getting Started
----------------

**Note**: this section (and this README.rst) are a work-in-progress and will be continuously updated and refined. Feedback and contributions are welcome!

Prerequisites
==============

1. `CMake utility`_ (version 3.6.2 or higher)

1. Fortran compiler:
    *  Intel Fortran compiler (``ifort``)
    *  GNU Fortran (``gfortran``)

1. MPI library: Any should work, but OpenMPI is a good starting point.

1. To use the Python front end (optional)
    *  Python 2.7.11 or higher
    *  ``f2py`` (part of ``numpy``)
    *  ``f90wrap`` (see the `f90wrap GitHub page`_)


From GitHub to simulation
==========================

1. Clone the repository
    ``git clone git@github.com:sseyler/HERMESHD.git``

1. Change to toplevel HERMESHD directory (should contain ``CMakeLists.txt``)
    ``cd /path/to/HERMESHD``

1. Generate a ``Makefile`` using ``cmake``
    .. code-block:: console

        mkdir build && cd build
        cmake ..
        make
        mpirun -n 6 hermes


-----------------
Acknowledgements
-----------------

This work was funded by a `2016 Blue Waters Graduate Fellowship`_ for `Developing a hybrid continuum-particle method for simulating large-scale heterogeneous biomolecular systems`_. Sean Seyler also expresses his eternal gratitude to Professor Oliver Beckstein who, as his doctoral advisor, was very generous in allowing him to pursue this project in his final year of graduate school.

-----------
References
-----------

.. Articles
.. --------

.. [1] C.E. Seyler & M.R. Martin.
   Relaxation model for extended magnetohydrodynamics: Comparison
   to magnetohydrodynamics for dense Z-pinches. *Phys. Plasmas* **18**,
   012703 (2011). doi:`10.1063/1.3543799`_.

.. _`10.1063/1.3543799`: http://dx.doi.org/10.1063/1.3543799

.. [2] X. Zhao, Y. Yang & C.E. Seyler.
   A positivity-preserving semi-implicit discontinuous Galerkin scheme
   for solving extended magnetohydrodynamics equations. *J. Comput. Phys.*
   **278**, 400–415 (2014). doi:`10.1016/j.jcp.2014.08.044`_.

.. _`10.1016/j.jcp.2014.08.044`: http://dx.doi.org/10.1016/j.jcp.2014.08.044

.. [3] S.L. Seyler, C.E. Seyler & O. Beckstein.
    *Fluctuating Hydrodynamics in the 13-moment Approximation for
    Simulating Biomacromolecular Nanomachines*. Talk, APS March Meeting 2018.
    url:`meetings.aps.org/Meeting/MAR18/Session/S51.5`_.

.. _`meetings.aps.org/Meeting/MAR18/Session/S51.5`: https://meetings.aps.org/Meeting/MAR18/Session/S51.5

.. _`2016 Blue Waters Graduate Fellowship`: https://bluewaters.ncsa.illinois.edu/fellowships/2016

.. _`Developing a hybrid continuum-particle method for simulating large-scale heterogeneous biomolecular systems`: https://bluewaters.ncsa.illinois.edu/science-teams?page=detail&psn=bafh

.. _`my talk presented at APS March Meeting 2018`: https://meetings.aps.org/Meeting/MAR18/Session/S51.5

.. _`CMake utility`: https://cmake.org/

.. _`f90wrap GitHub page`: https://github.com/jameskermode/f90wrap
