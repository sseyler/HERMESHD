======================
PERSEUS Hydrodynamics
======================

PERSEUSHydro is a Discontinuous Galerkin (DG) numerical code written in Fortran and Python for simulating 3D compressible hydrodynamics.

:Authors:      Charles Seyler (primary), Sean Seyler
:Organization: Cornell University, Arizona State University
:Contact:      ces7@cornell.edu, slseyler@asu.edu
:Year:         2016
:License:      GNU Public License, version 3 (or higher)
:Copyright:    © 2016 Charles Seyler

Background
===========

The original PERSEUS (Physics of the Extended-mhd Relaxation System using an Efficient Upwind Scheme) numerical code is a 3D finite volume (FV) method developed by Matt Martin and Charles Seyler in 2011 for solving the extended magnetohydrodynamics (XMHD) equations [1]_. PERSEUS finds applications in simulating High Energy Density (HED) plasmas, such as dense z-pinches, where a wide range of dynamical length scales and densities are encountered. In 2014, Xuan Zhao, Nat Hamlin, and Charles Seyler developed a discontinuous Galerkin method extending the original PERSEUS FV algorithm, further improving its accuracy and computational efficiency in HED applications [2]_. PERSEUSHydro is a specialized extension of the PERSEUS XMHD DG algorithm that solves the compressible Euler equations in three dimensions.

Overview
=========

* Variables—there are 5 independent field variables

  * density (``rh``)
  * velocity (``vx``, ``vy``, ``vz``)
  * energy (``en``)

* Units—a value of unity for each variable or parameter corresponds to the following dimensional units

  * length (``L0``)
  * time (``t0``)
  * number density (``n0``)
  * velocity (``v0``)
  * temperature (``te0``)

Extension to Fluctuating Hydrodynamics
=======================================

More to come soon!

References
===========

.. Articles
.. --------

.. [1] Seyler, C. E. & Martin, M. R.
   Relaxation model for extended magnetohydrodynamics: Comparison
   to magnetohydrodynamics for dense Z-pinches. Phys. Plasmas 18,
   012703 (2011). doi:`10.1063/1.3543799`_.

.. _`10.1063/1.3543799`: http://dx.doi.org/10.1063/1.3543799

.. [2] Zhao, X., Yang, Y. & Seyler, C. E.
   A positivity-preserving semi-implicit discontinuous Galerkin scheme
   for solving extended magnetohydrodynamics equations. J. Comput. Phys.
   278, 400–415 (2014). doi:`10.1016/j.jcp.2014.08.044`_.

.. _`10.1016/j.jcp.2014.08.044`: http://dx.doi.org/10.1016/j.jcp.2014.08.044
