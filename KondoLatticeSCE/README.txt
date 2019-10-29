  ------------------------------------------------------
     KondoLatticeSCE - Phase diagram calculation code
  ------------------------------------------------------

Author: eightbitastronomy (pseudonymously)
Contact: eightbitastronomy@protonmail.com

  This code represents a portion of the code used for my doctoral
thesis on the Kondo lattice model in 3 dimensions. Because I
do not personally regard it as "good" code, I have provided it for
two reasons.
  The first is to share the sm library. The library files, sm.h
and sm.c, contain functions for forking processes to make
calculations and to return the calculation results via shared
memory to the parent process. While there exist some powerful
tools for running parallelized computations on compute nodes, the
sm library might be more useful/appropriate for a person
attempting to make full use of modern, multi-core desktop.
  The second is for any scholarly interest/concern from a
physicist who might wish to see the code used to produce the
data in my thesis. Such a person might also wish to see
the code which I have not provided here, such as that used for
band calculations. Please contact me if this is the case.
  As for the code provided, it is enough to produce a diagram for
normal paramagnetic, Kondo paramagnetic, and magnetically-ordered
phases.

Copyright information:
  KondoLatticeSCE code is copyrighted by eightbitastronomy, 2018
  Mersenne twister code is copyrighted by its authors. Please see
    the dSFMT folder and its contents for relevant copyright
    information.

License information:
  KondoLatticeSCE code is provided under the GPL v. 3. This
    license and its terms are found in the accompanying file,
    COPYING.txt.
  Mersenne twister code is provided under a BSD license. For
    more information, please see the dSFMT folder and its
    contents.

28 Oct 2019
