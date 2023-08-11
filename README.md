Update 08/11/23: Added version to support Quantum Espresso, not only VASP.

This is a simple script to plot the dynamical structure factor using Phonopy force_sets/force_constants

The calculation of dynamical structure factor is as same as Phonopy website:
http://phonopy.github.io/phonopy/dynamic-structure-factor.html

If BORN file is within the same folder, it seems it can include LO-TO splitting automatically. But I need to double-check

SQE folder is for SQE simulation
Thermal_diffuse folder utilizes SQE script to integrate over energy to calculate thermal diffuse patterns

It seems there is some band folding issue after phonopy version 2.16.0. Please use 2.15.1 or 2.15.0 for this script.
