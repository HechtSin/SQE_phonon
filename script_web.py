#!/usr/bin/env python

import numpy as np
from phonopy import load
from phonopy.spectrum.dynamic_structure_factor import atomic_form_factor_WK1995
from phonopy.phonon.degeneracy import degenerate_sets
from phonopy.units import THzToEv


def get_AFF_func(f_params):
    def func(symbol, Q):
        return atomic_form_factor_WK1995(Q, f_params[symbol])
    return func


def run(phonon,
        Qpoints,
        temperature,
        atomic_form_factor_func=None,
        scattering_lengths=None):
    # Transformation to the Q-points in reciprocal primitive basis vectors
    Q_prim = np.dot(Qpoints, phonon.primitive_matrix)
    # Q_prim must be passed to the phonopy dynamical structure factor code.
    phonon.run_dynamic_structure_factor(
        Q_prim,
        temperature,
        atomic_form_factor_func=atomic_form_factor_func,
        scattering_lengths=scattering_lengths,
        freq_min=1e-3)
    dsf = phonon.dynamic_structure_factor
    q_cartesian = np.dot(dsf.qpoints,
                         np.linalg.inv(phonon.primitive.get_cell()).T)
    distances = np.sqrt((q_cartesian ** 2).sum(axis=1))

    print("# [1] Distance from Gamma point,")
    print("# [2-4] Q-points in cubic reciprocal space, ")
    print("# [5-8] 4 band frequencies in meV (becaues of degeneracy), ")
    print("# [9-12] 4 dynamic structure factors.")
    print("# For degenerate bands, dynamic structure factors are summed.")
    print("")

    # Use as iterator
    for Q, d, f, S in zip(Qpoints, distances, dsf.frequencies,
                          dsf.dynamic_structure_factors):
        bi_sets = degenerate_sets(f)  # to treat for band degeneracy
        text = "%f  " % d
        text += "%f %f %f  " % tuple(Q)
        text += " ".join(["%f" % (f[bi].sum() * THzToEv * 1000 / len(bi))
                          for bi in bi_sets])
        text += "  "
        text += " ".join(["%f" % (S[bi].sum()) for bi in bi_sets])
        print(text)


if __name__ == '__main__':
    phonon = load(supercell_matrix=[2, 2, 2],
                  primitive_matrix='F',
                  unitcell_filename="POSCAR")

    # Q-points in reduced coordinates wrt cubic reciprocal space
    Qpoints = [[2.970000, -2.970000, 2.970000],
               [2.950000, 2.950000, -2.950000],
               [2.930000, -2.930000, 2.930000],
               [2.905000, -2.905000, 2.905000],
               [2.895000, -2.895000, 2.895000],
               [2.880000, -2.880000, 2.880000],
               [2.850000, -2.850000, 2.850000],
               [2.810000, -2.810000, 2.810000],
               [2.735000, -2.735000, 2.735000],
               [2.660000, -2.660000, 2.660000],
               [2.580000, -2.580000, 2.580000],
               [2.500000, -2.500000, 2.500000]]

    # Mesh sampling phonon calculation is needed for Debye-Waller factor.
    # This must be done with is_mesh_symmetry=False and with_eigenvectors=True.
    mesh = [11, 11, 11]
    phonon.run_mesh(mesh,
                    is_mesh_symmetry=False,
                    with_eigenvectors=True)
    temperature = 300

    IXS = True

    if IXS:
        # For IXS, atomic form factor is needed and given as a function as
        # a parameter.
        # D. Waasmaier and A. Kirfel, Acta Cryst. A51, 416 (1995)
        # f(Q) = \sum_i a_i \exp((-b_i Q^2) + c
        # Q is in angstron^-1
        # a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, c
        f_params = {'Na': [3.148690, 2.594987, 4.073989, 6.046925,
                           0.767888, 0.070139, 0.995612, 14.1226457,
                           0.968249, 0.217037, 0.045300],  # 1+
                    'Cl': [1.061802, 0.144727, 7.139886, 1.171795,
                           6.524271, 19.467656, 2.355626, 60.320301,
                           35.829404, 0.000436, -34.916604]}  # 1-
        AFF_func = get_AFF_func(f_params)
        run(phonon,
            Qpoints,
            temperature,
            atomic_form_factor_func=AFF_func)
    else:
        # For INS, scattering length has to be given.
        # The following values is obtained at (Coh b)
        # https://www.nist.gov/ncnr/neutron-scattering-lengths-list
        run(phonon,
            Qpoints,
            temperature,
            scattering_lengths={'Na': 3.63, 'Cl': 9.5770})
