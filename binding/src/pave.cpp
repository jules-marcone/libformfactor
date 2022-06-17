//  ************************************************************************************************
//
//! @file      demo/pave_model.cpp
//! @brief     Computes the scattering intensity for a given scattering vector
//!
//! @homepage  https://github.com/jules-marcone/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Jules Marcone, June 2022
//
//  ************************************************************************************************

#include <ff/Cuboid.h>
#include <binding/pave.h>
#include <iostream>

constexpr double twopi = 6.28318530718;

//! Prints list t vs |F(q(t))| for a logarithmic range of t values
double Iqabc(double qa, double qb, double qc,
    double sld,
    double solvent_sld,
    double edge_a,
    double edge_b,
    double edge_c)
{
    ff::cuboid::Pave pave(edge_a, edge_b, edge_c);
    C3 q(qa, qb, qc);

    // Amplitude AP from eqn. (13)
    const double AP = abs(pave.formfactor(q));

    // Multiply by contrast and volume
    const double s = (sld - solvent_sld) * (edge_a * edge_b * edge_c);

    // Convert from [1e-12 A-1] to [cm-1]
    return 1.0e-4 * (s * AP)*(s * AP);
}


