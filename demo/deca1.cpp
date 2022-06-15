//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      demobipy/bipy3.cpp
//! @brief     regular bipyramid with a triangular base
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Marianne Imperor-Clerc May 2022
//
//  ************************************************************************************************

#include <ff/Penta.h>
//#include <ff/Bipy.h>
#include <iostream>

constexpr double twopi = 6.28318530718;

//! Prints list t vs |F(q(t))| for a logarithmic range of t values

int main()
{
    std::cout << "# regular decahedron \n";
    ff::penta::Decahedron p1(1.); 
    ff::penta::Decahedron p2(2.); 
    ff::penta::Decahedron p3(3.);

//    ff::bipy::Bipy3 p1(1.); 
//    ff::bipy::Bipy3 p2(2.); 
//    ff::bipy::Bipy3 p3(3.); 

    std::cout << "# edge=1, 2 and 3 regular decahedron \n";
    std::cout << " volume = " << p1.volume() << " " << p2.volume() << " " << p3.volume()  << std::endl;
    std::cout << " radius = " << p1.radius() << " " << p2.radius() << " " << p3.radius()  << std::endl;

    std::cout << "# q vs |F(q)| for q in direction 111 \n";
    for (double t = 0.2; t < 1.; t *= 1.1) {
        C3 q(t/sqrt(3), t/sqrt(3), t/sqrt(3));
        std::cout << t << " " << std::abs(p1.formfactor(q)) << " " << std::abs(p2.formfactor(q)) << " " << std::abs(p3.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;   

    std::cout << "# q vs |F(q)| for q in direction 100 \n";
    for (double t = 0.2; t < 1.; t *= 1.1) {
        C3 q(t, 0., 0.);
        std::cout << t << " " << std::abs(p1.formfactor(q)) << " " << std::abs(p2.formfactor(q)) << " " << std::abs(p3.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;   

}
