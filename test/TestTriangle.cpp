#include "ff/PolyhedralComponents.h"
#include "catch.hpp"
#include <cstdio>

#undef M_PI_2
#define M_PI_2 1.57079632679489661923     /* pi/2 */

//! Ad-hoc test of triangle form factor.
//!
//! Used while preparing polyhedral form factor manuscript for publication - JWu, dec 2020.

TEST_CASE("FF:Triangle", "")
{
    const double a = 1.;
    const double as = a / 2;
    const double ac = a / sqrt(3) / 2;
    const double ah = a / sqrt(3);
    const std::vector<R3> V{{-ac, as, 0.}, {-ac, -as, 0.}, {ah, 0., 0.}};

    const PolyhedralFace T(V, false);
    CHECK(std::abs(sqrt(3) / 4-T.area()) < 1e-15);

    int failures = 0;
    const int M = 37;
    for (int j = 0; j < M; ++j) {
        const double phi = M_PI_2 * j / (M - 1);
        const C3 uQ{sin(phi), cos(phi), 0.};
        const int N = 2800 + j;
        for (int i = 0; i < N; ++i) {
            const double q = 1e-17 * pow(1.7e17, i / (N - 1.));
            const C3 Q = uQ * q;
            const double f1 = std::abs(T.ff_2D_direct(Q));
            const double f2 = std::abs(T.ff_2D_expanded(Q));
            const double relerr = std::abs(f1 - f2) / f2;
            if (relerr > 7e-16) {
                printf("ERR1 %9.6f %16.11e %21.16e %21.16e %10.4e\n", phi, q, f1, f2, relerr);
                ++failures;
            }
            if (q > 1e-7)
                continue;
            const double relerr2 = std::abs(f1 - T.area()) / f2;
            if (relerr2 > 7e-16) {
                printf("ERR2 %9.6f %16.11e %21.16e %21.16e %10.4e\n", phi, q, f1, f2, relerr2);
                ++failures;
            }
        }
    }
    CHECK(0 == failures);
}
