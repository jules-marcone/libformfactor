#include "ff/Factorial.h"
#include "catch.hpp"
#include <memory>

namespace {
constexpr auto ReciprocalFactorialArray = ff_aux::generateReciprocalFactorialArray<171>();

} // namespace


TEST_CASE("FactorialTest", "")
{
    CHECK(ReciprocalFactorialArray.size() > 150);
    CHECK(ReciprocalFactorialArray[0] == Approx(1.));
    CHECK(ReciprocalFactorialArray[1] == Approx(1.));
    CHECK(ReciprocalFactorialArray[2] == Approx(0.5));
    CHECK(ReciprocalFactorialArray[3] == Approx(1.0 / 6));
    /* the following disabled because tgamma is too unprecise under
       old versions of glibc (at leat up to 2.12, but less than 2.22)
    for( size_t k=4; k<precomputed.factorial.size(); ++k )
        EXPECT_NEAR(precomputed.factorial[k], tgamma(k+1.), 12*eps*tgamma(k+1.) );
    */
    CHECK(ReciprocalFactorialArray[150] == Approx(1.75027620692601519e-263).epsilon(1e-14));
}
