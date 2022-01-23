#include <ff/PolyhedralTopology.h>
#include <ff/PolyhedralComponents.h>
#include <ff/Polyhedron.h>
#include <iostream>

const ff::PolyhedralTopology tetrahedron_topology = {
    {{{2, 1, 0}, false}, {{0, 1, 3}, false}, {{1, 2, 3}, false}, {{2, 0, 3}, false}}, false};

std::vector<R3> tetrahedron_vertices(const double edge)
{
    const double as = edge / 2;
    const double ac = edge / sqrt(3) / 2;
    const double ah = edge / sqrt(3);
    const double height = sqrt(2. / 3) * edge;
    const double zcom = height / 4; // center of mass

    return {
        {-ac, as, -zcom},
        {-ac, -as, -zcom},
        {ah, 0., -zcom},
        {0, 0., height - zcom}};
}

int main() {
    ff::Polyhedron tetrahedron(tetrahedron_topology, tetrahedron_vertices(1.));
    for (double t=0.2; t<200;  t *= 1.01) {
        C3 q{t, 1.2*t, .3*t};
        std::cout << t << " " << std::abs(tetrahedron.formfactor(q)) << std::endl;
    }
}
