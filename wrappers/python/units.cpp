#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/units.h"

void wrap_units(pybind11::module & module)
{
    using namespace pybind11;
    using namespace sycomore::units;

    auto units = module.def_submodule("units");

#define SYCOMORE_WRAPPED_UNIT(name) units.attr(#name) = name;
#define SYCOMORE_WRAPPED_UNITS(name) \
    SYCOMORE_WRAPPED_UNIT(Y##name) \
    SYCOMORE_WRAPPED_UNIT(Z##name) \
    SYCOMORE_WRAPPED_UNIT(E##name) \
    SYCOMORE_WRAPPED_UNIT(P##name) \
    SYCOMORE_WRAPPED_UNIT(T##name) \
    SYCOMORE_WRAPPED_UNIT(G##name) \
    SYCOMORE_WRAPPED_UNIT(M##name) \
    SYCOMORE_WRAPPED_UNIT(k##name) \
    SYCOMORE_WRAPPED_UNIT(h##name) \
    SYCOMORE_WRAPPED_UNIT(da##name) \
    SYCOMORE_WRAPPED_UNIT(name) \
    SYCOMORE_WRAPPED_UNIT(d##name) \
    SYCOMORE_WRAPPED_UNIT(c##name) \
    SYCOMORE_WRAPPED_UNIT(m##name) \
    SYCOMORE_WRAPPED_UNIT(u##name) \
    SYCOMORE_WRAPPED_UNIT(n##name) \
    SYCOMORE_WRAPPED_UNIT(p##name) \
    SYCOMORE_WRAPPED_UNIT(f##name) \
    SYCOMORE_WRAPPED_UNIT(a##name) \
    SYCOMORE_WRAPPED_UNIT(z##name) \
    SYCOMORE_WRAPPED_UNIT(y##name)

    SYCOMORE_WRAPPED_UNITS(m);

    SYCOMORE_WRAPPED_UNITS(g);

    SYCOMORE_WRAPPED_UNITS(s);
    SYCOMORE_WRAPPED_UNIT(min);
    SYCOMORE_WRAPPED_UNIT(h);

    SYCOMORE_WRAPPED_UNITS(A);
    SYCOMORE_WRAPPED_UNITS(K);
    SYCOMORE_WRAPPED_UNITS(mol);
    SYCOMORE_WRAPPED_UNITS(cd);

    SYCOMORE_WRAPPED_UNIT(deg);
    SYCOMORE_WRAPPED_UNIT(rad);

    SYCOMORE_WRAPPED_UNIT(sr);

    SYCOMORE_WRAPPED_UNITS(Hz);
    SYCOMORE_WRAPPED_UNITS(N);
    SYCOMORE_WRAPPED_UNITS(Pa);
    SYCOMORE_WRAPPED_UNITS(J);
    SYCOMORE_WRAPPED_UNITS(W);
    SYCOMORE_WRAPPED_UNITS(C);
    SYCOMORE_WRAPPED_UNITS(V);
    SYCOMORE_WRAPPED_UNITS(F);
    SYCOMORE_WRAPPED_UNITS(Ohm);
    SYCOMORE_WRAPPED_UNITS(S);
    SYCOMORE_WRAPPED_UNITS(Wb);
    SYCOMORE_WRAPPED_UNITS(T);
    SYCOMORE_WRAPPED_UNITS(H);
    SYCOMORE_WRAPPED_UNITS(lm);
    SYCOMORE_WRAPPED_UNITS(lx);
    SYCOMORE_WRAPPED_UNITS(Bq);
    SYCOMORE_WRAPPED_UNITS(Gy);
    SYCOMORE_WRAPPED_UNITS(Sv);
    SYCOMORE_WRAPPED_UNITS(kat);
}
