#include <pybind11/pybind11.h>

#include "../type_casters.h"

void wrap_isochromat_Model(pybind11::module &);
void wrap_isochromat_Operator(pybind11::module &);

void wrap_isochromat(pybind11::module & m)
{
    auto isochromat = m.def_submodule("isochromat");
    
    wrap_isochromat_Model(isochromat);
    wrap_isochromat_Operator(isochromat);
}
