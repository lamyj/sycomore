#include <pybind11/pybind11.h>

void wrap_como_Model(pybind11::module &);

void wrap_como(pybind11::module & m)
{
    auto como = m.def_submodule("como");

    wrap_como_Model(como);
}
