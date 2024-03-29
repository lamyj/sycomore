#include <pybind11/pybind11.h>

#include "../type_casters.h"

void wrap_epg_Base(pybind11::module &);
void wrap_epg_Discrete(pybind11::module &);
void wrap_epg_Discrete3D(pybind11::module &);
void wrap_epg_Model(pybind11::module &);
void wrap_epg_operators(pybind11::module &);
void wrap_epg_Regular(pybind11::module &);

void wrap_epg(pybind11::module & m)
{
    auto epg = m.def_submodule("epg");
    
    wrap_epg_Model(epg);
    wrap_epg_Base(epg);
    wrap_epg_Discrete(epg);
    wrap_epg_Discrete3D(epg);
    wrap_epg_operators(epg);
    wrap_epg_Regular(epg);
}
