#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "sycomore/epg/Model.h"

#include "../type_casters.h"

void wrap_epg_Model(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;
    
    auto model = class_<Model>(
        m, "Model",
        "Model a single- or two-pool system expressed in EPG formalism");
    enum_<Model::Kind>(model, "Kind")
        .value("SinglePool", Model::SinglePool)
        .value("Exchange", Model::Exchange)
        .value("MagnetizationTransfer", Model::MagnetizationTransfer)
        .export_values();
}
